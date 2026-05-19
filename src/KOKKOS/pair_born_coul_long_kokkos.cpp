// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "pair_born_coul_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairBornCoulLongKokkos<DeviceType>::PairBornCoulLongKokkos(LAMMPS *lmp) : PairBornCoulLong(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairBornCoulLongKokkos<DeviceType>::~PairBornCoulLongKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
    memoryKK->destroy_kokkos(k_cut_ljsq,cut_ljsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairBornCoulLongKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_cut_ljsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  newton_pair = force->newton_pair;

  g_ewald_kk = static_cast<KK_FLOAT>(g_ewald);

  copymode = 1;
  EV_FLOAT ev;
  if (ncoultablebits)
    ev = pair_compute<PairBornCoulLongKokkos<DeviceType>,CoulLongTable<1>>
      (this,(NeighListKokkos<DeviceType>*)list);
  else
    ev = pair_compute<PairBornCoulLongKokkos<DeviceType>,CoulLongTable<0>>
      (this,(NeighListKokkos<DeviceType>*)list);

  if (eflag) {
    eng_vdwl += static_cast<double>(ev.evdwl);
    eng_coul += static_cast<double>(ev.ecoul);
  }

  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   compute Born pair force between atoms i and j
------------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBornCoulLongKokkos<DeviceType>::
compute_fpair(const KK_FLOAT &rsq, const int & /*i*/, const int & /*j*/,
              const int &itype, const int &jtype) const {
  const KK_FLOAT cut_ljsq = STACKPARAMS ? m_params[itype][jtype].cut_ljsq : params(itype,jtype).cut_ljsq;
  if (rsq >= cut_ljsq) return static_cast<KK_FLOAT>(0.0);
  const KK_FLOAT r2inv  = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv  = r2inv*r2inv*r2inv;
  const KK_FLOAT r      = sqrt(rsq);
  const KK_FLOAT rhoinv = STACKPARAMS ? m_params[itype][jtype].rhoinv : params(itype,jtype).rhoinv;
  const KK_FLOAT sigma  = STACKPARAMS ? m_params[itype][jtype].sigma  : params(itype,jtype).sigma;
  const KK_FLOAT born1  = STACKPARAMS ? m_params[itype][jtype].born1  : params(itype,jtype).born1;
  const KK_FLOAT born2  = STACKPARAMS ? m_params[itype][jtype].born2  : params(itype,jtype).born2;
  const KK_FLOAT born3  = STACKPARAMS ? m_params[itype][jtype].born3  : params(itype,jtype).born3;
  const KK_FLOAT rexp   = exp((sigma - r) * rhoinv);
  const KK_FLOAT forceborn = born1*r*rexp - born2*r6inv + born3*r2inv*r6inv;
  return forceborn*r2inv;
}

/* ----------------------------------------------------------------------
   compute Coulomb pair force between atoms i and j
------------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBornCoulLongKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT &rsq, const int & /*i*/, const int &j,
              const int & /*itype*/, const int & /*jtype*/,
              const KK_FLOAT &factor_coul, const KK_FLOAT &qtmp) const {
  if (Specialisation::DoTable && rsq > tabinnersq_kk) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction = ((KK_FLOAT)rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_ftable[itable] + fraction*d_dftable[itable];
    KK_FLOAT forcecoul = qtmp*q[j] * table;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) {
      const KK_FLOAT table2 = d_ctable[itable] + fraction*d_dctable[itable];
      const KK_FLOAT prefactor = qtmp*q[j] * table2;
      forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
    }
    return forcecoul/rsq;
  } else {
    const KK_FLOAT r      = sqrt(rsq);
    const KK_FLOAT grij   = g_ewald_kk * r;
    const KK_FLOAT expm2  = exp(-grij*grij);
    const KK_FLOAT t      = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
    const KK_FLOAT erfc   = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+t*(static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+t*static_cast<KK_FLOAT>(A5))))) * expm2;
    const KK_FLOAT prefactor = qqrd2e * qtmp*q[j]/r;
    KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2);
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
    return forcecoul/rsq;
  }
}

/* ----------------------------------------------------------------------
   compute Born pair energy between atoms i and j
------------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBornCoulLongKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int & /*i*/, const int & /*j*/,
              const int &itype, const int &jtype) const {
  const KK_FLOAT cut_ljsq = STACKPARAMS ? m_params[itype][jtype].cut_ljsq : params(itype,jtype).cut_ljsq;
  if (rsq >= cut_ljsq) return static_cast<KK_FLOAT>(0.0);
  const KK_FLOAT r2inv  = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv  = r2inv*r2inv*r2inv;
  const KK_FLOAT r      = sqrt(rsq);
  const KK_FLOAT rhoinv = STACKPARAMS ? m_params[itype][jtype].rhoinv : params(itype,jtype).rhoinv;
  const KK_FLOAT sigma  = STACKPARAMS ? m_params[itype][jtype].sigma  : params(itype,jtype).sigma;
  const KK_FLOAT a      = STACKPARAMS ? m_params[itype][jtype].a      : params(itype,jtype).a;
  const KK_FLOAT born2  = STACKPARAMS ? m_params[itype][jtype].born2  : params(itype,jtype).born2;
  const KK_FLOAT born3  = STACKPARAMS ? m_params[itype][jtype].born3  : params(itype,jtype).born3;
  const KK_FLOAT offset = STACKPARAMS ? m_params[itype][jtype].offset : params(itype,jtype).offset;
  const KK_FLOAT rexp   = exp((sigma - r) * rhoinv);
  // born2 = 6*c, born3 = 8*d
  return a*rexp - (born2/static_cast<KK_FLOAT>(6.0))*r6inv
         + (born3/static_cast<KK_FLOAT>(8.0))*r6inv*r2inv - offset;
}

/* ----------------------------------------------------------------------
   compute Coulomb pair energy between atoms i and j
------------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairBornCoulLongKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT &rsq, const int & /*i*/, const int &j,
              const int & /*itype*/, const int & /*jtype*/,
              const KK_FLOAT &factor_coul, const KK_FLOAT &qtmp) const {
  if (Specialisation::DoTable && rsq > tabinnersq_kk) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction = ((KK_FLOAT)rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_etable[itable] + fraction*d_detable[itable];
    KK_FLOAT ecoul = qtmp*q[j] * table;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) {
      const KK_FLOAT table2 = d_ctable[itable] + fraction*d_dctable[itable];
      const KK_FLOAT prefactor = qtmp*q[j] * table2;
      ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
    }
    return ecoul;
  } else {
    const KK_FLOAT r      = sqrt(rsq);
    const KK_FLOAT grij   = g_ewald_kk * r;
    const KK_FLOAT expm2  = exp(-grij*grij);
    const KK_FLOAT t      = static_cast<KK_FLOAT>(1.0) / (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
    const KK_FLOAT erfc   = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+t*(static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+t*static_cast<KK_FLOAT>(A5))))) * expm2;
    const KK_FLOAT prefactor = qqrd2e * qtmp*q[j]/r;
    KK_FLOAT ecoul = prefactor*erfc;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
    return ecoul;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBornCoulLongKokkos<DeviceType>::allocate()
{
  PairBornCoulLong::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  memory->destroy(cut_ljsq);
  memoryKK->create_kokkos(k_cut_ljsq,cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();

  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_born_coul_long**,Kokkos::LayoutRight,DeviceType>("PairBornCoulLong::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

template<class DeviceType>
void PairBornCoulLongKokkos<DeviceType>::init_tables(double cut_coul, double *cut_respa)
{
  Pair::init_tables(cut_coul,cut_respa);

  typedef typename AT::t_kkfloat_1d table_type;
  typedef HAT::t_kkfloat_1d host_table_type;

  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;

  tabinnersq_kk = static_cast<KK_FLOAT>(tabinnersq);

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(rtable[i]);
  Kokkos::deep_copy(d_table,h_table);
  d_rtable = d_table;
  }
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(drtable[i]);
  Kokkos::deep_copy(d_table,h_table);
  d_drtable = d_table;
  }
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(ftable[i]);
  Kokkos::deep_copy(d_table,h_table);
  d_ftable = d_table;
  }
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(dftable[i]);
  Kokkos::deep_copy(d_table,h_table);
  d_dftable = d_table;
  }
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(ctable[i]);
  Kokkos::deep_copy(d_table,h_table);
  d_ctable = d_table;
  }
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(dctable[i]);
  Kokkos::deep_copy(d_table,h_table);
  d_dctable = d_table;
  }
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(etable[i]);
  Kokkos::deep_copy(d_table,h_table);
  d_etable = d_table;
  }
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) h_table(i) = static_cast<KK_FLOAT>(detable[i]);
  Kokkos::deep_copy(d_table,h_table);
  d_detable = d_table;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairBornCoulLongKokkos<DeviceType>::init_style()
{
  PairBornCoulLong::init_style();

  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(cut_coulsq));

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa) error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairBornCoulLongKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairBornCoulLong::init_one(i,j);
  double cut_ljsqm = cut_ljsq[i][j];

  k_params.view_host()(i,j).a        = static_cast<KK_FLOAT>(a[i][j]);
  k_params.view_host()(i,j).rhoinv   = static_cast<KK_FLOAT>(rhoinv[i][j]);
  k_params.view_host()(i,j).sigma    = static_cast<KK_FLOAT>(sigma[i][j]);
  k_params.view_host()(i,j).born1    = static_cast<KK_FLOAT>(born1[i][j]);
  k_params.view_host()(i,j).born2    = static_cast<KK_FLOAT>(born2[i][j]);
  k_params.view_host()(i,j).born3    = static_cast<KK_FLOAT>(born3[i][j]);
  k_params.view_host()(i,j).offset   = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i,j).cut_ljsq  = static_cast<KK_FLOAT>(cut_ljsqm);
  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq);
  k_params.view_host()(i,j).cutsq   = static_cast<KK_FLOAT>(cutone*cutone);

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cut_ljsqm);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(cut_coulsq);
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_cut_ljsq.view_host()(i,j) = k_cut_ljsq.view_host()(j,i) = cut_ljsqm;
  k_cut_ljsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairBornCoulLongKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairBornCoulLongKokkos<LMPHostType>;
#endif
}
