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
   Contributing author: Ray Shan (SNL)
------------------------------------------------------------------------- */

#include "pair_coul_long_cs_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

// the core/shell Ewald correction uses a longer erfc series than the A1..A5
// form of the parent style, and a minimal separation so that r = 0
// core/shell pairs stay finite until the special-bond factor removes them;
// both are taken verbatim from the CPU style

// the B series goes with its own EWALD_P, not the value EwaldConst pairs
// with A1..A5

static constexpr double EWALD_P_CS = 9.95473818e-1;

static constexpr double B0 = -0.1335096380159268;
static constexpr double B1 = -2.57839507e-1;
static constexpr double B2 = -1.37203639e-1;
static constexpr double B3 = -8.88822059e-3;
static constexpr double B4 = -5.80844129e-3;
static constexpr double B5 =  1.14652755e-1;

static constexpr double EPSILON = 1.0e-20;
static constexpr double EPS_EWALD = 1.0e-6;
static constexpr double EPS_EWALD_SQR = 1.0e-12;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCoulLongCSKokkos<DeviceType>::PairCoulLongCSKokkos(LAMMPS *lmp):PairCoulLongCS(lmp)
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
PairCoulLongCSKokkos<DeviceType>::~PairCoulLongCSKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongCSKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

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
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  newton_pair = force->newton_pair;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  if (ncoultablebits)
    ev = pair_compute<PairCoulLongCSKokkos<DeviceType>,CoulLongTable<1> >
      (this,(NeighListKokkos<DeviceType>*)list);
  else
    ev = pair_compute<PairCoulLongCSKokkos<DeviceType>,CoulLongTable<0> >
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
   compute coulomb pair force between atoms i and j
   ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS,  class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulLongCSKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT g_ewald_kk = static_cast<KK_FLOAT>(g_ewald);
  const KK_FLOAT tabinnersq_kk = static_cast<KK_FLOAT>(tabinnersq);

  // r = 0 must stay finite here; the special-bond factor removes the pair

  const KK_FLOAT rsq_cs = rsq + static_cast<KK_FLOAT>(EPSILON);

  if (Specialisation::DoTable && rsq_cs > tabinnersq_kk) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq_cs;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction = ((KK_FLOAT)rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_ftable[itable] + fraction*d_dftable[itable];
    KK_FLOAT forcecoul = qtmp*q[j] * table;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) {
      const KK_FLOAT ctable = d_ctable[itable] + fraction*d_dctable[itable];
      const KK_FLOAT prefactor = qtmp*q[j] * ctable;
      forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
    }
    return forcecoul/rsq_cs;
  } else {
    const KK_FLOAT r = Kokkos::sqrt(rsq_cs);

    if (factor_coul < static_cast<KK_FLOAT>(1.0)) {

      // a bonded core/shell pair needs the minimal separation EPS_EWALD added
      // to both the prefactor and the Ewald argument to keep the approximation
      // valid, and the 1/r^2 scaling adjusted to match

      const KK_FLOAT reps = r + static_cast<KK_FLOAT>(EPS_EWALD);
      const KK_FLOAT grij = g_ewald_kk * reps;
      const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
      const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
        (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P_CS)*grij);
      const KK_FLOAT u = static_cast<KK_FLOAT>(1.0) - t;
      const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(1.0)+u*(static_cast<KK_FLOAT>(B0)+
        u*(static_cast<KK_FLOAT>(B1)+u*(static_cast<KK_FLOAT>(B2)+u*(static_cast<KK_FLOAT>(B3)+
        u*(static_cast<KK_FLOAT>(B4)+u*static_cast<KK_FLOAT>(B5))))))) * expm2;
      const KK_FLOAT prefactor = qqrd2e * qtmp*q[j] / reps;
      const KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2 -
                                              (static_cast<KK_FLOAT>(1.0)-factor_coul));
      return forcecoul / (rsq_cs + static_cast<KK_FLOAT>(EPS_EWALD_SQR));

    } else {

      const KK_FLOAT grij = g_ewald_kk * r;
      const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
      const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
        (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P_CS)*grij);
      const KK_FLOAT u = static_cast<KK_FLOAT>(1.0) - t;
      const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(1.0)+u*(static_cast<KK_FLOAT>(B0)+
        u*(static_cast<KK_FLOAT>(B1)+u*(static_cast<KK_FLOAT>(B2)+u*(static_cast<KK_FLOAT>(B3)+
        u*(static_cast<KK_FLOAT>(B4)+u*static_cast<KK_FLOAT>(B5))))))) * expm2;
      const KK_FLOAT prefactor = qqrd2e * qtmp*q[j] / r;
      const KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2);
      return forcecoul / rsq_cs;
    }
  }
}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulLongCSKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT g_ewald_kk = static_cast<KK_FLOAT>(g_ewald);
  const KK_FLOAT tabinnersq_kk = static_cast<KK_FLOAT>(tabinnersq);

  const KK_FLOAT rsq_cs = rsq + static_cast<KK_FLOAT>(EPSILON);

  if (Specialisation::DoTable && rsq_cs > tabinnersq_kk) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq_cs;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction = ((KK_FLOAT)rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_etable[itable] + fraction*d_detable[itable];
    KK_FLOAT ecoul = qtmp*q[j] * table;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) {
      const KK_FLOAT ctable = d_ctable[itable] + fraction*d_dctable[itable];
      const KK_FLOAT prefactor = qtmp*q[j] * ctable;
      ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
    }
    return ecoul;
  } else {
    const KK_FLOAT r = Kokkos::sqrt(rsq_cs);
    const KK_FLOAT reps = (factor_coul < static_cast<KK_FLOAT>(1.0)) ?
      r + static_cast<KK_FLOAT>(EPS_EWALD) : r;
    {
      const KK_FLOAT grij = g_ewald_kk * reps;
      const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
      const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
        (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P_CS)*grij);
      const KK_FLOAT u = static_cast<KK_FLOAT>(1.0) - t;
      const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(1.0)+u*(static_cast<KK_FLOAT>(B0)+
        u*(static_cast<KK_FLOAT>(B1)+u*(static_cast<KK_FLOAT>(B2)+u*(static_cast<KK_FLOAT>(B3)+
        u*(static_cast<KK_FLOAT>(B4)+u*static_cast<KK_FLOAT>(B5))))))) * expm2;
      const KK_FLOAT prefactor = qqrd2e * qtmp*q[j] / reps;
      KK_FLOAT ecoul = prefactor * erfc;
      if (factor_coul < static_cast<KK_FLOAT>(1.0))
        ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
      return ecoul;
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongCSKokkos<DeviceType>::allocate()
{
  PairCoulLongCS::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  d_cut_ljsq = typename AT::t_kkfloat_2d("pair:cut_ljsq",n+1,n+1);
  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType>("PairCoulLongCS::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

template<class DeviceType>
void PairCoulLongCSKokkos<DeviceType>::init_tables(double cut_coul, double *cut_respa)
{
  Pair::init_tables(cut_coul,cut_respa);

  typedef typename AT::t_kkfloat_1d table_type;
  typedef HAT::t_kkfloat_1d host_table_type;

  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;


  // Copy rtable and drtable
  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) {
    h_table(i) = static_cast<KK_FLOAT>(rtable[i]);
  }
  Kokkos::deep_copy(d_table,h_table);
  d_rtable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);
  for (int i = 0; i < ntable; i++) {
    h_table(i) = static_cast<KK_FLOAT>(drtable[i]);
  }
  Kokkos::deep_copy(d_table,h_table);
  d_drtable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy ftable and dftable
  for (int i = 0; i < ntable; i++) {
    h_table(i) = static_cast<KK_FLOAT>(ftable[i]);
  }
  Kokkos::deep_copy(d_table,h_table);
  d_ftable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for (int i = 0; i < ntable; i++) {
    h_table(i) = static_cast<KK_FLOAT>(dftable[i]);
  }
  Kokkos::deep_copy(d_table,h_table);
  d_dftable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy ctable and dctable
  for (int i = 0; i < ntable; i++) {
    h_table(i) = static_cast<KK_FLOAT>(ctable[i]);
  }
  Kokkos::deep_copy(d_table,h_table);
  d_ctable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for (int i = 0; i < ntable; i++) {
    h_table(i) = static_cast<KK_FLOAT>(dctable[i]);
  }
  Kokkos::deep_copy(d_table,h_table);
  d_dctable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  // Copy etable and detable
  for (int i = 0; i < ntable; i++) {
    h_table(i) = static_cast<KK_FLOAT>(etable[i]);
  }
  Kokkos::deep_copy(d_table,h_table);
  d_etable = d_table;
  }

  {
  host_table_type h_table("HostTable",ntable);
  table_type d_table("DeviceTable",ntable);

  for (int i = 0; i < ntable; i++) {
    h_table(i) = static_cast<KK_FLOAT>(detable[i]);
  }
  Kokkos::deep_copy(d_table,h_table);
  d_detable = d_table;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongCSKokkos<DeviceType>::init_style()
{
  PairCoulLongCS::init_style();

  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(cut_coulsq));
  Kokkos::deep_copy(d_cut_ljsq,static_cast<KK_FLOAT>(cut_coulsq));

  // error if rRESPA with inner levels

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  // adjust neighbor list request for KOKKOS

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
double PairCoulLongCSKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairCoulLongCS::init_one(i,j);

  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq);

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(cut_coulsq);
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cut_coulsq);
  }

  k_cutsq.view_host()(i,j) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairCoulLongCSKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairCoulLongCSKokkos<LMPHostType>;
#endif
}

