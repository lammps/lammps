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

/* ----------------------------------------------------------------------
   Contributing author: Claude Opus 5.5 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_coul_esp_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCoulEspKokkos<DeviceType>::PairCoulEspKokkos(LAMMPS *lmp):PairCoulEsp(lmp)
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
PairCoulEspKokkos<DeviceType>::~PairCoulEspKokkos()
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
void PairCoulEspKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  cut_coul_kk = static_cast<KK_FLOAT>(cut_coul);

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  if (ncoultablebits)
    ev = pair_compute<PairCoulEspKokkos<DeviceType>,CoulLongTable<1> >
      (this,(NeighListKokkos<DeviceType>*)list);
  else
    ev = pair_compute<PairCoulEspKokkos<DeviceType>,CoulLongTable<0> >
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
KK_FLOAT PairCoulEspKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& itype, const int& jtype,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT tabinnersq_kk = static_cast<KK_FLOAT>(tabinnersq);
  const KK_FLOAT scale_kk = (STACKPARAMS && itype<MAX_TYPES_STACKPARAMS+1 && jtype<MAX_TYPES_STACKPARAMS+1) ?
    m_params[itype][jtype].scale : params(itype,jtype).scale;
  if (Specialisation::DoTable && rsq > tabinnersq_kk) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction = ((KK_FLOAT)rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_ftable[itable] + fraction*d_dftable[itable];
    KK_FLOAT forcecoul = scale_kk*qtmp*q[j] * table;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) {
      const KK_FLOAT table = d_ctable[itable] + fraction*d_dctable[itable];
      const KK_FLOAT prefactor = scale_kk*qtmp*q[j] * table;
      forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
    }
    return forcecoul/rsq;
  } else {
    // polynomial approximation of the ESP short-range kernel in the scaled
    // distance 2r/rc - 1, with the coefficients provided by kspace style esp

    const KK_FLOAT r = Kokkos::sqrt(rsq);
    const KK_FLOAT r_scal = static_cast<KK_FLOAT>(2.0) * r / cut_coul_kk - static_cast<KK_FLOAT>(1.0);
    KK_FLOAT force_poly_appx = d_force_poly_coeff(0);
    KK_FLOAT force_poly_r = static_cast<KK_FLOAT>(1.0);
    for (int index = 1; index < num_of_force_poly; index++) {
      force_poly_r *= r_scal;
      force_poly_appx += d_force_poly_coeff(index) * force_poly_r;
    }
    const KK_FLOAT prefactor = qqrd2e * scale_kk * qtmp*q[j]/r;
    KK_FLOAT forcecoul = prefactor * force_poly_appx;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;

    return forcecoul/rsq;
  }
}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulEspKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& itype, const int& jtype,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT tabinnersq_kk = static_cast<KK_FLOAT>(tabinnersq);
  const KK_FLOAT scale_kk = (STACKPARAMS && itype<MAX_TYPES_STACKPARAMS+1 && jtype<MAX_TYPES_STACKPARAMS+1) ?
    m_params[itype][jtype].scale : params(itype,jtype).scale;
  if (Specialisation::DoTable && rsq > tabinnersq_kk) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;
    const KK_FLOAT fraction = ((KK_FLOAT)rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_etable[itable] + fraction*d_detable[itable];
    KK_FLOAT ecoul = scale_kk*qtmp*q[j] * table;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) {
      const KK_FLOAT table = d_ctable[itable] + fraction*d_dctable[itable];
      const KK_FLOAT prefactor = scale_kk*qtmp*q[j] * table;
      ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
    }
    return ecoul;
  } else {
    const KK_FLOAT r = Kokkos::sqrt(rsq);
    const KK_FLOAT r_scal = static_cast<KK_FLOAT>(2.0) * r / cut_coul_kk - static_cast<KK_FLOAT>(1.0);
    KK_FLOAT energy_poly_appx = d_energy_poly_coeff(0);
    KK_FLOAT energy_poly_r = static_cast<KK_FLOAT>(1.0);
    for (int index = 1; index < num_of_energy_poly; index++) {
      energy_poly_r *= r_scal;
      energy_poly_appx += d_energy_poly_coeff(index) * energy_poly_r;
    }
    const KK_FLOAT prefactor = qqrd2e * scale_kk * qtmp*q[j]/r;
    KK_FLOAT ecoul = prefactor * energy_poly_appx;
    if (factor_coul < static_cast<KK_FLOAT>(1.0)) ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;
    return ecoul;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulEspKokkos<DeviceType>::allocate()
{
  PairCoulEsp::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  d_cut_ljsq = typename AT::t_kkfloat_2d("pair:cut_ljsq",n+1,n+1);
  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);

  k_params = Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType>("PairCoulEsp::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

template<class DeviceType>
void PairCoulEspKokkos<DeviceType>::init_tables(double cut_coul, double *cut_respa)
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
void PairCoulEspKokkos<DeviceType>::init_style()
{
  PairCoulEsp::init_style();

  // copy the polynomial coefficients of the short-range kernel, which the
  // base class just took from kspace style esp, to the device

  d_force_poly_coeff = typename AT::t_kkfloat_1d("pair:force_poly_coeff", num_of_force_poly);
  auto h_force_poly = Kokkos::create_mirror_view(d_force_poly_coeff);
  for (int i = 0; i < num_of_force_poly; ++i)
    h_force_poly(i) = static_cast<KK_FLOAT>(force_poly_coeff[i]);
  Kokkos::deep_copy(d_force_poly_coeff, h_force_poly);

  d_energy_poly_coeff = typename AT::t_kkfloat_1d("pair:energy_poly_coeff", num_of_energy_poly);
  auto h_energy_poly = Kokkos::create_mirror_view(d_energy_poly_coeff);
  for (int i = 0; i < num_of_energy_poly; ++i)
    h_energy_poly(i) = static_cast<KK_FLOAT>(energy_poly_coeff[i]);
  Kokkos::deep_copy(d_energy_poly_coeff, h_energy_poly);

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
double PairCoulEspKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairCoulEsp::init_one(i,j);

  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq);
  k_params.view_host()(i,j).scale = static_cast<KK_FLOAT>(scale[i][j]);

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
template class PairCoulEspKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairCoulEspKokkos<LMPHostType>;
#endif
}

