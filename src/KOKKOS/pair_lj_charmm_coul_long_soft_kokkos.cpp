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

#include "pair_lj_charmm_coul_long_soft_kokkos.h"

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

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCharmmCoulLongSoftKokkos<DeviceType>::PairLJCharmmCoulLongSoftKokkos(LAMMPS *lmp):PairLJCharmmCoulLongSoft(lmp)
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
PairLJCharmmCoulLongSoftKokkos<DeviceType>::~PairLJCharmmCoulLongSoftKokkos()
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
void PairLJCharmmCoulLongSoftKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);
  special_coul[0] = static_cast<KK_FLOAT>(force->special_coul[0]);
  special_coul[1] = static_cast<KK_FLOAT>(force->special_coul[1]);
  special_coul[2] = static_cast<KK_FLOAT>(force->special_coul[2]);
  special_coul[3] = static_cast<KK_FLOAT>(force->special_coul[3]);
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  g_ewald_kk = static_cast<KK_FLOAT>(g_ewald);
  cut_ljsq_kk = static_cast<KK_FLOAT>(cut_ljsq);
  cut_lj_innersq_kk = static_cast<KK_FLOAT>(cut_lj_innersq);
  denom_lj_inv_kk = static_cast<KK_FLOAT>(1.0 / denom_lj);
  newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<PairLJCharmmCoulLongSoftKokkos<DeviceType>,void >
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

}

/* ----------------------------------------------------------------------
   compute LJ 12-6 pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmCoulLongSoftKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  if (rsq >= (STACKPARAMS?m_params[itype][jtype].cut_ljsq:params(itype,jtype).cut_ljsq))
    return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT lj1 = (STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1);
  const KK_FLOAT lj2 = (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2);
  const KK_FLOAT lj3 = (STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3);
  const KK_FLOAT eps = (STACKPARAMS?m_params[itype][jtype].epsilon:params(itype,jtype).epsilon);

  const KK_FLOAT r4sig6 = rsq*rsq / lj2;
  const KK_FLOAT denlj = lj3 + rsq*r4sig6;

  KK_FLOAT forcelj = lj1 * eps * (static_cast<KK_FLOAT>(48.0)*r4sig6/(denlj*denlj*denlj) -
                                  static_cast<KK_FLOAT>(24.0)*r4sig6/(denlj*denlj));

  if (rsq > cut_lj_innersq_kk) {
    const KK_FLOAT switch1 = (cut_ljsq_kk-rsq) * (cut_ljsq_kk-rsq) *
      (cut_ljsq_kk + static_cast<KK_FLOAT>(2.0)*rsq -
       static_cast<KK_FLOAT>(3.0)*cut_lj_innersq_kk) * denom_lj_inv_kk;
    const KK_FLOAT switch2 = static_cast<KK_FLOAT>(12.0) * (cut_ljsq_kk-rsq) *
      (rsq-cut_lj_innersq_kk) * denom_lj_inv_kk;
    const KK_FLOAT philj = lj1 * static_cast<KK_FLOAT>(4.0) * eps *
      (static_cast<KK_FLOAT>(1.0)/(denlj*denlj) - static_cast<KK_FLOAT>(1.0)/denlj);
    forcelj = forcelj*switch1 + philj*switch2;
  }

  return forcelj;
}

/* ----------------------------------------------------------------------
   compute coulomb pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmCoulLongSoftKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& itype, const int& jtype,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  if (rsq >= (STACKPARAMS?m_params[itype][jtype].cut_coulsq:params(itype,jtype).cut_coulsq))
    return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT lj1 = (STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1);
  const KK_FLOAT lj4 = (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4);

  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT grij = g_ewald_kk * r;
  const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
  const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
    t*(static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+t*static_cast<KK_FLOAT>(A5))))) * expm2;

  const KK_FLOAT denc = Kokkos::sqrt(lj4 + rsq);
  const KK_FLOAT prefactor = qqrd2e * lj1 * qtmp * q(j) / (denc*denc*denc);

  KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2);
  if (factor_coul < static_cast<KK_FLOAT>(1.0))
    forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;

  return forcecoul;
}

/* ----------------------------------------------------------------------
   compute LJ 12-6 pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmCoulLongSoftKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  if (rsq >= (STACKPARAMS?m_params[itype][jtype].cut_ljsq:params(itype,jtype).cut_ljsq))
    return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT lj1 = (STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1);
  const KK_FLOAT lj2 = (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2);
  const KK_FLOAT lj3 = (STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3);
  const KK_FLOAT eps = (STACKPARAMS?m_params[itype][jtype].epsilon:params(itype,jtype).epsilon);

  const KK_FLOAT r4sig6 = rsq*rsq / lj2;
  const KK_FLOAT denlj = lj3 + rsq*r4sig6;

  KK_FLOAT evdwl = lj1 * static_cast<KK_FLOAT>(4.0) * eps *
    (static_cast<KK_FLOAT>(1.0)/(denlj*denlj) - static_cast<KK_FLOAT>(1.0)/denlj);

  if (rsq > cut_lj_innersq_kk) {
    const KK_FLOAT switch1 = (cut_ljsq_kk-rsq) * (cut_ljsq_kk-rsq) *
      (cut_ljsq_kk + static_cast<KK_FLOAT>(2.0)*rsq -
       static_cast<KK_FLOAT>(3.0)*cut_lj_innersq_kk) * denom_lj_inv_kk;
    evdwl *= switch1;
  }

  return evdwl;
}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmCoulLongSoftKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& itype, const int& jtype,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  if (rsq >= (STACKPARAMS?m_params[itype][jtype].cut_coulsq:params(itype,jtype).cut_coulsq))
    return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT lj1 = (STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1);
  const KK_FLOAT lj4 = (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4);

  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT grij = g_ewald_kk * r;
  const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
  const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
    t*(static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+t*static_cast<KK_FLOAT>(A5))))) * expm2;

  const KK_FLOAT denc = Kokkos::sqrt(lj4 + rsq);
  const KK_FLOAT prefactor = qqrd2e * lj1 * qtmp * q(j) / denc;

  KK_FLOAT ecoul = prefactor*erfc;
  if (factor_coul < static_cast<KK_FLOAT>(1.0))
    ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;

  return ecoul;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmCoulLongSoftKokkos<DeviceType>::allocate()
{
  PairLJCharmmCoulLongSoft::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  d_cut_ljsq = typename AT::t_kkfloat_2d("pair:cut_ljsq",n+1,n+1);
  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq",n+1,n+1);
  k_params = Kokkos::DualView<params_lj_coul_soft**,Kokkos::LayoutRight,DeviceType>("PairLJCharmmCoulLongSoft::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmCoulLongSoftKokkos<DeviceType>::init_style()
{
  PairLJCharmmCoulLongSoft::init_style();

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

  // the LJ and Coulomb cutoffs are single global values, not per type pair

  Kokkos::deep_copy(d_cut_ljsq,static_cast<KK_FLOAT>(cut_ljsq));
  Kokkos::deep_copy(d_cut_coulsq,static_cast<KK_FLOAT>(cut_coulsq));
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLJCharmmCoulLongSoftKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJCharmmCoulLongSoft::init_one(i,j);
  double cut_ljsqm = cut_ljsq;
  double cut_coulsqm = cut_coulsq;

  k_params.view_host()(i,j).lj1 = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i,j).lj2 = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3 = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4 = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).epsilon = static_cast<KK_FLOAT>(epsilon[i][j]);
  k_params.view_host()(i,j).cut_ljsq = static_cast<KK_FLOAT>(cut_ljsqm);
  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsqm);

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cut_ljsqm);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(cut_coulsqm);
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}



namespace LAMMPS_NS {
template class PairLJCharmmCoulLongSoftKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCharmmCoulLongSoftKokkos<LMPHostType>;
#endif
}

