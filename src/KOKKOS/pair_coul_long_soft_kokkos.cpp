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

#include "pair_coul_long_soft_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCoulLongSoftKokkos<DeviceType>::PairCoulLongSoftKokkos(LAMMPS *lmp) : PairCoulLongSoft(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairCoulLongSoftKokkos<DeviceType>::~PairCoulLongSoftKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq, cutsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongSoftKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_cut_ljsq.template sync<DeviceType>();
  k_cut_coulsq.template sync<DeviceType>();
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
  newton_pair = force->newton_pair;
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

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<PairCoulLongSoftKokkos<DeviceType>,void >
    (this,(NeighListKokkos<DeviceType>*)list);

  if (eflag) eng_coul += static_cast<double>(ev.ecoul);
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

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulLongSoftKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j, const int& itype,
              const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT lam1 = (STACKPARAMS?m_params[itype][jtype].lam1:params(itype,jtype).lam1);
  const KK_FLOAT lam2 = (STACKPARAMS?m_params[itype][jtype].lam2:params(itype,jtype).lam2);

  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT grij = g_ewald_kk * r;
  const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
  const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
    t*(static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+t*static_cast<KK_FLOAT>(A5))))) * expm2;

  const KK_FLOAT denc = Kokkos::sqrt(lam2 + rsq);
  const KK_FLOAT prefactor = qqrd2e * lam1 * qtmp * q(j) / (denc*denc*denc);

  // already a force/r quantity, as in the base style

  KK_FLOAT forcecoul = prefactor * (erfc + static_cast<KK_FLOAT>(EWALD_F)*grij*expm2);
  if (factor_coul < static_cast<KK_FLOAT>(1.0))
    forcecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;

  return forcecoul;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairCoulLongSoftKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int&j, const int& itype,
              const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
  const KK_FLOAT lam1 = (STACKPARAMS?m_params[itype][jtype].lam1:params(itype,jtype).lam1);
  const KK_FLOAT lam2 = (STACKPARAMS?m_params[itype][jtype].lam2:params(itype,jtype).lam2);

  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT grij = g_ewald_kk * r;
  const KK_FLOAT expm2 = Kokkos::exp(-grij*grij);
  const KK_FLOAT t = static_cast<KK_FLOAT>(1.0) /
    (static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(EWALD_P)*grij);
  const KK_FLOAT erfc = t * (static_cast<KK_FLOAT>(A1)+t*(static_cast<KK_FLOAT>(A2)+
    t*(static_cast<KK_FLOAT>(A3)+t*(static_cast<KK_FLOAT>(A4)+t*static_cast<KK_FLOAT>(A5))))) * expm2;

  const KK_FLOAT denc = Kokkos::sqrt(lam2 + rsq);
  const KK_FLOAT prefactor = qqrd2e * lam1 * qtmp * q(j) / denc;

  KK_FLOAT ecoul = prefactor*erfc;
  if (factor_coul < static_cast<KK_FLOAT>(1.0))
    ecoul -= (static_cast<KK_FLOAT>(1.0)-factor_coul)*prefactor;

  return ecoul;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongSoftKokkos<DeviceType>::allocate()
{
  PairCoulLongSoft::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  k_cut_ljsq = DAT::tdual_kkfloat_2d("pair:cut_ljsq",n+1,n+1);
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();
  k_cut_coulsq = DAT::tdual_kkfloat_2d("pair:cut_coulsq",n+1,n+1);
  d_cut_coulsq = k_cut_coulsq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType>("PairCoulLongSoft::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairCoulLongSoftKokkos<DeviceType>::init_style()
{
  PairCoulLongSoft::init_style();

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
double PairCoulLongSoftKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairCoulLongSoft::init_one(i,j);

  k_params.view_host()(i,j).lam1 = static_cast<KK_FLOAT>(lam1[i][j]);
  k_params.view_host()(i,j).lam2 = static_cast<KK_FLOAT>(lam2[i][j]);
  // the Coulomb cutoff is a single global value, not one per type pair

  k_params.view_host()(i,j).cutsq = static_cast<KK_FLOAT>(cut_coulsq);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);

  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cut_coulsq);
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = static_cast<KK_FLOAT>(cut_coulsq);
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = static_cast<KK_FLOAT>(cut_coulsq);
  }
  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cut_coulsq;
  k_cutsq.modify_host();
  k_cut_ljsq.view_host()(i,j) = k_cut_ljsq.view_host()(j,i) = static_cast<KK_FLOAT>(cut_coulsq);
  k_cut_ljsq.modify_host();
  k_cut_coulsq.view_host()(i,j) = k_cut_coulsq.view_host()(j,i) = static_cast<KK_FLOAT>(cut_coulsq);
  k_cut_coulsq.modify_host();
  k_params.modify_host();

  return cutone;
}



namespace LAMMPS_NS {
template class PairCoulLongSoftKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairCoulLongSoftKokkos<LMPHostType>;
#endif
}

