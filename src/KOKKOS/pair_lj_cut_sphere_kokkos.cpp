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
   (Kokkos version) LAMMPS development team
------------------------------------------------------------------------- */

#include "pair_lj_cut_sphere_kokkos.h"

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

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutSphereKokkos<DeviceType>::PairLJCutSphereKokkos(LAMMPS *lmp)
  : PairLJCutSphere(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK | RADIUS_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutSphereKokkos<DeviceType>::~PairLJCutSphereKokkos()
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
void PairLJCutSphereKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  m_mix_flag   = mix_flag;
  m_offset_flag = offset_flag;

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  newton_pair = force->newton_pair;

  EV_FLOAT ev;

  copymode = 1;

  ev = pair_compute<PairLJCutSphereKokkos<DeviceType>,void>
    (this,(NeighListKokkos<DeviceType>*)list);

  if (eflag) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
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
   LJ force using per-atom radius; mixing rule applied on-the-fly
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCutSphereKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& i, const int& j,
              const int& itype, const int& jtype) const
{
  const KK_FLOAT ri = radius(i);
  const KK_FLOAT rj = radius(j);

  // compute mixed diameter sigma from per-atom radii
  KK_FLOAT sigma;
  if (m_mix_flag == GEOMETRIC)
    sigma = static_cast<KK_FLOAT>(2.0) * Kokkos::sqrt(ri * rj);
  else   // ARITHMETIC (SIXTHPOWER is forbidden for this pair style)
    sigma = ri + rj;

  // actual per-atom cutoff (cutsq in params is the max-possible cutoff)
  const KK_FLOAT cut_ij = STACKPARAMS ? m_params[itype][jtype].cut : params(itype,jtype).cut;
  const KK_FLOAT rcutsq = cut_ij * cut_ij * sigma * sigma;
  if (rsq >= rcutsq) return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT r2inv  = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv  = r2inv * r2inv * r2inv;
  const KK_FLOAT sigma2 = sigma * sigma;
  const KK_FLOAT sigma6 = sigma2 * sigma2 * sigma2;
  const KK_FLOAT epsilon_val = STACKPARAMS ? m_params[itype][jtype].epsilon : params(itype,jtype).epsilon;
  const KK_FLOAT forcelj = r6inv * static_cast<KK_FLOAT>(24.0) * epsilon_val *
    (static_cast<KK_FLOAT>(2.0) * sigma6 * sigma6 * r6inv - sigma6);
  return forcelj * r2inv;
}

/* ----------------------------------------------------------------------
   LJ energy using per-atom radius
------------------------------------------------------------------------- */

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCutSphereKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& i, const int& j,
               const int& itype, const int& jtype) const
{
  const KK_FLOAT ri = radius(i);
  const KK_FLOAT rj = radius(j);

  KK_FLOAT sigma;
  if (m_mix_flag == GEOMETRIC)
    sigma = static_cast<KK_FLOAT>(2.0) * Kokkos::sqrt(ri * rj);
  else
    sigma = ri + rj;

  const KK_FLOAT cut_ij = STACKPARAMS ? m_params[itype][jtype].cut : params(itype,jtype).cut;
  const KK_FLOAT rcutsq = cut_ij * cut_ij * sigma * sigma;
  if (rsq >= rcutsq) return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT r2inv  = static_cast<KK_FLOAT>(1.0) / rsq;
  const KK_FLOAT r6inv  = r2inv * r2inv * r2inv;
  const KK_FLOAT sigma2 = sigma * sigma;
  const KK_FLOAT sigma6 = sigma2 * sigma2 * sigma2;
  const KK_FLOAT epsilon_val = STACKPARAMS ? m_params[itype][jtype].epsilon : params(itype,jtype).epsilon;
  KK_FLOAT evdwl = r6inv * static_cast<KK_FLOAT>(4.0) * epsilon_val *
    (sigma6 * sigma6 * r6inv - sigma6);

  if (m_offset_flag && (rcutsq > static_cast<KK_FLOAT>(0.0))) {
    const KK_FLOAT rcutsq3 = rcutsq * rcutsq * rcutsq;
    const KK_FLOAT ratio6  = sigma6 / rcutsq3;
    evdwl -= static_cast<KK_FLOAT>(4.0) * epsilon_val * (ratio6 * ratio6 - ratio6);
  }
  return evdwl;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutSphereKokkos<DeviceType>::allocate()
{
  PairLJCutSphere::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  k_cut_ljsq = DAT::tdual_kkfloat_2d("pair:cut_ljsq",n+1,n+1);
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_lj_cut_sphere**,Kokkos::LayoutRight,DeviceType>(
    "PairLJCutSphere::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutSphereKokkos<DeviceType>::init_style()
{
  PairLJCutSphere::init_style();

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
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
   NOTE: init_one returns the MAX possible cutoff (cut[i][j]*2*mix(rmax_i,rmax_j))
         which is stored in cutsq.  The actual per-atom cutoff check happens
         inside compute_fpair using params.cut (the plain multiplier).
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLJCutSphereKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJCutSphere::init_one(i,j);
  double cutsqm = cutone*cutone;

  k_params.view_host()(i,j).epsilon = static_cast<KK_FLOAT>(epsilon[i][j]);
  k_params.view_host()(i,j).cut     = static_cast<KK_FLOAT>(cut[i][j]);
  k_params.view_host()(i,j).cutsq   = static_cast<KK_FLOAT>(cutsqm);

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i]    = m_cutsq[i][j]    = cutsqm;
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = cutsqm;
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutsqm;
  k_cutsq.modify_host();
  k_cut_ljsq.view_host()(i,j) = k_cut_ljsq.view_host()(j,i) = cutsqm;
  k_cut_ljsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJCutSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutSphereKokkos<LMPHostType>;
#endif
}
