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

#include "pair_lj_expand_sphere_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJExpandSphereKokkos<DeviceType>::PairLJExpandSphereKokkos(LAMMPS *lmp) : PairLJExpandSphere(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | RADIUS_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJExpandSphereKokkos<DeviceType>::~PairLJExpandSphereKokkos()
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
void PairLJExpandSphereKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  type = atomKK->k_type.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = static_cast<KK_FLOAT>(force->special_lj[0]);
  special_lj[1] = static_cast<KK_FLOAT>(force->special_lj[1]);
  special_lj[2] = static_cast<KK_FLOAT>(force->special_lj[2]);
  special_lj[3] = static_cast<KK_FLOAT>(force->special_lj[3]);

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev = pair_compute<PairLJExpandSphereKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

  if (eflag_global) eng_vdwl += static_cast<double>(ev.evdwl);
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

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJExpandSphereKokkos<DeviceType>::
compute_fpair(const KK_FLOAT &rsq, const int &i, const int &j, const int &itype, const int &jtype) const {
  // cutsq is the maximum possible cutoff for the type pair; the real
  // cutoff depends on the two particle radii and is applied here

  const KK_FLOAT cut = (STACKPARAMS?m_params[itype][jtype].cut:params(itype,jtype).cut);

  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT rshift = r - radius(i) - radius(j);
  if (rshift >= cut) return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT lj1 = (STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1);
  const KK_FLOAT lj2 = (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2);

  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / (rshift*rshift);
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;

  return r6inv*(lj1*r6inv - lj2) * rshift / r;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJExpandSphereKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int &i, const int &j, const int &itype, const int &jtype) const {
  const KK_FLOAT cut = (STACKPARAMS?m_params[itype][jtype].cut:params(itype,jtype).cut);

  const KK_FLOAT r = Kokkos::sqrt(rsq);
  const KK_FLOAT rshift = r - radius(i) - radius(j);
  if (rshift >= cut) return static_cast<KK_FLOAT>(0.0);

  const KK_FLOAT lj3 = (STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3);
  const KK_FLOAT lj4 = (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4);

  const KK_FLOAT rshiftsq = rshift*rshift;
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rshiftsq;
  const KK_FLOAT r6inv = r2inv*r2inv*r2inv;

  KK_FLOAT evdwl = r6inv*(lj3*r6inv - lj4);

  if (offset_flag && (rshiftsq > static_cast<KK_FLOAT>(0.0))) {

    // the shift depends on the two radii, so it cannot be precomputed
    // per type pair the way offset[][] is in a plain lj/cut style

    const KK_FLOAT sigma = (STACKPARAMS?m_params[itype][jtype].sigma:params(itype,jtype).sigma);
    const KK_FLOAT eps = (STACKPARAMS?m_params[itype][jtype].epsilon:params(itype,jtype).epsilon);
    const KK_FLOAT ratio = sigma / (cut + radius(i) + radius(j));
    const KK_FLOAT ratio3 = ratio*ratio*ratio;
    const KK_FLOAT ratio6 = ratio3*ratio3;
    evdwl -= static_cast<KK_FLOAT>(4.0)*eps*(ratio6*ratio6 - ratio6);
  }

  return evdwl;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJExpandSphereKokkos<DeviceType>::allocate()
{
  PairLJExpandSphere::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType>("PairLJExpandSphere::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJExpandSphereKokkos<DeviceType>::init_style()
{
  PairLJExpandSphere::init_style();

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
double PairLJExpandSphereKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJExpandSphere::init_one(i,j);

  k_params.view_host()(i,j).lj1 = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i,j).lj2 = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3 = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4 = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).cut = static_cast<KK_FLOAT>(cut[i][j]);
  k_params.view_host()(i,j).epsilon = static_cast<KK_FLOAT>(epsilon[i][j]);
  k_params.view_host()(i,j).sigma = static_cast<KK_FLOAT>(sigma[i][j]);
  k_params.view_host()(i,j).cutsq = static_cast<KK_FLOAT>(cutone*cutone);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJExpandSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJExpandSphereKokkos<LMPHostType>;
#endif
}

