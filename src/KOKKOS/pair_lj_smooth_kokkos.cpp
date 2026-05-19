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

#include "pair_lj_smooth_kokkos.h"

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

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJSmoothKokkos<DeviceType>::PairLJSmoothKokkos(LAMMPS *lmp) : PairLJSmooth(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJSmoothKokkos<DeviceType>::~PairLJSmoothKokkos()
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
void PairLJSmoothKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  copymode = 1;
  EV_FLOAT ev = pair_compute<PairLJSmoothKokkos<DeviceType>,void>(this,(NeighListKokkos<DeviceType>*)list);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  copymode = 0;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJSmoothKokkos<DeviceType>::
compute_fpair(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT cut_inner_sq = STACKPARAMS ? m_params[itype][jtype].cut_inner_sq : params(itype,jtype).cut_inner_sq;
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
  KK_FLOAT forcelj;
  if (rsq < cut_inner_sq) {
    const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
    const KK_FLOAT lj1 = STACKPARAMS ? m_params[itype][jtype].lj1 : params(itype,jtype).lj1;
    const KK_FLOAT lj2 = STACKPARAMS ? m_params[itype][jtype].lj2 : params(itype,jtype).lj2;
    forcelj = r6inv * (lj1*r6inv - lj2);
  } else {
    const KK_FLOAT cut_inner = STACKPARAMS ? m_params[itype][jtype].cut_inner : params(itype,jtype).cut_inner;
    const KK_FLOAT ljsw1 = STACKPARAMS ? m_params[itype][jtype].ljsw1 : params(itype,jtype).ljsw1;
    const KK_FLOAT ljsw2 = STACKPARAMS ? m_params[itype][jtype].ljsw2 : params(itype,jtype).ljsw2;
    const KK_FLOAT ljsw3 = STACKPARAMS ? m_params[itype][jtype].ljsw3 : params(itype,jtype).ljsw3;
    const KK_FLOAT ljsw4 = STACKPARAMS ? m_params[itype][jtype].ljsw4 : params(itype,jtype).ljsw4;
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT t = r - cut_inner;
    const KK_FLOAT tsq = t*t;
    const KK_FLOAT fskin = ljsw1 + ljsw2*t + ljsw3*tsq + ljsw4*tsq*t;
    forcelj = fskin*r;
  }
  return forcelj*r2inv;
}

template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJSmoothKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT &rsq, const int &, const int &, const int &itype, const int &jtype) const {
  const KK_FLOAT cut_inner_sq = STACKPARAMS ? m_params[itype][jtype].cut_inner_sq : params(itype,jtype).cut_inner_sq;
  const KK_FLOAT offset = STACKPARAMS ? m_params[itype][jtype].offset : params(itype,jtype).offset;
  if (rsq < cut_inner_sq) {
    const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;
    const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
    const KK_FLOAT lj3 = STACKPARAMS ? m_params[itype][jtype].lj3 : params(itype,jtype).lj3;
    const KK_FLOAT lj4 = STACKPARAMS ? m_params[itype][jtype].lj4 : params(itype,jtype).lj4;
    return r6inv * (lj3*r6inv - lj4) - offset;
  } else {
    const KK_FLOAT cut_inner = STACKPARAMS ? m_params[itype][jtype].cut_inner : params(itype,jtype).cut_inner;
    const KK_FLOAT ljsw0 = STACKPARAMS ? m_params[itype][jtype].ljsw0 : params(itype,jtype).ljsw0;
    const KK_FLOAT ljsw1 = STACKPARAMS ? m_params[itype][jtype].ljsw1 : params(itype,jtype).ljsw1;
    const KK_FLOAT ljsw2 = STACKPARAMS ? m_params[itype][jtype].ljsw2 : params(itype,jtype).ljsw2;
    const KK_FLOAT ljsw3 = STACKPARAMS ? m_params[itype][jtype].ljsw3 : params(itype,jtype).ljsw3;
    const KK_FLOAT ljsw4 = STACKPARAMS ? m_params[itype][jtype].ljsw4 : params(itype,jtype).ljsw4;
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT t = r - cut_inner;
    const KK_FLOAT tsq = t*t;
    return ljsw0 - ljsw1*t - ljsw2*tsq/static_cast<KK_FLOAT>(2.0)
           - ljsw3*tsq*t/static_cast<KK_FLOAT>(3.0)
           - ljsw4*tsq*tsq/static_cast<KK_FLOAT>(4.0) - offset;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJSmoothKokkos<DeviceType>::allocate()
{
  PairLJSmooth::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_lj_smooth**,Kokkos::LayoutRight,DeviceType>("PairLJSmooth::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJSmoothKokkos<DeviceType>::init_style()
{
  PairLJSmooth::init_style();

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
double PairLJSmoothKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJSmooth::init_one(i,j);

  k_params.view_host()(i,j).cut_inner_sq = static_cast<KK_FLOAT>(cut_inner_sq[i][j]);
  k_params.view_host()(i,j).lj1          = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i,j).lj2          = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3          = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4          = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).ljsw0        = static_cast<KK_FLOAT>(ljsw0[i][j]);
  k_params.view_host()(i,j).ljsw1        = static_cast<KK_FLOAT>(ljsw1[i][j]);
  k_params.view_host()(i,j).ljsw2        = static_cast<KK_FLOAT>(ljsw2[i][j]);
  k_params.view_host()(i,j).ljsw3        = static_cast<KK_FLOAT>(ljsw3[i][j]);
  k_params.view_host()(i,j).ljsw4        = static_cast<KK_FLOAT>(ljsw4[i][j]);
  k_params.view_host()(i,j).cut_inner    = static_cast<KK_FLOAT>(cut_inner[i][j]);
  k_params.view_host()(i,j).offset       = static_cast<KK_FLOAT>(offset[i][j]);
  k_params.view_host()(i,j).cutsq        = static_cast<KK_FLOAT>(cutone*cutone);
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
template class PairLJSmoothKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJSmoothKokkos<LMPHostType>;
#endif
}
