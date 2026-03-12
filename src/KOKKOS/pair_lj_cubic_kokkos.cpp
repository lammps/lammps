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
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include "pair_lj_cubic_kokkos.h"

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

#include "pair_lj_cubic_const.h"

using namespace LAMMPS_NS;
using namespace PairLJCubicConstants;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCubicKokkos<DeviceType>::PairLJCubicKokkos(LAMMPS *lmp):PairLJCubic(lmp)
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
PairLJCubicKokkos<DeviceType>::~PairLJCubicKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
    memoryKK->destroy_kokkos(k_cut_inner,cut_inner);
    memoryKK->destroy_kokkos(k_cut_inner_sq,cut_inner_sq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCubicKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_cut_inner.template sync<DeviceType>();
  k_cut_inner_sq.template sync<DeviceType>();
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
  newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev = pair_compute<PairLJCubicKokkos<DeviceType>,void >(this,(NeighListKokkos<DeviceType>*)list);

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

/* ----------------------------------------------------------------------
   compute LJ cubic pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCubicKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  KK_FLOAT forcelj;
  const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0) / rsq;

  if (rsq <= (STACKPARAMS?m_params[itype][jtype].cut_inner_sq:params(itype,jtype).cut_inner_sq)) {
    const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv *
      ((STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1)*r6inv -
       (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2));
  } else {
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT rmin = (STACKPARAMS?m_params[itype][jtype].sigma:params(itype,jtype).sigma) * RT6TWO;
    const KK_FLOAT t = (r - (STACKPARAMS?m_params[itype][jtype].cut_inner:params(itype,jtype).cut_inner)) / rmin;
    forcelj = (STACKPARAMS?m_params[itype][jtype].epsilon:params(itype,jtype).epsilon) *
      (static_cast<KK_FLOAT>(-DPHIDS) + static_cast<KK_FLOAT>(A3) * t * t / static_cast<KK_FLOAT>(2.0)) * r / rmin;
  }
  return forcelj*r2inv;
}

/* ----------------------------------------------------------------------
   compute LJ cubic pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCubicKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  KK_FLOAT englj;

  const KK_FLOAT r2inv = 1.0/rsq;

  if (rsq <= (STACKPARAMS?m_params[itype][jtype].cut_inner_sq:params(itype,jtype).cut_inner_sq)) {
    const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
    englj = r6inv *
      ((STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*r6inv -
       (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4));
  } else {
    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT rmin = (STACKPARAMS?m_params[itype][jtype].sigma:params(itype,jtype).sigma) * RT6TWO;
    const KK_FLOAT t = (r - (STACKPARAMS?m_params[itype][jtype].cut_inner:params(itype,jtype).cut_inner)) / rmin;
    englj = (STACKPARAMS?m_params[itype][jtype].epsilon:params(itype,jtype).epsilon) *
      (static_cast<KK_FLOAT>(PHIS) + static_cast<KK_FLOAT>(DPHIDS) * t - static_cast<KK_FLOAT>(A3) * t * t * t/ static_cast<KK_FLOAT>(6.0));
  }
  return englj;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCubicKokkos<DeviceType>::allocate()
{
  PairLJCubic::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  memory->destroy(cut_inner);
  memoryKK->create_kokkos(k_cut_inner,cut_inner,n+1,n+1,"pair:cut_inner");
  d_cut_inner = k_cut_inner.template view<DeviceType>();

  memory->destroy(cut_inner_sq);
  memoryKK->create_kokkos(k_cut_inner_sq,cut_inner_sq,n+1,n+1,"pair:cut_inner_sq");
  d_cut_inner_sq = k_cut_inner_sq.template view<DeviceType>();

  k_params = Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType>("PairLJCubic::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCubicKokkos<DeviceType>::init_style()
{
  PairLJCubic::init_style();

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
double PairLJCubicKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJCubic::init_one(i,j);
  double cut_inner_sqm = cut_inner_sq[i][j];

  k_params.view_host()(i,j).lj1 = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i,j).lj2 = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3 = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4 = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).cut_inner_sq = cut_inner_sqm;
  k_params.view_host()(i,j).cut_inner = static_cast<KK_FLOAT>(cut_inner[i][j]);
  k_params.view_host()(i,j).epsilon = static_cast<KK_FLOAT>(epsilon[i][j]);
  k_params.view_host()(i,j).sigma = static_cast<KK_FLOAT>(sigma[i][j]);

  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
    m_cut_inner_sq[j][i] = m_cut_inner_sq[i][j] = cut_inner_sqm;
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = static_cast<KK_FLOAT>(cutone*cutone);
  k_cutsq.modify_host();
  k_cut_inner_sq.view_host()(i,j) = k_cut_inner_sq.view_host()(j,i) = cut_inner_sqm;
  k_cut_inner_sq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJCubicKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCubicKokkos<LMPHostType>;
#endif
}

