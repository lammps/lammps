// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_lj_cut_dipole_cut_kokkos.h"

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

#define KOKKOS_CUDA_MAX_THREADS 256
#define KOKKOS_CUDA_MIN_BLOCKS 8

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutDipoleCutKokkos<DeviceType>::PairLJCutDipoleCutKokkos(LAMMPS *lmp):PairLJCutDipoleCut(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutDipoleCutKokkos<DeviceType>::~PairLJCutDipoleCutKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
    memoryKK->destroy_kokkos(k_cut_ljsq,cut_ljsq);
    memoryKK->destroy_kokkos(k_cut_coulsq,cut_coulsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutDipoleCutKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  special_coul[0] = force->special_coul[0];
  special_coul[1] = force->special_coul[1];
  special_coul[2] = force->special_coul[2];
  special_coul[3] = force->special_coul[3];
  qqrd2e = force->qqrd2e;
  newton_pair = force->newton_pair;

  // loop over neighbors of my atoms
/*
  EV_FLOAT ev = pair_compute<PairLJCutCoulCutKokkos<DeviceType>,void >
    (this,(NeighListKokkos<DeviceType>*)list);

*/
  EV_FLOAT ev;

  if (neighflag == HALF) {
    if (force->newton_pair) {
      if (vflag_atom) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,inum),*this);
      } else if (vflag_global) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,inum),*this, ev);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,inum),*this);
      }
    } else {
      if (vflag_atom) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,inum),*this);
      } else if (vflag_global) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,inum),*this, ev);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,inum),*this);
      }
    }
  } else { // HALFTHREAD
    if (force->newton_pair) {
      if (vflag_atom) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,inum),*this);
      } else if (vflag_global) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,inum),*this, ev);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,inum),*this);
      }
    } else {
      if (vflag_atom) {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,inum),*this);
      } else if (vflag_global) {
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,inum),*this, ev);
      } else {
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType>(0,inum),*this);        
      }
    }
  }

  if (eflag) {
    eng_vdwl += ev.evdwl;
    eng_coul += ev.ecoul;
  }

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
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairLJCutDipoleCutKokkos<DeviceType>::operator()(const int ii, EV_FLOAT &ev) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  const int i = list.d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  
  const int itype = c.type(i);

  const AtomNeighborsConst neighbors_i = list.get_neighbors_const(i);
  const int jnum = list.d_numneigh[i];

  F_FLOAT fx_i = 0.0;
  F_FLOAT fy_i = 0.0;
  F_FLOAT fz_i = 0.0;

  F_FLOAT torquex_i = 0.0;
  F_FLOAT torquey_i = 0.0;
  F_FLOAT torquez_i = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = neighbors_i(jj);
    const F_FLOAT factor_lj = c.special_lj[sbmask(j)];
    j &= NEIGHMASK;

    const X_FLOAT delx = xtmp - c.x(j,0);
    const X_FLOAT dely = ytmp - c.x(j,1);
    const X_FLOAT delz = ztmp - c.x(j,2);
    const int jtype = c.type(j);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    // forces & torques
/*
    F_FLOAT fx = delx*ccel + fs1;
    F_FLOAT fy = dely*ccel + fs2;
    F_FLOAT fz = delz*ccel + fs3;
    fx_i += fx;
    fy_i += fy;
    fz_i += fz;

    F_FLOAT tor1 = rinv * (dely*fs3 - delz*fs2);
    F_FLOAT tor2 = rinv * (delz*fs1 - delx*fs3);
    F_FLOAT tor3 = rinv * (delx*fs2 - dely*fs1);
    torquex_i -= irad*tor1;
    torquey_i -= irad*tor2;
    torquez_i -= irad*tor3;
*/
    if (NEWTON_PAIR || j < nlocal) {
      a_f(j,0) -= fx;
      a_f(j,1) -= fy;
      a_f(j,2) -= fz;
      a_torque(j,0) -= jrad*tor1;
      a_torque(j,1) -= jrad*tor2;
      a_torque(j,2) -= jrad*tor3;
    }

    if (EVFLAG == 2)
      ev_tally_xyz_atom<NEIGHFLAG, NEWTON_PAIR>(ev, i, j, fx_i, fy_i, fz_i, delx, dely, delz);
    if (EVFLAG == 1)
      ev_tally_xyz<NEWTON_PAIR>(ev, i, j, fx_i, fy_i, fz_i, delx, dely, delz);
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
  a_torque(i,0) += torquex_i;
  a_torque(i,1) += torquey_i;
  a_torque(i,2) += torquez_i;
}

/* ----------------------------------------------------------------------
   compute LJ 12-6 pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCutDipoleCutKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;
  F_FLOAT forcelj;

  forcelj = r6inv *
    ((STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1)*r6inv -
     (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2));

  return forcelj*r2inv;
}

/* ----------------------------------------------------------------------
   compute coulomb pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCutDipoleCutKokkos<DeviceType>::
compute_fcoul(const F_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/,
              const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT rinv = sqrt(r2inv);
  F_FLOAT forcecoul;

  forcecoul = qqrd2e*qtmp*q(j) *rinv;

  return factor_coul*forcecoul*r2inv;
}

/* ----------------------------------------------------------------------
   compute LJ 12-6 pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCutDipoleCutKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;

  return r6inv*
    ((STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*r6inv
     - (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4))
    -  (STACKPARAMS?m_params[itype][jtype].offset:params(itype,jtype).offset);

}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCutDipoleCutKokkos<DeviceType>::
compute_ecoul(const F_FLOAT& rsq, const int& /*i*/, const int&j,
              const int& /*itype*/, const int& /*jtype*/,
              const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT rinv = sqrt(r2inv);

  return factor_coul*qqrd2e*qtmp*q(j)*rinv;

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutDipoleCutKokkos<DeviceType>::allocate()
{
  PairLJCutDipoleCut::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  memory->destroy(cut_ljsq);
  memoryKK->create_kokkos(k_cut_ljsq,cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();
  memory->destroy(cut_coulsq);
  memoryKK->create_kokkos(k_cut_coulsq,cut_coulsq,n+1,n+1,"pair:cut_coulsq");
  d_cut_coulsq = k_cut_coulsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_lj_coul**,Kokkos::LayoutRight,DeviceType>("PairLJCutDipoleCut::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutDipoleCutKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Illegal pair_style command");

  PairLJCutDipoleCut::settings(1,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutDipoleCutKokkos<DeviceType>::init_style()
{
  PairLJCutDipoleCut::init_style();

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
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLJCutDipoleCutKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJCutDipoleCut::init_one(i,j);
  double cut_ljsqm = cut_ljsq[i][j];
  double cut_coulsqm = cut_coulsq[i][j];

  k_params.h_view(i,j).lj1 = lj1[i][j];
  k_params.h_view(i,j).lj2 = lj2[i][j];
  k_params.h_view(i,j).lj3 = lj3[i][j];
  k_params.h_view(i,j).lj4 = lj4[i][j];
  k_params.h_view(i,j).offset = offset[i][j];
  k_params.h_view(i,j).cut_ljsq = cut_ljsqm;
  k_params.h_view(i,j).cut_coulsq = cut_coulsqm;

  k_params.h_view(j,i) = k_params.h_view(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = cut_ljsqm;
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = cut_coulsqm;
  }

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  k_cut_ljsq.h_view(i,j) = k_cut_ljsq.h_view(j,i) = cut_ljsqm;
  k_cut_ljsq.template modify<LMPHostType>();
  k_cut_coulsq.h_view(i,j) = k_cut_coulsq.h_view(j,i) = cut_coulsqm;
  k_cut_coulsq.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();

  return cutone;
}



namespace LAMMPS_NS {
template class PairLJCutDipoleCutKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutDipoleCutKokkos<LMPHostType>;
#endif
}

