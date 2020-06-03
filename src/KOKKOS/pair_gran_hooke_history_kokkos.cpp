/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_gran_hooke_history_kokkos.h"
#include "kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "error.h"
#include "modify.h"
#include "fix_neigh_history_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairGranHookeHistoryKokkos<Space>::PairGranHookeHistoryKokkos(LAMMPS *lmp) : PairGranHookeHistory(lmp)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = X_MASK | V_MASK | OMEGA_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | MASK_MASK | ENERGY_MASK | VIRIAL_MASK | RMASS_MASK | RADIUS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairGranHookeHistoryKokkos<Space>::~PairGranHookeHistoryKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = NULL;
    vatom = NULL;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairGranHookeHistoryKokkos<Space>::init_style()
{
  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (history && fix_history == NULL) {
    char dnumstr[16];
    sprintf(dnumstr,"%d",3);
    char **fixarg = new char*[4];
    fixarg[0] = (char *) "NEIGH_HISTORY_HH";
    fixarg[1] = (char *) "all";
    if (execution_space == Device)
      fixarg[2] = (char *) "NEIGH_HISTORY/KK/DEVICE";
    else
      fixarg[2] = (char *) "NEIGH_HISTORY/KK/HOST";
    fixarg[3] = dnumstr;
    modify->replace_fix("NEIGH_HISTORY_HH_DUMMY",4,fixarg,1);
    delete [] fixarg;
    int ifix = modify->find_fix("NEIGH_HISTORY_HH");
    fix_history = (FixNeighHistory *) modify->fix[ifix];
    fix_history->pair = this;
    fix_historyKK = (FixNeighHistoryKokkos<Space> *)fix_history;
  }

  PairGranHookeHistory::init_style();

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = (Space == Host) &&
    !(Space == Device);
  neighbor->requests[irequest]->
    kokkos_device = (Space == Device);

  if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Must use half neighbor list with gran/hooke/history/kk");
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairGranHookeHistoryKokkos<Space>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  int shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = DualViewHelper<Space>::view(k_eatom);
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = DualViewHelper<Space>::view(k_vatom);
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  v = DualViewHelper<Space>::view(atomKK->k_v);
  omega = DualViewHelper<Space>::view(atomKK->k_omega);
  c_x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  torque = DualViewHelper<Space>::view(atomKK->k_torque);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  mask = DualViewHelper<Space>::view(atomKK->k_mask);
  tag = DualViewHelper<Space>::view(atomKK->k_tag);
  rmass = DualViewHelper<Space>::view(atomKK->k_rmass);
  radius = DualViewHelper<Space>::view(atomKK->k_radius);
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  int inum = list->inum;
  NeighListKokkos<Space>* k_list = static_cast<NeighListKokkos<Space>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  if (d_numneigh.extent(0) != d_numneigh_touch.extent(0))
    d_numneigh_touch = typename AT::t_int_1d("pair:numneigh_touch",d_numneigh.extent(0));
  if (d_neighbors.extent(0) != d_neighbors_touch.extent(0) ||
      d_neighbors.extent(1) != d_neighbors_touch.extent(1))
    d_neighbors_touch = typename AT::t_neighbors_2d("pair:neighbors_touch",d_neighbors.extent(0),d_neighbors.extent(1));

  d_firsttouch = fix_historyKK->d_firstflag;
  d_firstshear = fix_historyKK->d_firstvalue;

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryReduce>(0,inum),*this);

  EV_FLOAT ev;

  if (lmp->kokkos->neighflag == HALF) {
    if (force->newton_pair) {
      if (vflag_atom) {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,2,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,2,0>>(0,inum),*this);
        }
      } else if (vflag_global) {
        if (shearupdate) {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,1,1>>(0,inum),*this, ev);
        } else {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,1,0>>(0,inum),*this, ev);
        }
      } else {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,0,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,1,0,0>>(0,inum),*this);
        }
      }
    } else {
      if (vflag_atom) {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,2,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,2,0>>(0,inum),*this);
        }
      } else if (vflag_global) {
        if (shearupdate) {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,1,1>>(0,inum),*this, ev);
        } else {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,1,0>>(0,inum),*this, ev);
        }
      } else {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,0,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALF,0,0,0>>(0,inum),*this);
        }
      }
    }
  } else { // HALFTHREAD
    if (force->newton_pair) {
      if (vflag_atom) {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,2,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,2,0>>(0,inum),*this);
        }
      } else if (vflag_global) {
        if (shearupdate) {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,1,1>>(0,inum),*this, ev);
        } else {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,1,0>>(0,inum),*this, ev);
        }
      } else {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,0,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,1,0,0>>(0,inum),*this);
        }
      }
    } else {
      if (vflag_atom) {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,2,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,2,0>>(0,inum),*this);
        }
      } else if (vflag_global) {
        if (shearupdate) {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,1,1>>(0,inum),*this, ev);
        } else {
          Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,1,0>>(0,inum),*this, ev);
        }
      } else {
        if (shearupdate) {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,0,1>>(0,inum),*this);
        } else {
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairGranHookeHistoryCompute<HALFTHREAD,0,0,0>>(0,inum),*this);
        }
      }
    }
  }

  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_atom) {
    DualViewHelper<Space>::modify(k_vatom);
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute<Space>(this);

  copymode = 0;
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void PairGranHookeHistoryKokkos<Space>::operator()(TagPairGranHookeHistoryReduce, const int ii) const {
  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const KK_FLOAT imass = rmass[i];
  const KK_FLOAT irad = radius[i];
  const int jnum = d_numneigh[i];
  int count = 0;

  for (int jj = 0; jj < jnum; jj++) {
    const int j = d_neighbors(i,jj) & NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const KK_FLOAT jmass = rmass[j];
    const KK_FLOAT jrad = radius[j];
    const KK_FLOAT radsum = irad + jrad;

    // check for touching neighbors

    if (rsq >= radsum * radsum) {
      d_firsttouch(i,jj) = 0;
      d_firstshear(i,3*jj) = 0;
      d_firstshear(i,3*jj+1) = 0;
      d_firstshear(i,3*jj+2) = 0;
    } else {
      d_firsttouch(i,jj) = 1;
      d_neighbors_touch(i,count++) = jj;
    }
  }
  d_numneigh_touch[i] = count;
}

template<ExecutionSpace Space>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int SHEARUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranHookeHistoryKokkos<Space>::operator()(TagPairGranHookeHistoryCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,SHEARUPDATE>, const int ii, EV_FLOAT &ev) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<typename AT::t_float_1d_3::data_type, typename AT::t_float_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<typename AT::t_float_1d_3::data_type, typename AT::t_float_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const KK_FLOAT imass = rmass[i];
  const KK_FLOAT irad = radius[i];
  const int jnum = d_numneigh_touch[i];

  KK_FLOAT fx_i = 0.0;
  KK_FLOAT fy_i = 0.0;
  KK_FLOAT fz_i = 0.0;

  KK_FLOAT torquex_i = 0.0;
  KK_FLOAT torquey_i = 0.0;
  KK_FLOAT torquez_i = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    const int m = d_neighbors_touch(i, jj);
    const int j = d_neighbors(i, m) & NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const KK_FLOAT jmass = rmass[j];
    const KK_FLOAT jrad = radius[j];
    const KK_FLOAT radsum = irad + jrad;

    // check for touching neighbors

    const KK_FLOAT r = sqrt(rsq);
    const KK_FLOAT rinv = 1.0/r;
    const KK_FLOAT rsqinv = 1/rsq;

    // relative translational velocity

    KK_FLOAT vr1 = v(i,0) - v(j,0);
    KK_FLOAT vr2 = v(i,1) - v(j,1);
    KK_FLOAT vr3 = v(i,2) - v(j,2);

    // normal component

    KK_FLOAT vnnr = vr1*delx + vr2*dely + vr3*delz;
    KK_FLOAT vn1 = delx*vnnr * rsqinv;
    KK_FLOAT vn2 = dely*vnnr * rsqinv;
    KK_FLOAT vn3 = delz*vnnr * rsqinv;

    // tangential component

    KK_FLOAT vt1 = vr1 - vn1;
    KK_FLOAT vt2 = vr2 - vn2;
    KK_FLOAT vt3 = vr3 - vn3;

    // relative rotational velocity

    KK_FLOAT wr1 = (irad*omega(i,0) + jrad*omega(j,0)) * rinv;
    KK_FLOAT wr2 = (irad*omega(i,1) + jrad*omega(j,1)) * rinv;
    KK_FLOAT wr3 = (irad*omega(i,2) + jrad*omega(j,2)) * rinv;

    KK_FLOAT meff = imass*jmass / (imass+jmass);
    if (mask[i] & freeze_group_bit) meff = jmass;
    if (mask[j] & freeze_group_bit) meff = imass;

    KK_FLOAT damp = meff*gamman*vnnr*rsqinv;
    KK_FLOAT ccel = kn*(radsum-r)*rinv - damp;

    // relative velocities

    KK_FLOAT vtr1 = vt1 - (delz*wr2-dely*wr3);
    KK_FLOAT vtr2 = vt2 - (delx*wr3-delz*wr1);
    KK_FLOAT vtr3 = vt3 - (dely*wr1-delx*wr2);
    KK_FLOAT vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
    vrel = sqrt(vrel);

    // shear history effects

    KK_FLOAT shear1 = d_firstshear(i,3*m);
    KK_FLOAT shear2 = d_firstshear(i,3*m+1);
    KK_FLOAT shear3 = d_firstshear(i,3*m+2);
    if (SHEARUPDATE) {
      shear1 += vtr1*dt;
      shear2 += vtr2*dt;
      shear3 += vtr3*dt;
    }
    KK_FLOAT shrmag = sqrt(shear1*shear1 + shear2*shear2 +
                          shear3*shear3);

    // rotate shear displacements

    KK_FLOAT rsht = shear1*delx + shear2*dely + shear3*delz;
    rsht *= rsqinv;
    if (SHEARUPDATE) {
      shear1 -= rsht*delx;
      shear2 -= rsht*dely;
      shear3 -= rsht*delz;
    }

    // tangential forces = shear + tangential velocity damping

    KK_FLOAT fs1 = - (kt*shear1 + meff*gammat*vtr1);
    KK_FLOAT fs2 = - (kt*shear2 + meff*gammat*vtr2);
    KK_FLOAT fs3 = - (kt*shear3 + meff*gammat*vtr3);

    // rescale frictional displacements and forces if needed

    KK_FLOAT fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
    KK_FLOAT fn = xmu * fabs(ccel*r);

    if (fs > fn) {
      if (shrmag != 0.0) {
        shear1 = (fn/fs) * (shear1 + meff*gammat*vtr1/kt) -
          meff*gammat*vtr1/kt;
        shear2 = (fn/fs) * (shear2 + meff*gammat*vtr2/kt) -
          meff*gammat*vtr2/kt;
        shear3 = (fn/fs) * (shear3 + meff*gammat*vtr3/kt) -
          meff*gammat*vtr3/kt;
        fs1 *= fn/fs;
        fs2 *= fn/fs;
        fs3 *= fn/fs;
      } else fs1 = fs2 = fs3 = 0.0;
    }

    if (SHEARUPDATE) {
      d_firstshear(i,3*m) = shear1;
      d_firstshear(i,3*m+1) = shear2;
      d_firstshear(i,3*m+2) = shear3;
    }

    // forces & torques

    KK_FLOAT fx = delx*ccel + fs1;
    KK_FLOAT fy = dely*ccel + fs2;
    KK_FLOAT fz = delz*ccel + fs3;
    fx_i += fx;
    fy_i += fy;
    fz_i += fz;

    KK_FLOAT tor1 = rinv * (dely*fs3 - delz*fs2);
    KK_FLOAT tor2 = rinv * (delz*fs1 - delx*fs3);
    KK_FLOAT tor3 = rinv * (delx*fs2 - dely*fs1);
    torquex_i -= irad*tor1;
    torquey_i -= irad*tor2;
    torquez_i -= irad*tor3;

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


template<ExecutionSpace Space>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int SHEARUPDATE>
KOKKOS_INLINE_FUNCTION
void PairGranHookeHistoryKokkos<Space>::operator()(TagPairGranHookeHistoryCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,SHEARUPDATE>, const int ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG,SHEARUPDATE>(TagPairGranHookeHistoryCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,SHEARUPDATE>(), ii, ev);
}

template<ExecutionSpace Space>
template<int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairGranHookeHistoryKokkos<Space>::ev_tally_xyz(EV_FLOAT &ev, int i, int j,
                                                          KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                                                          KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const
{
  KK_FLOAT v[6];

  v[0] = delx*fx;
  v[1] = dely*fy;
  v[2] = delz*fz;
  v[3] = delx*fy;
  v[4] = delx*fz;
  v[5] = dely*fz;

  if (NEWTON_PAIR) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  } else {
    if (i < nlocal) {
      ev.v[0] += 0.5*v[0];
      ev.v[1] += 0.5*v[1];
      ev.v[2] += 0.5*v[2];
      ev.v[3] += 0.5*v[3];
      ev.v[4] += 0.5*v[4];
      ev.v[5] += 0.5*v[5];
    }
    if (j < nlocal) {
      ev.v[0] += 0.5*v[0];
      ev.v[1] += 0.5*v[1];
      ev.v[2] += 0.5*v[2];
      ev.v[3] += 0.5*v[3];
      ev.v[4] += 0.5*v[4];
      ev.v[5] += 0.5*v[5];
    }
  }
}

template<ExecutionSpace Space>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairGranHookeHistoryKokkos<Space>::ev_tally_xyz_atom(EV_FLOAT &ev, int i, int j,
                                                               KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                                                               KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const
{
  Kokkos::View<typename AT::t_float_1d_6::data_type, typename AT::t_float_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = DualViewHelper<Space>::view(k_vatom);

  KK_FLOAT v[6];

  v[0] = delx*fx;
  v[1] = dely*fy;
  v[2] = delz*fz;
  v[3] = delx*fy;
  v[4] = delx*fz;
  v[5] = dely*fz;

  if (NEWTON_PAIR || i < nlocal) {
    v_vatom(i,0) += 0.5*v[0];
    v_vatom(i,1) += 0.5*v[1];
    v_vatom(i,2) += 0.5*v[2];
    v_vatom(i,3) += 0.5*v[3];
    v_vatom(i,4) += 0.5*v[4];
    v_vatom(i,5) += 0.5*v[5];
  }
  if (NEWTON_PAIR || j < nlocal) {
    v_vatom(j,0) += 0.5*v[0];
    v_vatom(j,1) += 0.5*v[1];
    v_vatom(j,2) += 0.5*v[2];
    v_vatom(j,3) += 0.5*v[3];
    v_vatom(j,4) += 0.5*v[4];
    v_vatom(j,5) += 0.5*v[5];
  }
}

namespace LAMMPS_NS {
template class PairGranHookeHistoryKokkos<Device>;
template class PairGranHookeHistoryKokkos<Host>;
}
