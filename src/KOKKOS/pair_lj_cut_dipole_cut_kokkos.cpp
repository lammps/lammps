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
   Contributing author: Trung Nguyen (U Chicago)
------------------------------------------------------------------------- */

#include "pair_lj_cut_dipole_cut_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutDipoleCutKokkos<DeviceType>::PairLJCutDipoleCutKokkos(LAMMPS *lmp) : PairLJCutDipoleCut(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TORQUE_MASK | TYPE_MASK | Q_MASK | MU_MASK | ENERGY_MASK | VIRIAL_MASK;
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
  mu = atomKK->k_mu.view<DeviceType>();
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

  // get the neighbor list and neighbors used in operator()

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  int inum = list->inum;

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;

  // compute kernel NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS
  if (atom->ntypes > MAX_TYPES_STACKPARAMS) { // STACKPARAMS==false
    if (evflag) { // EVFLAG==1
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALF,1,1,false> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALF,0,1,false> >(0,inum),*this,ev);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALFTHREAD,1,1,false> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALFTHREAD,0,1,false> >(0,inum),*this,ev);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<FULL,1,1,false> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<FULL,0,1,false> >(0,inum),*this,ev);
      }
    } else {  // EVFLAG==0
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALF,1,0,false> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALF,0,0,false> >(0,inum),*this);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALFTHREAD,1,0,false> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALFTHREAD,0,0,false> >(0,inum),*this);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<FULL,1,0,false> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<FULL,0,0,false> >(0,inum),*this);
      }
    }
  } else { // STACKPARAMS==true
    if (evflag) { // EVFLAG==1
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALF,1,1,true> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALF,0,1,true> >(0,inum),*this,ev);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALFTHREAD,1,1,true> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALFTHREAD,0,1,true> >(0,inum),*this,ev);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<FULL,1,1,true> >(0,inum),*this,ev);
        else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<FULL,0,1,true> >(0,inum),*this,ev);
      }
    } else {
      if (neighflag == HALF) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALF,1,0,true> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALF,0,0,true> >(0,inum),*this);
      } else if (neighflag == HALFTHREAD) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALFTHREAD,1,0,true> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<HALFTHREAD,0,0,true> >(0,inum),*this);
      } else if (neighflag == FULL) {
        if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<FULL,1,0,true> >(0,inum),*this);
        else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJCutDipoleCutKernel<FULL,0,0,true> >(0,inum),*this);
      }
    }
  }

  if (eflag_global) {
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

  copymode = 0;
}

/* ----------------------------------------------------------------------
  needs torque as with AtomVecSphereKokkos
  needs energy calculation as well
  ---------------------------------------------------------------------- */
template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
KOKKOS_INLINE_FUNCTION
void PairLJCutDipoleCutKokkos<DeviceType>::operator()(TagPairLJCutDipoleCutKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int ii, EV_FLOAT &ev) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const X_FLOAT mui = mu(i,3);
  const int itype = type(i);
  const F_FLOAT qtmp = q[i];

  const int jnum = d_numneigh[i];

  F_FLOAT fx_i = 0.0;
  F_FLOAT fy_i = 0.0;
  F_FLOAT fz_i = 0.0;
  F_FLOAT torquex_i = 0.0;
  F_FLOAT torquey_i = 0.0;
  F_FLOAT torquez_i = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    const F_FLOAT factor_lj = special_lj[sbmask(j)];
    const F_FLOAT factor_coul = special_coul[sbmask(j)];
    j &= NEIGHMASK;

    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const X_FLOAT muj = mu(j,3);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    X_FLOAT cutsq_ij = STACKPARAMS?m_cutsq[itype][jtype]:d_cutsq(itype,jtype);

    if (rsq < cutsq_ij) {
      const F_FLOAT r2inv = 1.0/rsq;
      const F_FLOAT r6inv = r2inv*r2inv*r2inv;
      F_FLOAT forcelj = 0;
      F_FLOAT evdwl = 0;
      F_FLOAT ecoul = 0;
      F_FLOAT forcecoulx = 0;
      F_FLOAT forcecouly = 0;
      F_FLOAT forcecoulz = 0;
      F_FLOAT tixcoul = 0;
      F_FLOAT tiycoul = 0;
      F_FLOAT tizcoul = 0;
      F_FLOAT tjxcoul = 0;
      F_FLOAT tjycoul = 0;
      F_FLOAT tjzcoul = 0;
      F_FLOAT fx = 0;
      F_FLOAT fy = 0;
      F_FLOAT fz = 0;

      // lj term

      X_FLOAT cut_ljsq_ij = STACKPARAMS?m_cut_ljsq[itype][jtype]:d_cut_ljsq(itype,jtype);
      if (rsq < cut_ljsq_ij) {
        forcelj = r6inv * ((STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1)*r6inv -
                           (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2));
        forcelj *= factor_lj*r2inv;
        if (eflag_global) {
          evdwl = r6inv * ((STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*r6inv -
                          (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4)) -
                          (STACKPARAMS?m_params[itype][jtype].offset:params(itype,jtype).offset);
          evdwl *= factor_lj;
          ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*evdwl;
        }
      } // cutsq_ljsq_ij

      // coul term

      X_FLOAT cut_coulsq_ij = STACKPARAMS?m_cut_coulsq[itype][jtype]:d_cut_coulsq(itype,jtype);

      if (rsq < cut_coulsq_ij) {

        const F_FLOAT r2inv = 1.0/rsq;
        const F_FLOAT rinv = sqrt(r2inv);
        const F_FLOAT qj = q[j];

        X_FLOAT r3inv = r2inv*rinv;

        // charge-charge
        if (qtmp != 0.0 && qj != 0.0) {
          X_FLOAT pre1 = qtmp*qj*r3inv;
          forcecoulx += pre1*delx;
          forcecouly += pre1*dely;
          forcecoulz += pre1*delz;
        }

        // dipole-dipole

        F_FLOAT pdotp, pidotr, pjdotr;
        F_FLOAT r5inv = r3inv*r2inv;

       if (mui > 0.0 && muj > 0.0) {

          F_FLOAT r7inv = r5inv*r2inv;

          pdotp = mu(i,0)*mu(j,0) + mu(i,1)*mu(j,1) + mu(i,2)*mu(j,2);
          pidotr = mu(i,0)*delx + mu(i,1)*dely + mu(i,2)*delz;
          pjdotr = mu(j,0)*delx + mu(j,1)*dely + mu(j,2)*delz;

          X_FLOAT pre1 = 3.0*r5inv*pdotp - 15.0*r7inv*pidotr*pjdotr;
          F_FLOAT pre2 = 3.0*r5inv*pjdotr;
          F_FLOAT pre3 = 3.0*r5inv*pidotr;
          F_FLOAT pre4 = -1.0*r3inv;

          forcecoulx += pre1*delx + pre2*mu(i,0) + pre3*mu(j,0);
          forcecouly += pre1*dely + pre2*mu(i,1) + pre3*mu(j,1);
          forcecoulz += pre1*delz + pre2*mu(i,2) + pre3*mu(j,2);

          F_FLOAT crossx = pre4 * (mu(i,1)*mu(j,2) - mu(i,2)*mu(j,1));
          F_FLOAT crossy = pre4 * (mu(i,2)*mu(j,0) - mu(i,0)*mu(j,2));
          F_FLOAT crossz = pre4 * (mu(i,0)*mu(j,1) - mu(i,1)*mu(j,0));

          tixcoul += crossx + pre2 * (mu(i,1)*delz - mu(i,2)*dely);
          tiycoul += crossy + pre2 * (mu(i,2)*delx - mu(i,0)*delz);
          tizcoul += crossz + pre2 * (mu(i,0)*dely - mu(i,1)*delx);
          tjxcoul += -crossx + pre3 * (mu(j,1)*delz - mu(j,2)*dely);
          tjycoul += -crossy + pre3 * (mu(j,2)*delx - mu(j,0)*delz);
          tjzcoul += -crossz + pre3 * (mu(j,0)*dely - mu(j,1)*delx);
        }

        // dipole-charge

        if (mui > 0 && qj != 0) {
          pidotr = mu(i,0)*delx + mu(i,1)*dely + mu(i,2)*delz;
          F_FLOAT pre1 = 3.0*qj*r5inv * pidotr;
          F_FLOAT pre2 = qj*r3inv;

          forcecoulx += pre2*mu(i,0) - pre1*delx;
          forcecouly += pre2*mu(i,1) - pre1*dely;
          forcecoulz += pre2*mu(i,2) - pre1*delz;
          tixcoul += pre2 * (mu(i,1)*delz - mu(i,2)*dely);
          tiycoul += pre2 * (mu(i,2)*delx - mu(i,0)*delz);
          tizcoul += pre2 * (mu(i,0)*dely - mu(i,1)*delx);
        }

        // charge-dipole

        if (qtmp != 0 && muj > 0) {
          pjdotr = mu(j,0)*delx + mu(j,1)*dely + mu(j,2)*delz;
          X_FLOAT pre1 = 3.0*qtmp*r5inv * pjdotr;
          X_FLOAT pre2 = qtmp*r3inv;

          forcecoulx += pre1*delx - pre2*mu(j,0);
          forcecouly += pre1*dely - pre2*mu(j,1);
          forcecoulz += pre1*delz - pre2*mu(j,2);
          tjxcoul += -pre2 * (mu(j,1)*delz - mu(j,2)*dely);
          tjycoul += -pre2 * (mu(j,2)*delx - mu(j,0)*delz);
          tjzcoul += -pre2 * (mu(j,0)*dely - mu(j,1)*delx);
        }

        F_FLOAT fq = factor_coul*qqrd2e;
        fx = fq*forcecoulx + delx*forcelj;
        fy = fq*forcecouly + dely*forcelj;
        fz = fq*forcecoulz + delz*forcelj;

        // force & torque accumulation

        fx_i += fx;
        fy_i += fy;
        fz_i += fz;
        torquex_i += fq*tixcoul;
        torquey_i += fq*tiycoul;
        torquez_i += fq*tizcoul;

        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
          a_f(j,0) -= fx;
          a_f(j,1) -= fy;
          a_f(j,2) -= fz;
          a_torque(j,0) += fq*tjxcoul;
          a_torque(j,1) += fq*tjycoul;
          a_torque(j,2) += fq*tjzcoul;
        }

        if (EVFLAG && eflag_global) {
          ecoul = qtmp*qj*rinv;
          if (mu(i,3) > 0.0 && mu(j,3) > 0.0)
            ecoul += r3inv*pdotp - 3.0*r5inv*pidotr*pjdotr;
          if (mu(i,3) > 0.0 && qj != 0.0)
            ecoul += -qj*r3inv*pidotr;
          if (mu(j,3) > 0.0 && qtmp != 0.0)
            ecoul += qtmp*r3inv*pjdotr;
          ecoul *= factor_coul*qqrd2e;
          ev.ecoul += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*ecoul;
        }
      } // cutsq_coulsq_ij

      if (EVFLAG && (eflag_atom || vflag_either))
        ev_tally_xyz<NEIGHFLAG, NEWTON_PAIR>(ev, i, j, ecoul+evdwl, fx, fy, fz, delx, dely, delz);
    } // cutsq_ij
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
  a_torque(i,0) += torquex_i;
  a_torque(i,1) += torquey_i;
  a_torque(i,2) += torquez_i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool STACKPARAMS>
KOKKOS_INLINE_FUNCTION
void PairLJCutDipoleCutKokkos<DeviceType>::operator()(TagPairLJCutDipoleCutKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>, const int ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairLJCutDipoleCutKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairLJCutDipoleCutKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT & ev, int i, int j, const F_FLOAT &epair,
                                                        F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
                                                        X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const
{
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();

  if (eflag_atom) {
    const E_FLOAT epairhalf = 0.5 * epair;
    if (NEIGHFLAG == FULL || newton_pair || i < nlocal)
      v_eatom[i] += epairhalf;
    if (NEIGHFLAG != FULL && (newton_pair || j < nlocal))
      v_eatom[j] += epairhalf;
  }

  if (vflag_either) {
    const F_FLOAT v0 = delx*fx;
    const F_FLOAT v1 = dely*fy;
    const F_FLOAT v2 = delz*fz;
    const F_FLOAT v3 = delx*fy;
    const F_FLOAT v4 = delx*fz;
    const F_FLOAT v5 = dely*fz;

    if (vflag_global) {
      if (NEIGHFLAG != FULL) {
        if (NEWTON_PAIR) { // neigh half, newton on
          ev.v[0] += v0;
          ev.v[1] += v1;
          ev.v[2] += v2;
          ev.v[3] += v3;
          ev.v[4] += v4;
          ev.v[5] += v5;
        } else { // neigh half, newton off
          if (i < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
          if (j < nlocal) {
            ev.v[0] += 0.5*v0;
            ev.v[1] += 0.5*v1;
            ev.v[2] += 0.5*v2;
            ev.v[3] += 0.5*v3;
            ev.v[4] += 0.5*v4;
            ev.v[5] += 0.5*v5;
          }
        }
      } else { //neigh full
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
      }
    }

    if (vflag_atom) {

      if (NEIGHFLAG == FULL || NEWTON_PAIR || i < nlocal) {
        v_vatom(i,0) += 0.5*v0;
        v_vatom(i,1) += 0.5*v1;
        v_vatom(i,2) += 0.5*v2;
        v_vatom(i,3) += 0.5*v3;
        v_vatom(i,4) += 0.5*v4;
        v_vatom(i,5) += 0.5*v5;
      }
      if (NEIGHFLAG != FULL && (NEWTON_PAIR || j < nlocal)) {
        v_vatom(j,0) += 0.5*v0;
        v_vatom(j,1) += 0.5*v1;
        v_vatom(j,2) += 0.5*v2;
        v_vatom(j,3) += 0.5*v3;
        v_vatom(j,4) += 0.5*v4;
        v_vatom(j,5) += 0.5*v5;
      }
    }
  }
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairLJCutDipoleCutKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}


namespace LAMMPS_NS {
template class PairLJCutDipoleCutKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutDipoleCutKokkos<LMPHostType>;
#endif
}

