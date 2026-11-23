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
------------------------------------------------------------------------- */

#include "pair_lj_sf_dipole_sf_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJSFDipoleSFKokkos<DeviceType>::PairLJSFDipoleSFKokkos(LAMMPS *lmp) : PairLJSFDipoleSF(lmp)
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
PairLJSFDipoleSFKokkos<DeviceType>::~PairLJSFDipoleSFKokkos()
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
void PairLJSFDipoleSFKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  mu = atomKK->k_mu.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  for (int m = 0; m < 4; ++m) {
    special_lj[m] = static_cast<KK_FLOAT>(force->special_lj[m]);
    special_coul[m] = static_cast<KK_FLOAT>(force->special_coul[m]);
  }
  qqrd2e = static_cast<KK_FLOAT>(force->qqrd2e);
  newton_pair = force->newton_pair;

  // get the neighbor list and neighbors used in operator()

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  if (neighflag == HALF) {
    if (newton_pair) launch<HALF,1>(ev);
    else launch<HALF,0>(ev);
  } else if (neighflag == HALFTHREAD) {
    if (newton_pair) launch<HALFTHREAD,1>(ev);
    else launch<HALFTHREAD,0>(ev);
  } else if (neighflag == FULL) {
    if (newton_pair) launch<FULL,1>(ev);
    else launch<FULL,0>(ev);
  }

  if (eflag_global) {
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
void PairLJSFDipoleSFKokkos<DeviceType>::launch(EV_FLOAT &ev)
{
  if (evflag)
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLJSFDipoleSFKernel<NEIGHFLAG,NEWTON_PAIR,1> >(0,inum),*this,ev);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLJSFDipoleSFKernel<NEIGHFLAG,NEWTON_PAIR,0> >(0,inum),*this);
}

/* ----------------------------------------------------------------------
   shifted-force Lennard-Jones plus shifted-force charge and point dipole
   interactions, same as PairLJSFDipoleSF::compute()
------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairLJSFDipoleSFKokkos<DeviceType>::operator()(TagPairLJSFDipoleSFKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int ii, EV_FLOAT &ev) const {

  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;

  const int i = d_ilist[ii];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);
  const KK_FLOAT qtmp = q(i);
  const KK_FLOAT mui = mu(i,3);
  const int itype = type(i);
  const int jnum = d_numneigh[i];

  KK_ACC_FLOAT fx_i = 0.0;
  KK_ACC_FLOAT fy_i = 0.0;
  KK_ACC_FLOAT fz_i = 0.0;
  KK_ACC_FLOAT torquex_i = 0.0;
  KK_ACC_FLOAT torquey_i = 0.0;
  KK_ACC_FLOAT torquez_i = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    const KK_FLOAT factor_lj = special_lj[sbmask(j)];
    const KK_FLOAT factor_coul = special_coul[sbmask(j)];
    j &= NEIGHMASK;

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const int jtype = type(j);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (rsq >= d_cutsq(itype,jtype)) continue;

    const params_lj_sf_dipole &p = params(itype,jtype);
    const KK_FLOAT qj = q(j);
    const KK_FLOAT muj = mu(j,3);
    const KK_FLOAT r2inv = static_cast<KK_FLOAT>(1.0)/rsq;
    const KK_FLOAT rinv = Kokkos::sqrt(r2inv);

    // atom can have both a charge and dipole
    // i,j = charge-charge, dipole-dipole, dipole-charge, or charge-dipole

    KK_FLOAT forcecoulx = 0.0, forcecouly = 0.0, forcecoulz = 0.0;
    KK_FLOAT tixcoul = 0.0, tiycoul = 0.0, tizcoul = 0.0;
    KK_FLOAT tjxcoul = 0.0, tjycoul = 0.0, tjzcoul = 0.0;
    KK_FLOAT ecoul = 0.0;
    KK_FLOAT evdwl = 0.0;

    if (rsq < p.cut_coulsq) {
      const KK_FLOAT rcutcoul2inv = static_cast<KK_FLOAT>(1.0)/p.cut_coulsq;
      KK_FLOAT edd = 0.0, edq = 0.0, eqd = 0.0;

      if ((qtmp != static_cast<KK_FLOAT>(0.0)) && (qj != static_cast<KK_FLOAT>(0.0))) {
        const KK_FLOAT pre1 = qtmp*qj*rinv*(r2inv-rcutcoul2inv);
        forcecoulx += pre1*delx;
        forcecouly += pre1*dely;
        forcecoulz += pre1*delz;
      }

      if ((mui > static_cast<KK_FLOAT>(0.0)) && (muj > static_cast<KK_FLOAT>(0.0))) {
        const KK_FLOAT r3inv = r2inv*rinv;
        const KK_FLOAT r5inv = r3inv*r2inv;

        const KK_FLOAT pdotp = mu(i,0)*mu(j,0) + mu(i,1)*mu(j,1) + mu(i,2)*mu(j,2);
        const KK_FLOAT pidotr = mu(i,0)*delx + mu(i,1)*dely + mu(i,2)*delz;
        const KK_FLOAT pjdotr = mu(j,0)*delx + mu(j,1)*dely + mu(j,2)*delz;

        const KK_FLOAT afac = static_cast<KK_FLOAT>(1.0) - rsq*rsq * rcutcoul2inv*rcutcoul2inv;
        const KK_FLOAT pre1 = afac * (pdotp - static_cast<KK_FLOAT>(3.0) * r2inv * pidotr * pjdotr);
        const KK_FLOAT aforcecoulx = pre1*delx;
        const KK_FLOAT aforcecouly = pre1*dely;
        const KK_FLOAT aforcecoulz = pre1*delz;

        const KK_FLOAT bfac = static_cast<KK_FLOAT>(1.0) -
          static_cast<KK_FLOAT>(4.0)*rsq*Kokkos::sqrt(rsq*rcutcoul2inv)*rcutcoul2inv +
          static_cast<KK_FLOAT>(3.0)*rsq*rsq*rcutcoul2inv*rcutcoul2inv;
        const KK_FLOAT presf = static_cast<KK_FLOAT>(2.0) * r2inv * pidotr * pjdotr;
        const KK_FLOAT bforcecoulx = bfac * (pjdotr*mu(i,0)+pidotr*mu(j,0)-presf*delx);
        const KK_FLOAT bforcecouly = bfac * (pjdotr*mu(i,1)+pidotr*mu(j,1)-presf*dely);
        const KK_FLOAT bforcecoulz = bfac * (pjdotr*mu(i,2)+pidotr*mu(j,2)-presf*delz);

        forcecoulx += static_cast<KK_FLOAT>(3.0) * r5inv * (aforcecoulx + bforcecoulx);
        forcecouly += static_cast<KK_FLOAT>(3.0) * r5inv * (aforcecouly + bforcecouly);
        forcecoulz += static_cast<KK_FLOAT>(3.0) * r5inv * (aforcecoulz + bforcecoulz);

        const KK_FLOAT pre2 = static_cast<KK_FLOAT>(3.0) * bfac * r5inv * pjdotr;
        const KK_FLOAT pre3 = static_cast<KK_FLOAT>(3.0) * bfac * r5inv * pidotr;
        const KK_FLOAT pre4 = -bfac * r3inv;

        const KK_FLOAT crossx = pre4 * (mu(i,1)*mu(j,2) - mu(i,2)*mu(j,1));
        const KK_FLOAT crossy = pre4 * (mu(i,2)*mu(j,0) - mu(i,0)*mu(j,2));
        const KK_FLOAT crossz = pre4 * (mu(i,0)*mu(j,1) - mu(i,1)*mu(j,0));

        tixcoul += crossx + pre2 * (mu(i,1)*delz - mu(i,2)*dely);
        tiycoul += crossy + pre2 * (mu(i,2)*delx - mu(i,0)*delz);
        tizcoul += crossz + pre2 * (mu(i,0)*dely - mu(i,1)*delx);
        tjxcoul += -crossx + pre3 * (mu(j,1)*delz - mu(j,2)*dely);
        tjycoul += -crossy + pre3 * (mu(j,2)*delx - mu(j,0)*delz);
        tjzcoul += -crossz + pre3 * (mu(j,0)*dely - mu(j,1)*delx);

        if (EVFLAG && eflag)
          edd = bfac * (r3inv*pdotp - static_cast<KK_FLOAT>(3.0)*r5inv*pidotr*pjdotr);
      }

      if ((mui > static_cast<KK_FLOAT>(0.0)) && (qj != static_cast<KK_FLOAT>(0.0))) {
        const KK_FLOAT r3inv = r2inv*rinv;
        const KK_FLOAT r5inv = r3inv*r2inv;
        const KK_FLOAT pidotr = mu(i,0)*delx + mu(i,1)*dely + mu(i,2)*delz;
        const KK_FLOAT pre1 = static_cast<KK_FLOAT>(3.0) * qj * r5inv * pidotr * (static_cast<KK_FLOAT>(1.0)-rsq*rcutcoul2inv);
        const KK_FLOAT pqfac = static_cast<KK_FLOAT>(1.0) - static_cast<KK_FLOAT>(3.0)*rsq*rcutcoul2inv +
          static_cast<KK_FLOAT>(2.0)*rsq*Kokkos::sqrt(rsq*rcutcoul2inv)*rcutcoul2inv;
        const KK_FLOAT pre2 = qj * r3inv * pqfac;

        forcecoulx += pre2*mu(i,0) - pre1*delx;
        forcecouly += pre2*mu(i,1) - pre1*dely;
        forcecoulz += pre2*mu(i,2) - pre1*delz;
        tixcoul += pre2 * (mu(i,1)*delz - mu(i,2)*dely);
        tiycoul += pre2 * (mu(i,2)*delx - mu(i,0)*delz);
        tizcoul += pre2 * (mu(i,0)*dely - mu(i,1)*delx);

        if (EVFLAG && eflag) edq = -qj * r3inv * pqfac * pidotr;
      }

      if ((muj > static_cast<KK_FLOAT>(0.0)) && (qtmp != static_cast<KK_FLOAT>(0.0))) {
        const KK_FLOAT r3inv = r2inv*rinv;
        const KK_FLOAT r5inv = r3inv*r2inv;
        const KK_FLOAT pjdotr = mu(j,0)*delx + mu(j,1)*dely + mu(j,2)*delz;
        const KK_FLOAT pre1 = static_cast<KK_FLOAT>(3.0) * qtmp * r5inv * pjdotr * (static_cast<KK_FLOAT>(1.0)-rsq*rcutcoul2inv);
        const KK_FLOAT qpfac = static_cast<KK_FLOAT>(1.0) - static_cast<KK_FLOAT>(3.0)*rsq*rcutcoul2inv +
          static_cast<KK_FLOAT>(2.0)*rsq*Kokkos::sqrt(rsq*rcutcoul2inv)*rcutcoul2inv;
        const KK_FLOAT pre2 = qtmp * r3inv * qpfac;

        forcecoulx += pre1*delx - pre2*mu(j,0);
        forcecouly += pre1*dely - pre2*mu(j,1);
        forcecoulz += pre1*delz - pre2*mu(j,2);
        tjxcoul += -pre2 * (mu(j,1)*delz - mu(j,2)*dely);
        tjycoul += -pre2 * (mu(j,2)*delx - mu(j,0)*delz);
        tjzcoul += -pre2 * (mu(j,0)*dely - mu(j,1)*delx);

        if (EVFLAG && eflag) eqd = qtmp * r3inv * qpfac * pjdotr;
      }

      // sum the energy terms in the same order as the base style

      if (EVFLAG && eflag) {
        ecoul = static_cast<KK_FLOAT>(1.0) - Kokkos::sqrt(rsq/p.cut_coulsq);
        ecoul *= ecoul;
        ecoul *= qtmp * qj * rinv;
        if ((mui > static_cast<KK_FLOAT>(0.0)) && (muj > static_cast<KK_FLOAT>(0.0))) ecoul += edd;
        if ((mui > static_cast<KK_FLOAT>(0.0)) && (qj != static_cast<KK_FLOAT>(0.0))) ecoul += edq;
        if ((muj > static_cast<KK_FLOAT>(0.0)) && (qtmp != static_cast<KK_FLOAT>(0.0))) ecoul += eqd;
        ecoul *= factor_coul*qqrd2e*p.scale;
      }
    }

    // LJ interaction

    KK_FLOAT forcelj = 0.0;
    if (rsq < p.cut_ljsq) {
      const KK_FLOAT r6inv = r2inv*r2inv*r2inv;
      const KK_FLOAT forceljcut = r6inv*(p.lj1*r6inv-p.lj2)*r2inv;

      const KK_FLOAT rcutlj2inv = static_cast<KK_FLOAT>(1.0) / p.cut_ljsq;
      const KK_FLOAT rcutlj6inv = rcutlj2inv * rcutlj2inv * rcutlj2inv;
      const KK_FLOAT forceljsf = (p.lj1*rcutlj6inv - p.lj2) * rcutlj6inv * rcutlj2inv;

      forcelj = factor_lj * (forceljcut - forceljsf);

      if (EVFLAG && eflag) {
        evdwl = r6inv*(p.lj3*r6inv-p.lj4) +
          rcutlj6inv*(static_cast<KK_FLOAT>(6.0)*p.lj3*rcutlj6inv-static_cast<KK_FLOAT>(3.0)*p.lj4)*rsq*rcutlj2inv +
          rcutlj6inv*(static_cast<KK_FLOAT>(-7.0)*p.lj3*rcutlj6inv+static_cast<KK_FLOAT>(4.0)*p.lj4);
        evdwl *= factor_lj;
      }
    }

    // total force

    const KK_FLOAT fq = factor_coul*qqrd2e*p.scale;
    const KK_FLOAT fx = fq*forcecoulx + delx*forcelj;
    const KK_FLOAT fy = fq*forcecouly + dely*forcelj;
    const KK_FLOAT fz = fq*forcecoulz + delz*forcelj;

    // force & torque accumulation

    fx_i += static_cast<KK_ACC_FLOAT>(fx);
    fy_i += static_cast<KK_ACC_FLOAT>(fy);
    fz_i += static_cast<KK_ACC_FLOAT>(fz);
    torquex_i += static_cast<KK_ACC_FLOAT>(fq*tixcoul);
    torquey_i += static_cast<KK_ACC_FLOAT>(fq*tiycoul);
    torquez_i += static_cast<KK_ACC_FLOAT>(fq*tizcoul);

    if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
      a_f(j,0) -= static_cast<KK_ACC_FLOAT>(fx);
      a_f(j,1) -= static_cast<KK_ACC_FLOAT>(fy);
      a_f(j,2) -= static_cast<KK_ACC_FLOAT>(fz);
      a_torque(j,0) += static_cast<KK_ACC_FLOAT>(fq*tjxcoul);
      a_torque(j,1) += static_cast<KK_ACC_FLOAT>(fq*tjycoul);
      a_torque(j,2) += static_cast<KK_ACC_FLOAT>(fq*tjzcoul);
    }

    if (EVFLAG) {
      if (eflag_global) {
        const KK_ACC_FLOAT efac = (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(j<nlocal)))?static_cast<KK_ACC_FLOAT>(1.0):static_cast<KK_ACC_FLOAT>(0.5));
        ev.evdwl += efac*static_cast<KK_ACC_FLOAT>(evdwl);
        ev.ecoul += efac*static_cast<KK_ACC_FLOAT>(ecoul);
      }
      if (eflag_atom || vflag_either)
        ev_tally_xyz<NEIGHFLAG, NEWTON_PAIR>(ev, i, j, ecoul+evdwl, fx, fy, fz, delx, dely, delz);
    }
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
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairLJSFDipoleSFKokkos<DeviceType>::operator()(TagPairLJSFDipoleSFKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairLJSFDipoleSFKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairLJSFDipoleSFKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT & ev, int i, int j, const KK_FLOAT &epair,
                                                        KK_FLOAT fx, KK_FLOAT fy, KK_FLOAT fz,
                                                        KK_FLOAT delx, KK_FLOAT dely, KK_FLOAT delz) const
{
  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = d_vatom;

  if (eflag_atom) {
    const KK_ACC_FLOAT epairhalf = static_cast<KK_ACC_FLOAT>(0.5) * static_cast<KK_ACC_FLOAT>(epair);
    if (NEIGHFLAG == FULL || newton_pair || i < nlocal)
      v_eatom[i] += epairhalf;
    if (NEIGHFLAG != FULL && (newton_pair || j < nlocal))
      v_eatom[j] += epairhalf;
  }

  if (vflag_either) {
    const KK_ACC_FLOAT v0 = static_cast<KK_ACC_FLOAT>(delx*fx);
    const KK_ACC_FLOAT v1 = static_cast<KK_ACC_FLOAT>(dely*fy);
    const KK_ACC_FLOAT v2 = static_cast<KK_ACC_FLOAT>(delz*fz);
    const KK_ACC_FLOAT v3 = static_cast<KK_ACC_FLOAT>(delx*fy);
    const KK_ACC_FLOAT v4 = static_cast<KK_ACC_FLOAT>(delx*fz);
    const KK_ACC_FLOAT v5 = static_cast<KK_ACC_FLOAT>(dely*fz);

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
            ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5)*v0;
            ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5)*v1;
            ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5)*v2;
            ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5)*v3;
            ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5)*v4;
            ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5)*v5;
          }
          if (j < nlocal) {
            ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5)*v0;
            ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5)*v1;
            ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5)*v2;
            ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5)*v3;
            ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5)*v4;
            ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5)*v5;
          }
        }
      } else { //neigh full
        ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5)*v0;
        ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5)*v1;
        ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5)*v2;
        ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5)*v3;
        ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5)*v4;
        ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5)*v5;
      }
    }

    if (vflag_atom) {

      if (NEIGHFLAG == FULL || NEWTON_PAIR || i < nlocal) {
        v_vatom(i,0) += static_cast<KK_ACC_FLOAT>(0.5)*v0;
        v_vatom(i,1) += static_cast<KK_ACC_FLOAT>(0.5)*v1;
        v_vatom(i,2) += static_cast<KK_ACC_FLOAT>(0.5)*v2;
        v_vatom(i,3) += static_cast<KK_ACC_FLOAT>(0.5)*v3;
        v_vatom(i,4) += static_cast<KK_ACC_FLOAT>(0.5)*v4;
        v_vatom(i,5) += static_cast<KK_ACC_FLOAT>(0.5)*v5;
      }
      if (NEIGHFLAG != FULL && (NEWTON_PAIR || j < nlocal)) {
        v_vatom(j,0) += static_cast<KK_ACC_FLOAT>(0.5)*v0;
        v_vatom(j,1) += static_cast<KK_ACC_FLOAT>(0.5)*v1;
        v_vatom(j,2) += static_cast<KK_ACC_FLOAT>(0.5)*v2;
        v_vatom(j,3) += static_cast<KK_ACC_FLOAT>(0.5)*v3;
        v_vatom(j,4) += static_cast<KK_ACC_FLOAT>(0.5)*v4;
        v_vatom(j,5) += static_cast<KK_ACC_FLOAT>(0.5)*v5;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJSFDipoleSFKokkos<DeviceType>::allocate()
{
  PairLJSFDipoleSF::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_lj_sf_dipole**,Kokkos::LayoutRight,DeviceType>("PairLJSFDipoleSF::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJSFDipoleSFKokkos<DeviceType>::init_style()
{
  PairLJSFDipoleSF::init_style();

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
double PairLJSFDipoleSFKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJSFDipoleSF::init_one(i,j);

  k_params.view_host()(i,j).cut_ljsq = static_cast<KK_FLOAT>(cut_ljsq[i][j]);
  k_params.view_host()(i,j).cut_coulsq = static_cast<KK_FLOAT>(cut_coulsq[i][j]);
  k_params.view_host()(i,j).lj1 = static_cast<KK_FLOAT>(lj1[i][j]);
  k_params.view_host()(i,j).lj2 = static_cast<KK_FLOAT>(lj2[i][j]);
  k_params.view_host()(i,j).lj3 = static_cast<KK_FLOAT>(lj3[i][j]);
  k_params.view_host()(i,j).lj4 = static_cast<KK_FLOAT>(lj4[i][j]);
  k_params.view_host()(i,j).scale = static_cast<KK_FLOAT>(scale[i][j]);
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  k_params.modify_host();

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();

  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
int PairLJSFDipoleSFKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairLJSFDipoleSFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJSFDipoleSFKokkos<LMPHostType>;
#endif
}
