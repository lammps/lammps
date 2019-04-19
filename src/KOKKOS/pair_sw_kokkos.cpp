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

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_sw_kokkos.h"
#include "kokkos.h"
#include "pair_kokkos.h"
#include "atom_kokkos.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairSWKokkos<DeviceType>::PairSWKokkos(LAMMPS *lmp) : PairSW(lmp)
{
  respa_enable = 0;


  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TAG_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairSWKokkos<DeviceType>::~PairSWKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = NULL;
    vatom = NULL;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairSWKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  nall = atom->nlocal + atom->nghost;

  const int inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  copymode = 1;

  EV_FLOAT ev;
  EV_FLOAT ev_all;

  // build short neighbor list

  int max_neighs = d_neighbors.extent(1);

  if ((d_neighbors_short.extent(1) != max_neighs) ||
     (d_neighbors_short.extent(0) != ignum)) {
    d_neighbors_short = Kokkos::View<int**,DeviceType>("SW::neighbors_short",ignum,max_neighs);
  }
  if (d_numneigh_short.extent(0)!=ignum)
    d_numneigh_short = Kokkos::View<int*,DeviceType>("SW::numneighs_short",ignum);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairSWComputeShortNeigh>(0,neighflag==FULL?ignum:inum), *this);

  // loop over neighbor list of my atoms

  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairSWComputeHalf<HALF,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairSWComputeHalf<HALF,0> >(0,inum),*this);
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairSWComputeHalf<HALFTHREAD,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairSWComputeHalf<HALFTHREAD,0> >(0,inum),*this);
    ev_all += ev;
  } else if (neighflag == FULL) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairSWComputeFullA<FULL,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairSWComputeFullA<FULL,0> >(0,inum),*this);
    ev_all += ev;

    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairSWComputeFullB<FULL,1> >(0,ignum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairSWComputeFullB<FULL,0> >(0,ignum),*this);
    ev_all += ev;
  }

  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev_all.evdwl;
  if (vflag_global) {
    virial[0] += ev_all.v[0];
    virial[1] += ev_all.v[1];
    virial[2] += ev_all.v[2];
    virial[3] += ev_all.v[3];
    virial[4] += ev_all.v[4];
    virial[5] += ev_all.v[5];
  }

  if (eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f            = decltype(dup_f)();
    dup_eatom        = decltype(dup_eatom)();
    dup_vatom        = decltype(dup_vatom)();
  }
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::operator()(TagPairSWComputeShortNeigh, const int& ii) const {
    const int i = d_ilist[ii];
    const X_FLOAT xtmp = x(i,0);
    const X_FLOAT ytmp = x(i,1);
    const X_FLOAT ztmp = x(i,2);

    const int jnum = d_numneigh[i];
    int inside = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      const X_FLOAT delx = xtmp - x(j,0);
      const X_FLOAT dely = ytmp - x(j,1);
      const X_FLOAT delz = ztmp - x(j,2);
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cutmax*cutmax) {
        d_neighbors_short(i,inside) = j;
        inside++;
      }
    }
    d_numneigh_short(i) = inside;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::operator()(TagPairSWComputeHalf<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  F_FLOAT delr1[3],delr2[3],fj[3],fk[3];
  F_FLOAT evdwl = 0.0;
  F_FLOAT fpair = 0.0;

  const int i = d_ilist[ii];
  const tagint itag = tag[i];
  const int itype = d_map[type[i]];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);

  // two-body interactions, skip half of them

  const int jnum = d_numneigh_short[i];

  F_FLOAT fxtmpi = 0.0;
  F_FLOAT fytmpi = 0.0;
  F_FLOAT fztmpi = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    const tagint jtag = tag[j];

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2) < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1) < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    const int jtype = d_map[type[j]];

    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    const int ijparam = d_elem2param(itype,jtype,jtype);
    if (rsq >= d_params[ijparam].cutsq) continue;

    twobody(d_params[ijparam],rsq,fpair,eflag,evdwl);

    fxtmpi += delx*fpair;
    fytmpi += dely*fpair;
    fztmpi += delz*fpair;
    a_f(j,0) -= delx*fpair;
    a_f(j,1) -= dely*fpair;
    a_f(j,2) -= delz*fpair;

    if (EVFLAG) {
      if (eflag) ev.evdwl += evdwl;
      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }

  const int jnumm1 = jnum - 1;

  for (int jj = 0; jj < jnumm1; jj++) {
    int j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    const int jtype = d_map[type[j]];
    const int ijparam = d_elem2param(itype,jtype,jtype);
    delr1[0] = x(j,0) - xtmp;
    delr1[1] = x(j,1) - ytmp;
    delr1[2] = x(j,2) - ztmp;
    const F_FLOAT rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
    if (rsq1 >= d_params[ijparam].cutsq) continue;

    F_FLOAT fxtmpj = 0.0;
    F_FLOAT fytmpj = 0.0;
    F_FLOAT fztmpj = 0.0;

    for (int kk = jj+1; kk < jnum; kk++) {
      int k = d_neighbors_short(i,kk);
      k &= NEIGHMASK;
      const int ktype = d_map[type[k]];
      const int ikparam = d_elem2param(itype,ktype,ktype);
      const int ijkparam = d_elem2param(itype,jtype,ktype);

      delr2[0] = x(k,0) - xtmp;
      delr2[1] = x(k,1) - ytmp;
      delr2[2] = x(k,2) - ztmp;
      const F_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

      if (rsq2 >= d_params[ikparam].cutsq) continue;

      threebody(d_params[ijparam],d_params[ikparam],d_params[ijkparam],
                rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

      fxtmpi -= fj[0] + fk[0];
      fytmpi -= fj[1] + fk[1];
      fztmpi -= fj[2] + fk[2];
      fxtmpj += fj[0];
      fytmpj += fj[1];
      fztmpj += fj[2];
      a_f(k,0) += fk[0];
      a_f(k,1) += fk[1];
      a_f(k,2) += fk[2];

      if (EVFLAG) {
        if (eflag) ev.evdwl += evdwl;
        if (vflag_either || eflag_atom) this->template ev_tally3<NEIGHFLAG>(ev,i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
    }

    a_f(j,0) += fxtmpj;
    a_f(j,1) += fytmpj;
    a_f(j,2) += fztmpj;
  }

  a_f(i,0) += fxtmpi;
  a_f(i,1) += fytmpi;
  a_f(i,2) += fztmpi;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::operator()(TagPairSWComputeHalf<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairSWComputeHalf<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::operator()(TagPairSWComputeFullA<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  F_FLOAT delr1[3],delr2[3],fj[3],fk[3];
  F_FLOAT evdwl = 0.0;
  F_FLOAT fpair = 0.0;

  const int i = d_ilist[ii];

  const tagint itag = tag[i];
  const int itype = d_map[type[i]];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);

  // two-body interactions

  const int jnum = d_numneigh_short[i];

  F_FLOAT fxtmpi = 0.0;
  F_FLOAT fytmpi = 0.0;
  F_FLOAT fztmpi = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    const tagint jtag = tag[j];

    const int jtype = d_map[type[j]];

    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    const int ijparam = d_elem2param(itype,jtype,jtype);

    if (rsq >= d_params[ijparam].cutsq) continue;

    twobody(d_params[ijparam],rsq,fpair,eflag,evdwl);

    fxtmpi += delx*fpair;
    fytmpi += dely*fpair;
    fztmpi += delz*fpair;

    if (EVFLAG) {
      if (eflag) ev.evdwl += 0.5*evdwl;
      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }

  const int jnumm1 = jnum - 1;

  for (int jj = 0; jj < jnumm1; jj++) {
    int j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    const int jtype = d_map[type[j]];
    const int ijparam = d_elem2param(itype,jtype,jtype);
    delr1[0] = x(j,0) - xtmp;
    delr1[1] = x(j,1) - ytmp;
    delr1[2] = x(j,2) - ztmp;
    const F_FLOAT rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

    if (rsq1 >= d_params[ijparam].cutsq) continue;

    for (int kk = jj+1; kk < jnum; kk++) {
      int k = d_neighbors_short(i,kk);
      k &= NEIGHMASK;
      const int ktype = d_map[type[k]];
      const int ikparam = d_elem2param(itype,ktype,ktype);
      const int ijkparam = d_elem2param(itype,jtype,ktype);

      delr2[0] = x(k,0) - xtmp;
      delr2[1] = x(k,1) - ytmp;
      delr2[2] = x(k,2) - ztmp;
      const F_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

      if (rsq2 >= d_params[ikparam].cutsq) continue;

      threebody(d_params[ijparam],d_params[ikparam],d_params[ijkparam],
                rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

      fxtmpi -= fj[0] + fk[0];
      fytmpi -= fj[1] + fk[1];
      fztmpi -= fj[2] + fk[2];

      if (EVFLAG) {
        if (eflag) ev.evdwl += evdwl;
        if (vflag_either || eflag_atom) this->template ev_tally3<NEIGHFLAG>(ev,i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
    }
  }

  f(i,0) += fxtmpi;
  f(i,1) += fytmpi;
  f(i,2) += fztmpi;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::operator()(TagPairSWComputeFullA<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairSWComputeFullA<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::operator()(TagPairSWComputeFullB<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  F_FLOAT delr1[3],delr2[3],fj[3],fk[3];
  F_FLOAT evdwl = 0.0;

  const int i = d_ilist[ii];

  const int itype = d_map[type[i]];
  const X_FLOAT xtmpi = x(i,0);
  const X_FLOAT ytmpi = x(i,1);
  const X_FLOAT ztmpi = x(i,2);

  const int jnum = d_numneigh_short[i];

  F_FLOAT fxtmpi = 0.0;
  F_FLOAT fytmpi = 0.0;
  F_FLOAT fztmpi = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    if (j >= nlocal) continue;
    const int jtype = d_map[type[j]];
    const int jiparam = d_elem2param(jtype,itype,itype);
    const X_FLOAT xtmpj = x(j,0);
    const X_FLOAT ytmpj = x(j,1);
    const X_FLOAT ztmpj = x(j,2);

    delr1[0] = xtmpi - xtmpj;
    delr1[1] = ytmpi - ytmpj;
    delr1[2] = ztmpi - ztmpj;
    const F_FLOAT rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

    if (rsq1 >= d_params[jiparam].cutsq) continue;

    const int j_jnum = d_numneigh_short[j];

    for (int kk = 0; kk < j_jnum; kk++) {
      int k = d_neighbors_short(j,kk);
      k &= NEIGHMASK;
      if (k == i) continue;
      const int ktype = d_map[type[k]];
      const int jkparam = d_elem2param(jtype,ktype,ktype);
      const int jikparam = d_elem2param(jtype,itype,ktype);

      delr2[0] = x(k,0) - xtmpj;
      delr2[1] = x(k,1) - ytmpj;
      delr2[2] = x(k,2) - ztmpj;
      const F_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

      if (rsq2 >= d_params[jkparam].cutsq) continue;

      if (vflag_atom)
        threebody(d_params[jiparam],d_params[jkparam],d_params[jikparam],
                  rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);
      else
        threebodyj(d_params[jiparam],d_params[jkparam],d_params[jikparam],
                  rsq1,rsq2,delr1,delr2,fj);

      fxtmpi += fj[0];
      fytmpi += fj[1];
      fztmpi += fj[2];

      if (EVFLAG)
        if (vflag_atom || eflag_atom) ev_tally3_atom(ev,i,evdwl,0.0,fj,fk,delr1,delr2);
    }
  }

  f(i,0) += fxtmpi;
  f(i,1) += fytmpi;
  f(i,2) += fztmpi;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::operator()(TagPairSWComputeFullB<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairSWComputeFullB<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSWKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairSW::coeff(narg,arg);

  // sync map

  int n = atom->ntypes;

  DAT::tdual_int_1d k_map = DAT::tdual_int_1d("pair:map",n+1);
  HAT::t_int_1d h_map = k_map.h_view;

  for (int i = 1; i <= n; i++)
    h_map[i] = map[i];

  k_map.template modify<LMPHostType>();
  k_map.template sync<DeviceType>();

  d_map = k_map.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairSWKokkos<DeviceType>::init_style()
{
  PairSW::init_style();

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;

  // always request a full neighbor list

  if (neighflag == FULL || neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    if (neighflag == FULL)
      neighbor->requests[irequest]->ghost = 1;
    else
      neighbor->requests[irequest]->ghost = 0;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with pair sw/kk");
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairSWKokkos<DeviceType>::setup_params()
{
  PairSW::setup_params();

  // sync elem2param and params

  tdual_int_3d k_elem2param = tdual_int_3d("pair:elem2param",nelements,nelements,nelements);
  t_host_int_3d h_elem2param = k_elem2param.h_view;

  tdual_param_1d k_params = tdual_param_1d("pair:params",nparams);
  t_host_param_1d h_params = k_params.h_view;

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      for (int k = 0; k < nelements; k++)
        h_elem2param(i,j,k) = elem2param[i][j][k];

  for (int m = 0; m < nparams; m++)
    h_params[m] = params[m];

  k_elem2param.template modify<LMPHostType>();
  k_elem2param.template sync<DeviceType>();
  k_params.template modify<LMPHostType>();
  k_params.template sync<DeviceType>();

  d_elem2param = k_elem2param.template view<DeviceType>();
  d_params = k_params.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::twobody(const Param& param, const F_FLOAT& rsq, F_FLOAT& fforce,
                     const int& eflag, F_FLOAT& eng) const
{
  F_FLOAT r,rinvsq,rp,rq,rainv,rainvsq,expsrainv;

  r = sqrt(rsq);
  rinvsq = 1.0/rsq;
  rp = pow(r,-param.powerp);
  rq = pow(r,-param.powerq);
  rainv = 1.0 / (r - param.cut);
  rainvsq = rainv*rainv*r;
  expsrainv = exp(param.sigma * rainv);
  fforce = (param.c1*rp - param.c2*rq +
            (param.c3*rp -param.c4*rq) * rainvsq) * expsrainv * rinvsq;
  if (eflag) eng = (param.c5*rp - param.c6*rq) * expsrainv;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::threebody(const Param& paramij, const Param& paramik, const Param& paramijk,
                       const F_FLOAT& rsq1, const F_FLOAT& rsq2,
                       F_FLOAT *delr1, F_FLOAT *delr2,
                       F_FLOAT *fj, F_FLOAT *fk, const int& eflag, F_FLOAT& eng) const
{
  F_FLOAT r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  F_FLOAT r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
  F_FLOAT rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2;
  F_FLOAT facang,facang12,csfacang,csfac1,csfac2;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  rainv1 = 1.0/(r1 - paramij.cut);
  gsrainv1 = paramij.sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(rsq2);
  rinvsq2 = 1.0/rsq2;
  rainv2 = 1.0/(r2 - paramik.cut);
  gsrainv2 = paramik.sigma_gamma * rainv2;
  gsrainvsq2 = gsrainv2*rainv2/r2;
  expgsrainv2 = exp(gsrainv2);

  rinv12 = 1.0/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
  delcs = cs - paramijk.costheta;
  delcssq = delcs*delcs;

  facexp = expgsrainv1*expgsrainv2;

  // facrad = sqrt(paramij.lambda_epsilon*paramik.lambda_epsilon) *
  //          facexp*delcssq;

  facrad = paramijk.lambda_epsilon * facexp*delcssq;
  frad1 = facrad*gsrainvsq1;
  frad2 = facrad*gsrainvsq2;
  facang = paramijk.lambda_epsilon2 * facexp*delcs;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;

  fj[0] = delr1[0]*(frad1+csfac1)-delr2[0]*facang12;
  fj[1] = delr1[1]*(frad1+csfac1)-delr2[1]*facang12;
  fj[2] = delr1[2]*(frad1+csfac1)-delr2[2]*facang12;

  csfac2 = rinvsq2*csfacang;

  fk[0] = delr2[0]*(frad2+csfac2)-delr1[0]*facang12;
  fk[1] = delr2[1]*(frad2+csfac2)-delr1[1]*facang12;
  fk[2] = delr2[2]*(frad2+csfac2)-delr1[2]*facang12;

  if (eflag) eng = facrad;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::threebodyj(const Param& paramij, const Param& paramik, const Param& paramijk,
                       const F_FLOAT& rsq1, const F_FLOAT& rsq2, F_FLOAT *delr1, F_FLOAT *delr2, F_FLOAT *fj) const
{
  F_FLOAT r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  F_FLOAT r2, rainv2, gsrainv2, expgsrainv2;
  F_FLOAT rinv12,cs,delcs,delcssq,facexp,facrad,frad1;
  F_FLOAT facang,facang12,csfacang,csfac1;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0/rsq1;
  rainv1 = 1.0/(r1 - paramij.cut);
  gsrainv1 = paramij.sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(rsq2);
  rainv2 = 1.0/(r2 - paramik.cut);
  gsrainv2 = paramik.sigma_gamma * rainv2;
  expgsrainv2 = exp(gsrainv2);

  rinv12 = 1.0/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
  delcs = cs - paramijk.costheta;
  delcssq = delcs*delcs;

  facexp = expgsrainv1*expgsrainv2;

  // facrad = sqrt(paramij.lambda_epsilon*paramik.lambda_epsilon) *
  //          facexp*delcssq;

  facrad = paramijk.lambda_epsilon * facexp*delcssq;
  frad1 = facrad*gsrainvsq1;
  facang = paramijk.lambda_epsilon2 * facexp*delcs;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;

  fj[0] = delr1[0]*(frad1+csfac1)-delr2[0]*facang12;
  fj[1] = delr1[1]*(frad1+csfac1)-delr2[1]*facang12;
  fj[2] = delr1[2]*(frad1+csfac1)-delr2[2]*facang12;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  auto v_vatom = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  if (eflag_atom) {
    const E_FLOAT epairhalf = 0.5 * epair;
    a_eatom[i] += epairhalf;
    if (NEIGHFLAG != FULL)
      a_eatom[j] += epairhalf;
  }

  if (VFLAG) {
    const E_FLOAT v0 = delx*delx*fpair;
    const E_FLOAT v1 = dely*dely*fpair;
    const E_FLOAT v2 = delz*delz*fpair;
    const E_FLOAT v3 = delx*dely*fpair;
    const E_FLOAT v4 = delx*delz*fpair;
    const E_FLOAT v5 = dely*delz*fpair;

    if (vflag_global) {
      if (NEIGHFLAG != FULL) {
        ev.v[0] += v0;
        ev.v[1] += v1;
        ev.v[2] += v2;
        ev.v[3] += v3;
        ev.v[4] += v4;
        ev.v[5] += v5;
      } else {
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
      }
    }

    if (vflag_atom) {
      a_vatom(i,0) += 0.5*v0;
      a_vatom(i,1) += 0.5*v1;
      a_vatom(i,2) += 0.5*v2;
      a_vatom(i,3) += 0.5*v3;
      a_vatom(i,4) += 0.5*v4;
      a_vatom(i,5) += 0.5*v5;

      if (NEIGHFLAG != FULL) {
        a_vatom(j,0) += 0.5*v0;
        a_vatom(j,1) += 0.5*v1;
        a_vatom(j,2) += 0.5*v2;
        a_vatom(j,3) += 0.5*v3;
        a_vatom(j,4) += 0.5*v4;
        a_vatom(j,5) += 0.5*v5;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW and hbond potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
          const F_FLOAT &evdwl, const F_FLOAT &ecoul,
                     F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drji, F_FLOAT *drki) const
{
  F_FLOAT epairthird,v[6];

  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  auto v_vatom = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  if (eflag_atom) {
    epairthird = THIRD * (evdwl + ecoul);
    a_eatom[i] += epairthird;
    if (NEIGHFLAG != FULL) {
      a_eatom[j] += epairthird;
      a_eatom[k] += epairthird;
    }
  }

  if (VFLAG) {
    v[0] = drji[0]*fj[0] + drki[0]*fk[0];
    v[1] = drji[1]*fj[1] + drki[1]*fk[1];
    v[2] = drji[2]*fj[2] + drki[2]*fk[2];
    v[3] = drji[0]*fj[1] + drki[0]*fk[1];
    v[4] = drji[0]*fj[2] + drki[0]*fk[2];
    v[5] = drji[1]*fj[2] + drki[1]*fk[2];

    if (vflag_global) {
      ev.v[0] += v[0];
      ev.v[1] += v[1];
      ev.v[2] += v[2];
      ev.v[3] += v[3];
      ev.v[4] += v[4];
      ev.v[5] += v[5];
    }

    if (vflag_atom) {
      a_vatom(i,0) += THIRD*v[0]; a_vatom(i,1) += THIRD*v[1];
      a_vatom(i,2) += THIRD*v[2]; a_vatom(i,3) += THIRD*v[3];
      a_vatom(i,4) += THIRD*v[4]; a_vatom(i,5) += THIRD*v[5];

      if (NEIGHFLAG != FULL) {
        a_vatom(j,0) += THIRD*v[0]; a_vatom(j,1) += THIRD*v[1];
        a_vatom(j,2) += THIRD*v[2]; a_vatom(j,3) += THIRD*v[3];
        a_vatom(j,4) += THIRD*v[4]; a_vatom(j,5) += THIRD*v[5];

        a_vatom(k,0) += THIRD*v[0]; a_vatom(k,1) += THIRD*v[1];
        a_vatom(k,2) += THIRD*v[2]; a_vatom(k,3) += THIRD*v[3];
        a_vatom(k,4) += THIRD*v[4]; a_vatom(k,5) += THIRD*v[5];
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally eng_vdwl and virial into global and per-atom accumulators
   called by SW and hbond potentials, newton_pair is always on
   virial = riFi + rjFj + rkFk = (rj-ri) Fj + (rk-ri) Fk = drji*fj + drki*fk
 ------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairSWKokkos<DeviceType>::ev_tally3_atom(EV_FLOAT &ev, const int &i,
          const F_FLOAT &evdwl, const F_FLOAT &ecoul,
                     F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drji, F_FLOAT *drki) const
{
  F_FLOAT epairthird,v[6];

  const int VFLAG = vflag_either;

  if (eflag_atom) {
    epairthird = THIRD * (evdwl + ecoul);
    d_eatom[i] += epairthird;
  }

  if (VFLAG) {
    v[0] = drji[0]*fj[0] + drki[0]*fk[0];
    v[1] = drji[1]*fj[1] + drki[1]*fk[1];
    v[2] = drji[2]*fj[2] + drki[2]*fk[2];
    v[3] = drji[0]*fj[1] + drki[0]*fk[1];
    v[4] = drji[0]*fj[2] + drki[0]*fk[2];
    v[5] = drji[1]*fj[2] + drki[1]*fk[2];

    if (vflag_atom) {
      d_vatom(i,0) += THIRD*v[0]; d_vatom(i,1) += THIRD*v[1];
      d_vatom(i,2) += THIRD*v[2]; d_vatom(i,3) += THIRD*v[3];
      d_vatom(i,4) += THIRD*v[4]; d_vatom(i,5) += THIRD*v[5];
    }
  }
}

namespace LAMMPS_NS {
template class PairSWKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairSWKokkos<LMPHostType>;
#endif
}

