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
   Contributing author: Anders Hafreager (UiO), andershaf@gmail.com
------------------------------------------------------------------------- */

#include "pair_vashishta_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairVashishtaKokkos<DeviceType>::PairVashishtaKokkos(LAMMPS *lmp) : PairVashishta(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TAG_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

template<class DeviceType>
PairVashishtaKokkos<DeviceType>::~PairVashishtaKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    eatom = nullptr;
    vatom = nullptr;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairVashishtaKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  nall = atom->nlocal + atom->nghost;

  inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

  copymode = 1;

  EV_FLOAT ev;
  EV_FLOAT ev_all;

  // build short neighbor list

  int max_neighs = d_neighbors.extent(1);

  if (((int)d_neighbors_short_2body.extent(1) < max_neighs) ||
     ((int)d_neighbors_short_2body.extent(0) < ignum)) {
    d_neighbors_short_2body = typename AT::t_int_2d_dl("Vashishta::neighbors_short_2body",ignum*1.2,max_neighs);
  }
  if ((int)d_numneigh_short_2body.extent(0) < ignum)
    d_numneigh_short_2body = typename AT::t_int_1d("Vashishta::numneighs_short_2body",ignum*1.2);

  if (((int)d_neighbors_short_3body.extent(1) < max_neighs) ||
     ((int)d_neighbors_short_3body.extent(0) < ignum)) {
    d_neighbors_short_3body = typename AT::t_int_2d_dl("Vashishta::neighbors_short_3body",ignum*1.2,max_neighs);
  }
  if ((int)d_numneigh_short_3body.extent(0) < ignum)
    d_numneigh_short_3body = typename AT::t_int_1d("Vashishta::numneighs_short_3body",ignum*1.2);

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairVashishtaComputeShortNeigh>(0,neighflag==FULL?ignum:inum), *this);



  // loop over neighbor list of my atoms

  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairVashishtaComputeHalf<HALF,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairVashishtaComputeHalf<HALF,0> >(0,inum),*this);
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairVashishtaComputeHalf<HALFTHREAD,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairVashishtaComputeHalf<HALFTHREAD,0> >(0,inum),*this);
    ev_all += ev;
  } else if (neighflag == FULL) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairVashishtaComputeFullA<FULL,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairVashishtaComputeFullA<FULL,0> >(0,inum),*this);
    ev_all += ev;

    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairVashishtaComputeFullB<FULL,1> >(0,ignum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairVashishtaComputeFullB<FULL,0> >(0,ignum),*this);
    ev_all += ev;
  }

  if (eflag_global) eng_vdwl += static_cast<double>(ev_all.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev_all.v[0]);
    virial[1] += static_cast<double>(ev_all.v[1]);
    virial[2] += static_cast<double>(ev_all.v[2]);
    virial[3] += static_cast<double>(ev_all.v[3]);
    virial[4] += static_cast<double>(ev_all.v[4]);
    virial[5] += static_cast<double>(ev_all.v[5]);
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
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::operator()(TagPairVashishtaComputeShortNeigh, const int& ii) const {
    const int i = d_ilist[ii];
    const int itype = d_map[type[i]];
    const KK_FLOAT xtmp = x(i,0);
    const KK_FLOAT ytmp = x(i,1);
    const KK_FLOAT ztmp = x(i,2);

    const int jnum = d_numneigh[i];
    int inside_2body = 0;
    int inside_3body = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const int jtype = d_map[type[j]];
      const int ijparam = d_elem3param(itype,jtype,jtype);

      const KK_FLOAT delx = xtmp - x(j,0);
      const KK_FLOAT dely = ytmp - x(j,1);
      const KK_FLOAT delz = ztmp - x(j,2);
      const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < static_cast<KK_FLOAT>(d_params[ijparam].cutsq)) {
        d_neighbors_short_2body(ii,inside_2body) = j;
        inside_2body++;
      }

      if (rsq < static_cast<KK_FLOAT>(d_params[ijparam].cutsq2)) {
        d_neighbors_short_3body(ii,inside_3body) = j;
        inside_3body++;
      }
    }
    d_numneigh_short_2body(ii) = inside_2body;
    d_numneigh_short_3body(ii) = inside_3body;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::operator()(TagPairVashishtaComputeHalf<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is atomic

  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;

  KK_FLOAT delr1[3],delr2[3];
  KK_ACC_FLOAT fj[3],fk[3];
  KK_FLOAT evdwl = 0.0;
  KK_FLOAT fpair = 0.0;

  const int i = d_ilist[ii];
  const tagint itag = tag[i];
  const int itype = d_map[type[i]];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  // two-body interactions, skip half of them

  const int jnum = d_numneigh_short_2body[ii];

  KK_ACC_FLOAT fxtmpi = 0.0;
  KK_ACC_FLOAT fytmpi = 0.0;
  KK_ACC_FLOAT fztmpi = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short_2body(ii,jj);
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

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    const int ijparam = d_elem3param(itype,jtype,jtype);

    twobody(d_params[ijparam],rsq,fpair,eflag,evdwl);

    fxtmpi += static_cast<KK_ACC_FLOAT>(delx*fpair);
    fytmpi += static_cast<KK_ACC_FLOAT>(dely*fpair);
    fztmpi += static_cast<KK_ACC_FLOAT>(delz*fpair);
    a_f(j,0) -= static_cast<KK_ACC_FLOAT>(delx*fpair);
    a_f(j,1) -= static_cast<KK_ACC_FLOAT>(dely*fpair);
    a_f(j,2) -= static_cast<KK_ACC_FLOAT>(delz*fpair);

    if (EVFLAG) {
      if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }

  const int jnumm1 = d_numneigh_short_3body[ii];

  for (int jj = 0; jj < jnumm1-1; jj++) {
    int j = d_neighbors_short_3body(ii,jj);
    j &= NEIGHMASK;
    const int jtype = d_map[type[j]];
    const int ijparam = d_elem3param(itype,jtype,jtype);
    delr1[0] = x(j,0) - xtmp;
    delr1[1] = x(j,1) - ytmp;
    delr1[2] = x(j,2) - ztmp;
    const KK_FLOAT rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

    KK_ACC_FLOAT fxtmpj = 0.0;
    KK_ACC_FLOAT fytmpj = 0.0;
    KK_ACC_FLOAT fztmpj = 0.0;

    for (int kk = jj+1; kk < jnumm1; kk++) {
      int k = d_neighbors_short_3body(ii,kk);
      k &= NEIGHMASK;
      const int ktype = d_map[type[k]];
      const int ikparam = d_elem3param(itype,ktype,ktype);
      const int ijkparam = d_elem3param(itype,jtype,ktype);

      delr2[0] = x(k,0) - xtmp;
      delr2[1] = x(k,1) - ytmp;
      delr2[2] = x(k,2) - ztmp;
      const KK_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

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
        if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
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
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::operator()(TagPairVashishtaComputeHalf<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairVashishtaComputeHalf<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::operator()(TagPairVashishtaComputeFullA<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  KK_FLOAT delr1[3],delr2[3];
  KK_ACC_FLOAT fj[3],fk[3];
  KK_FLOAT evdwl = 0.0;
  KK_FLOAT fpair = 0.0;

  const int i = d_ilist[ii];

  const int itype = d_map[type[i]];
  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  // two-body interactions

  const int jnum = d_numneigh_short_2body[ii];

  KK_ACC_FLOAT fxtmpi = 0.0;
  KK_ACC_FLOAT fytmpi = 0.0;
  KK_ACC_FLOAT fztmpi = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short_2body(ii,jj);
    j &= NEIGHMASK;

    const int jtype = d_map[type[j]];

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    const int ijparam = d_elem3param(itype,jtype,jtype);

    twobody(d_params[ijparam],rsq,fpair,eflag,evdwl);

    fxtmpi += static_cast<KK_ACC_FLOAT>(delx*fpair);
    fytmpi += static_cast<KK_ACC_FLOAT>(dely*fpair);
    fztmpi += static_cast<KK_ACC_FLOAT>(delz*fpair);

    if (EVFLAG) {
      if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(evdwl);
      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }

  const int jnumm1 = d_numneigh_short_3body[ii];

  for (int jj = 0; jj < jnumm1-1; jj++) {
    int j = d_neighbors_short_3body(ii,jj);
    j &= NEIGHMASK;
    const int jtype = d_map[type[j]];
    const int ijparam = d_elem3param(itype,jtype,jtype);
    delr1[0] = x(j,0) - xtmp;
    delr1[1] = x(j,1) - ytmp;
    delr1[2] = x(j,2) - ztmp;
    const KK_FLOAT rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

    for (int kk = jj+1; kk < jnumm1; kk++) {
      int k = d_neighbors_short_3body(ii,kk);
      k &= NEIGHMASK;
      const int ktype = d_map[type[k]];
      const int ikparam = d_elem3param(itype,ktype,ktype);
      const int ijkparam = d_elem3param(itype,jtype,ktype);

      delr2[0] = x(k,0) - xtmp;
      delr2[1] = x(k,1) - ytmp;
      delr2[2] = x(k,2) - ztmp;
      const KK_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

      threebody(d_params[ijparam],d_params[ikparam],d_params[ijkparam],
                rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

      fxtmpi -= fj[0] + fk[0];
      fytmpi -= fj[1] + fk[1];
      fztmpi -= fj[2] + fk[2];

      if (EVFLAG) {
        if (eflag) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
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
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::operator()(TagPairVashishtaComputeFullA<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairVashishtaComputeFullA<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::operator()(TagPairVashishtaComputeFullB<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  KK_FLOAT delr1[3],delr2[3];
  KK_ACC_FLOAT fj[3],fk[3];
  KK_FLOAT evdwl = 0.0;

  const int i = d_ilist[ii];

  const int itype = d_map[type[i]];
  const KK_FLOAT xtmpi = x(i,0);
  const KK_FLOAT ytmpi = x(i,1);
  const KK_FLOAT ztmpi = x(i,2);

  const int jnum = d_numneigh_short_3body[ii];

  KK_ACC_FLOAT fxtmpi = 0.0;
  KK_ACC_FLOAT fytmpi = 0.0;
  KK_ACC_FLOAT fztmpi = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short_3body(ii,jj);
    j &= NEIGHMASK;
    if (j >= nlocal) continue;
    const int jtype = d_map[type[j]];
    const int jiparam = d_elem3param(jtype,itype,itype);
    const KK_FLOAT xtmpj = x(j,0);
    const KK_FLOAT ytmpj = x(j,1);
    const KK_FLOAT ztmpj = x(j,2);

    delr1[0] = xtmpi - xtmpj;
    delr1[1] = ytmpi - ytmpj;
    delr1[2] = ztmpi - ztmpj;
    const KK_FLOAT rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

    const int j_jnum = d_numneigh_short_3body[j];

    for (int kk = 0; kk < j_jnum; kk++) {
      int k = d_neighbors_short_3body(j,kk);
      k &= NEIGHMASK;
      if (k == i) continue;
      const int ktype = d_map[type[k]];
      const int jkparam = d_elem3param(jtype,ktype,ktype);
      const int jikparam = d_elem3param(jtype,itype,ktype);

      delr2[0] = x(k,0) - xtmpj;
      delr2[1] = x(k,1) - ytmpj;
      delr2[2] = x(k,2) - ztmpj;
      const KK_FLOAT rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

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
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::operator()(TagPairVashishtaComputeFullB<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairVashishtaComputeFullB<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairVashishtaKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairVashishta::coeff(narg,arg);

  // sync map

  int n = atom->ntypes;

  DAT::tdual_int_1d k_map = DAT::tdual_int_1d("pair:map",n+1);
  HAT::t_int_1d h_map = k_map.view_host();

  for (int i = 1; i <= n; i++)
    h_map[i] = map[i];

  k_map.modify_host();
  k_map.template sync<DeviceType>();

  d_map = k_map.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairVashishtaKokkos<DeviceType>::init_style()
{
  PairVashishta::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  request->enable_full();
  if (neighflag == FULL) request->enable_ghost();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairVashishtaKokkos<DeviceType>::setup_params()
{
  PairVashishta::setup_params();

  // sync elem3param and params

  DAT::tdual_int_3d k_elem3param = DAT::tdual_int_3d("pair:elem3param",nelements,nelements,nelements);
  HAT::t_int_3d h_elem3param = k_elem3param.view_host();

  tdual_param_1d k_params = tdual_param_1d("pair:params",nparams);
  t_host_param_1d h_params = k_params.view_host();

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      for (int k = 0; k < nelements; k++)
        h_elem3param(i,j,k) = elem3param[i][j][k];

  for (int m = 0; m < nparams; m++)
    h_params[m] = params[m];

  k_elem3param.modify_host();
  k_elem3param.template sync<DeviceType>();
  k_params.modify_host();
  k_params.template sync<DeviceType>();

  d_elem3param = k_elem3param.template view<DeviceType>();
  d_params = k_params.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::twobody(const Param& param, const KK_FLOAT& rsq, KK_FLOAT& fforce,
                     const int& eflag, KK_FLOAT& eng) const
{
  KK_FLOAT r,rinvsq,r4inv,r6inv,reta,lam1r,lam4r,vc2,vc3;
  const KK_FLOAT eta_kk = static_cast<KK_FLOAT>(param.eta);
  const KK_FLOAT lam1inv_kk = static_cast<KK_FLOAT>(param.lam1inv);
  const KK_FLOAT lam4inv_kk = static_cast<KK_FLOAT>(param.lam4inv);
  const KK_FLOAT zizj_kk = static_cast<KK_FLOAT>(param.zizj);
  const KK_FLOAT mbigd_kk = static_cast<KK_FLOAT>(param.mbigd);
  const KK_FLOAT dvrc_kk = static_cast<KK_FLOAT>(param.dvrc);
  const KK_FLOAT big6w_kk = static_cast<KK_FLOAT>(param.big6w);
  const KK_FLOAT heta_kk = static_cast<KK_FLOAT>(param.heta);
  const KK_FLOAT bigh_kk = static_cast<KK_FLOAT>(param.bigh);
  const KK_FLOAT bigw_kk = static_cast<KK_FLOAT>(param.bigw);
  const KK_FLOAT c0_kk = static_cast<KK_FLOAT>(param.c0);
  r = Kokkos::sqrt(rsq);
  rinvsq = static_cast<KK_FLOAT>(1.0)/rsq;
  r4inv = rinvsq*rinvsq;
  r6inv = rinvsq*r4inv;
  reta = Kokkos::pow(r,-eta_kk);
  lam1r = r*lam1inv_kk;
  lam4r = r*lam4inv_kk;
  vc2 = zizj_kk * Kokkos::exp(-lam1r)/r;
  vc3 = mbigd_kk * r4inv*Kokkos::exp(-lam4r);

  fforce = (dvrc_kk*r
      - (static_cast<KK_FLOAT>(4.0)*vc3 + lam4r*vc3+big6w_kk*r6inv
         - heta_kk*reta - vc2 - lam1r*vc2)
      ) * rinvsq;

  if (eflag) eng = bigh_kk*reta + vc2 - vc3 - bigw_kk*r6inv - r*dvrc_kk + c0_kk;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::threebody(const Param& paramij, const Param& paramik, const Param& paramijk,
                       const KK_FLOAT& rsq1, const KK_FLOAT& rsq2,
                       KK_FLOAT *delr1, KK_FLOAT *delr2,
                       KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, const int& eflag, KK_FLOAT& eng) const
{
  KK_FLOAT r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  KK_FLOAT r2,rinvsq2,rainv2,gsrainv2,gsrainvsq2,expgsrainv2;
  KK_FLOAT rinv12,cs,delcs,delcssq,facexp,facrad,frad1,frad2,pcsinv,pcsinvsq,pcs;
  KK_FLOAT facang,facang12,csfacang,csfac1,csfac2;

  const KK_FLOAT r0ij_kk = static_cast<KK_FLOAT>(paramij.r0);
  const KK_FLOAT gammaij_kk = static_cast<KK_FLOAT>(paramij.gamma);
  const KK_FLOAT r0ik_kk = static_cast<KK_FLOAT>(paramik.r0);
  const KK_FLOAT gammaik_kk = static_cast<KK_FLOAT>(paramik.gamma);
  const KK_FLOAT costheta_kk = static_cast<KK_FLOAT>(paramijk.costheta);
  const KK_FLOAT bigc_kk = static_cast<KK_FLOAT>(paramijk.bigc);
  const KK_FLOAT bigb_kk = static_cast<KK_FLOAT>(paramijk.bigb);
  const KK_FLOAT big2b_kk = static_cast<KK_FLOAT>(paramijk.big2b);

  r1 = Kokkos::sqrt(rsq1);
  rinvsq1 = static_cast<KK_FLOAT>(1.0)/rsq1;
  rainv1 = static_cast<KK_FLOAT>(1.0)/(r1 - r0ij_kk);
  gsrainv1 = gammaij_kk * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = Kokkos::exp(gsrainv1);

  r2 = Kokkos::sqrt(rsq2);
  rinvsq2 = static_cast<KK_FLOAT>(1.0)/rsq2;
  rainv2 = static_cast<KK_FLOAT>(1.0)/(r2 - r0ik_kk);
  gsrainv2 = gammaik_kk * rainv2;
  gsrainvsq2 = gsrainv2*rainv2/r2;
  expgsrainv2 = Kokkos::exp(gsrainv2);

  rinv12 = static_cast<KK_FLOAT>(1.0)/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
  delcs = cs - costheta_kk;
  delcssq = delcs*delcs;
  pcsinv = bigc_kk*delcssq + static_cast<KK_FLOAT>(1.0);
  pcsinvsq = pcsinv*pcsinv;
  pcs = delcssq/pcsinv;

  facexp = expgsrainv1*expgsrainv2;

  facrad = bigb_kk * facexp * pcs;
  frad1 = facrad*gsrainvsq1;
  frad2 = facrad*gsrainvsq2;
  facang = big2b_kk * facexp * delcs/pcsinvsq;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;

  fj[0] = static_cast<KK_ACC_FLOAT>(delr1[0]*(frad1+csfac1)-delr2[0]*facang12);
  fj[1] = static_cast<KK_ACC_FLOAT>(delr1[1]*(frad1+csfac1)-delr2[1]*facang12);
  fj[2] = static_cast<KK_ACC_FLOAT>(delr1[2]*(frad1+csfac1)-delr2[2]*facang12);

  csfac2 = rinvsq2*csfacang;

  fk[0] = static_cast<KK_ACC_FLOAT>(delr2[0]*(frad2+csfac2)-delr1[0]*facang12);
  fk[1] = static_cast<KK_ACC_FLOAT>(delr2[1]*(frad2+csfac2)-delr1[1]*facang12);
  fk[2] = static_cast<KK_ACC_FLOAT>(delr2[2]*(frad2+csfac2)-delr1[2]*facang12);

  if (eflag) eng = facrad;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::threebodyj(const Param& paramij, const Param& paramik, const Param& paramijk,
                       const KK_FLOAT& rsq1, const KK_FLOAT& rsq2, KK_FLOAT *delr1, KK_FLOAT *delr2, KK_ACC_FLOAT *fj) const
{
  KK_FLOAT r1,rinvsq1,rainv1,gsrainv1,gsrainvsq1,expgsrainv1;
  KK_FLOAT r2,rainv2,gsrainv2,expgsrainv2;
  KK_FLOAT rinv12,cs,delcs,delcssq,facexp,facrad,frad1,pcsinv,pcsinvsq,pcs;
  KK_FLOAT facang,facang12,csfacang,csfac1;

  const KK_FLOAT r0ij_kk = static_cast<KK_FLOAT>(paramij.r0);
  const KK_FLOAT gammaij_kk = static_cast<KK_FLOAT>(paramij.gamma);
  const KK_FLOAT r0ik_kk = static_cast<KK_FLOAT>(paramik.r0);
  const KK_FLOAT gammaik_kk = static_cast<KK_FLOAT>(paramik.gamma);
  const KK_FLOAT costheta_kk = static_cast<KK_FLOAT>(paramijk.costheta);
  const KK_FLOAT bigc_kk = static_cast<KK_FLOAT>(paramijk.bigc);
  const KK_FLOAT bigb_kk = static_cast<KK_FLOAT>(paramijk.bigb);
  const KK_FLOAT big2b_kk = static_cast<KK_FLOAT>(paramijk.big2b);

  r1 = Kokkos::sqrt(rsq1);
  rinvsq1 = static_cast<KK_FLOAT>(1.0)/rsq1;
  rainv1 = static_cast<KK_FLOAT>(1.0)/(r1 - r0ij_kk);
  gsrainv1 = gammaij_kk * rainv1;
  gsrainvsq1 = gsrainv1*rainv1/r1;
  expgsrainv1 = Kokkos::exp(gsrainv1);

  r2 = Kokkos::sqrt(rsq2);
  rainv2 = static_cast<KK_FLOAT>(1.0)/(r2 - r0ik_kk);
  gsrainv2 = gammaik_kk * rainv2;
  expgsrainv2 = Kokkos::exp(gsrainv2);

  rinv12 = static_cast<KK_FLOAT>(1.0)/(r1*r2);
  cs = (delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * rinv12;
  delcs = cs - costheta_kk;
  delcssq = delcs*delcs;
  pcsinv = bigc_kk*delcssq + static_cast<KK_FLOAT>(1.0);
  pcsinvsq = pcsinv*pcsinv;
  pcs = delcssq/pcsinv;

  facexp = expgsrainv1*expgsrainv2;

  facrad = bigb_kk * facexp * pcs;
  frad1 = facrad*gsrainvsq1;
  facang = big2b_kk * facexp * delcs/pcsinvsq;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;

  fj[0] = static_cast<KK_ACC_FLOAT>(delr1[0]*(frad1+csfac1)-delr2[0]*facang12);
  fj[1] = static_cast<KK_ACC_FLOAT>(delr1[1]*(frad1+csfac1)-delr2[1]*facang12);
  fj[2] = static_cast<KK_ACC_FLOAT>(delr1[2]*(frad1+csfac1)-delr2[2]*facang12);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are atomic for half/thread neighbor list

  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = d_vatom;


  if (eflag_atom) {
    const KK_FLOAT epairhalf = static_cast<KK_FLOAT>(0.5) * epair;
    v_eatom[i] += static_cast<KK_ACC_FLOAT>(epairhalf);
    if (NEIGHFLAG != FULL)
      v_eatom[j] += static_cast<KK_ACC_FLOAT>(epairhalf);
  }

  if (VFLAG) {
    const KK_FLOAT v0 = delx*delx*fpair;
    const KK_FLOAT v1 = dely*dely*fpair;
    const KK_FLOAT v2 = delz*delz*fpair;
    const KK_FLOAT v3 = delx*dely*fpair;
    const KK_FLOAT v4 = delx*delz*fpair;
    const KK_FLOAT v5 = dely*delz*fpair;

    if (vflag_global) {
      if (NEIGHFLAG != FULL) {
        ev.v[0] += static_cast<KK_ACC_FLOAT>(v0);
        ev.v[1] += static_cast<KK_ACC_FLOAT>(v1);
        ev.v[2] += static_cast<KK_ACC_FLOAT>(v2);
        ev.v[3] += static_cast<KK_ACC_FLOAT>(v3);
        ev.v[4] += static_cast<KK_ACC_FLOAT>(v4);
        ev.v[5] += static_cast<KK_ACC_FLOAT>(v5);
      } else {
        ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v0);
        ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v1);
        ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v2);
        ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v3);
        ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v4);
        ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v5);
      }
    }

    if (vflag_atom) {
      v_vatom(i,0) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v0);
      v_vatom(i,1) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v1);
      v_vatom(i,2) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v2);
      v_vatom(i,3) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v3);
      v_vatom(i,4) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v4);
      v_vatom(i,5) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v5);

      if (NEIGHFLAG != FULL) {
        v_vatom(j,0) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v0);
        v_vatom(j,1) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v1);
        v_vatom(j,2) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v2);
        v_vatom(j,3) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v3);
        v_vatom(j,4) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v4);
        v_vatom(j,5) += static_cast<KK_ACC_FLOAT>(0.5)*static_cast<KK_ACC_FLOAT>(v5);
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
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
          const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                     KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const
{
  KK_FLOAT epairthird,v[6];

  const int VFLAG = vflag_either;

// The eatom and vatom arrays are atomic for half/thread neighbor list

  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = d_vatom;

  if (eflag_atom) {
    epairthird = static_cast<KK_FLOAT>(THIRD) * (evdwl + ecoul);
    v_eatom[i] += static_cast<KK_ACC_FLOAT>(epairthird);
    if (NEIGHFLAG != FULL) {
      v_eatom[j] += static_cast<KK_ACC_FLOAT>(epairthird);
      v_eatom[k] += static_cast<KK_ACC_FLOAT>(epairthird);
    }
  }

  if (VFLAG) {
    v[0] = drji[0]*static_cast<KK_FLOAT>(fj[0]) + drki[0]*static_cast<KK_FLOAT>(fk[0]);
    v[1] = drji[1]*static_cast<KK_FLOAT>(fj[1]) + drki[1]*static_cast<KK_FLOAT>(fk[1]);
    v[2] = drji[2]*static_cast<KK_FLOAT>(fj[2]) + drki[2]*static_cast<KK_FLOAT>(fk[2]);
    v[3] = drji[0]*static_cast<KK_FLOAT>(fj[1]) + drki[0]*static_cast<KK_FLOAT>(fk[1]);
    v[4] = drji[0]*static_cast<KK_FLOAT>(fj[2]) + drki[0]*static_cast<KK_FLOAT>(fk[2]);
    v[5] = drji[1]*static_cast<KK_FLOAT>(fj[2]) + drki[1]*static_cast<KK_FLOAT>(fk[2]);

    if (vflag_global) {
      ev.v[0] += static_cast<KK_ACC_FLOAT>(v[0]);
      ev.v[1] += static_cast<KK_ACC_FLOAT>(v[1]);
      ev.v[2] += static_cast<KK_ACC_FLOAT>(v[2]);
      ev.v[3] += static_cast<KK_ACC_FLOAT>(v[3]);
      ev.v[4] += static_cast<KK_ACC_FLOAT>(v[4]);
      ev.v[5] += static_cast<KK_ACC_FLOAT>(v[5]);
    }

    if (vflag_atom) {
      v_vatom(i,0) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[0]); v_vatom(i,1) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[1]);
      v_vatom(i,2) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[2]); v_vatom(i,3) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[3]);
      v_vatom(i,4) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[4]); v_vatom(i,5) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[5]);

      if (NEIGHFLAG != FULL) {
        v_vatom(j,0) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[0]); v_vatom(j,1) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[1]);
        v_vatom(j,2) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[2]); v_vatom(j,3) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[3]);
        v_vatom(j,4) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[4]); v_vatom(j,5) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[5]);

        v_vatom(k,0) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[0]); v_vatom(k,1) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[1]);
        v_vatom(k,2) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[2]); v_vatom(k,3) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[3]);
        v_vatom(k,4) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[4]); v_vatom(k,5) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[5]);
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
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void PairVashishtaKokkos<DeviceType>::ev_tally3_atom(EV_FLOAT & /*ev*/, const int &i,
          const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                     KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const
{
  KK_FLOAT epairthird,v[6];

  // only tally per-atom virial data here: the caller computes valid fj/fk data
  // only when vflag_atom is set, so vflag_either would be too broad

  const int VFLAG = vflag_atom;

  if (eflag_atom) {
    epairthird = static_cast<KK_FLOAT>(THIRD) * (evdwl + ecoul);
    d_eatom[i] += static_cast<KK_ACC_FLOAT>(epairthird);
  }

  if (VFLAG) {
    v[0] = drji[0]*static_cast<KK_FLOAT>(fj[0]) + drki[0]*static_cast<KK_FLOAT>(fk[0]);
    v[1] = drji[1]*static_cast<KK_FLOAT>(fj[1]) + drki[1]*static_cast<KK_FLOAT>(fk[1]);
    v[2] = drji[2]*static_cast<KK_FLOAT>(fj[2]) + drki[2]*static_cast<KK_FLOAT>(fk[2]);
    v[3] = drji[0]*static_cast<KK_FLOAT>(fj[1]) + drki[0]*static_cast<KK_FLOAT>(fk[1]);
    v[4] = drji[0]*static_cast<KK_FLOAT>(fj[2]) + drki[0]*static_cast<KK_FLOAT>(fk[2]);
    v[5] = drji[1]*static_cast<KK_FLOAT>(fj[2]) + drki[1]*static_cast<KK_FLOAT>(fk[2]);

    if (vflag_atom) {
      d_vatom(i,0) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[0]); d_vatom(i,1) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[1]);
      d_vatom(i,2) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[2]); d_vatom(i,3) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[3]);
      d_vatom(i,4) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[4]); d_vatom(i,5) += static_cast<KK_ACC_FLOAT>(THIRD)*static_cast<KK_ACC_FLOAT>(v[5]);
    }
  }
}

namespace LAMMPS_NS {
template class PairVashishtaKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairVashishtaKokkos<LMPHostType>;
#endif
}

