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
   Contributing author: Ray Shan (SNL) and Christian Trott (SNL)
------------------------------------------------------------------------- */

#include "pair_tersoff_kokkos.h"
#include <cmath>
#include "kokkos.h"
#include "atom_kokkos.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list_kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define KOKKOS_CUDA_MAX_THREADS 256
#define KOKKOS_CUDA_MIN_BLOCKS 8

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTersoffKokkos<DeviceType>::PairTersoffKokkos(LAMMPS *lmp) : PairTersoff(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairTersoffKokkos<DeviceType>::~PairTersoffKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTersoffKokkos<DeviceType>::allocate()
{
  PairTersoff::allocate();

  int n = atom->ntypes;

  k_params = Kokkos::DualView<params_ters***,Kokkos::LayoutRight,DeviceType>
          ("PairTersoff::paramskk",n+1,n+1,n+1);
  paramskk = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairTersoffKokkos<DeviceType>::init_style()
{
  PairTersoff::init_style();

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = std::is_same<DeviceType,LMPHostType>::value &&
    !std::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = std::is_same<DeviceType,LMPDeviceType>::value;

  if (neighflag == FULL)
    error->all(FLERR,"Cannot (yet) use full neighbor list style with tersoff/kk");

  if (neighflag == FULL || neighflag == HALF || neighflag == HALFTHREAD) {
  //if (neighflag == FULL || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    if (neighflag == FULL)
      neighbor->requests[irequest]->ghost = 1;
    else
      neighbor->requests[irequest]->ghost = 0;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with tersoff/kk");
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTersoffKokkos<DeviceType>::setup_params()
{
  PairTersoff::setup_params();

  int i,j,k,m;
  int n = atom->ntypes;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      for (k = 1; k <= n; k++) {
        m = elem2param[map[i]][map[j]][map[k]];
        k_params.h_view(i,j,k).powerm = params[m].powerm;
        k_params.h_view(i,j,k).gamma = params[m].gamma;
        k_params.h_view(i,j,k).lam3 = params[m].lam3;
        k_params.h_view(i,j,k).c = params[m].c;
        k_params.h_view(i,j,k).d = params[m].d;
        k_params.h_view(i,j,k).h = params[m].h;
        k_params.h_view(i,j,k).powern = params[m].powern;
        k_params.h_view(i,j,k).beta = params[m].beta;
        k_params.h_view(i,j,k).lam2 = params[m].lam2;
        k_params.h_view(i,j,k).bigb = params[m].bigb;
        k_params.h_view(i,j,k).bigr = params[m].bigr;
        k_params.h_view(i,j,k).bigd = params[m].bigd;
        k_params.h_view(i,j,k).lam1 = params[m].lam1;
        k_params.h_view(i,j,k).biga = params[m].biga;
        k_params.h_view(i,j,k).cutsq = params[m].cutsq;
        k_params.h_view(i,j,k).c1 = params[m].c1;
        k_params.h_view(i,j,k).c2 = params[m].c2;
        k_params.h_view(i,j,k).c3 = params[m].c3;
        k_params.h_view(i,j,k).c4 = params[m].c4;
      }

  k_params.template modify<LMPHostType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairTersoffKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;

  inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

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
    d_neighbors_short = Kokkos::View<int**,DeviceType>("Tersoff::neighbors_short",ignum,max_neighs);
  }
  if (d_numneigh_short.extent(0)!=ignum)
    d_numneigh_short = Kokkos::View<int*,DeviceType>("Tersoff::numneighs_short",ignum);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagPairTersoffComputeShortNeigh>(0,neighflag==FULL?ignum:inum), *this);

  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairTersoffComputeHalf<HALF,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairTersoffComputeHalf<HALF,0> >(0,inum),*this);
    ev_all += ev;
  } else if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairTersoffComputeHalf<HALFTHREAD,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairTersoffComputeHalf<HALFTHREAD,0> >(0,inum),*this);
    ev_all += ev;
  } else if (neighflag == FULL) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairTersoffComputeFullA<FULL,1> >(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairTersoffComputeFullA<FULL,0> >(0,inum),*this);
    ev_all += ev;

    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairTersoffComputeFullB<FULL,1> >(0,ignum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairTersoffComputeFullB<FULL,0> >(0,ignum),*this);
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
    dup_f     = decltype(dup_f)();
    dup_eatom = decltype(dup_eatom)();
    dup_vatom = decltype(dup_vatom)();
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::operator()(TagPairTersoffComputeShortNeigh, const int& ii) const {
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
void PairTersoffKokkos<DeviceType>::operator()(TagPairTersoffComputeHalf<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // The f array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  const int i = d_ilist[ii];
  if (i >= nlocal) return;
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);
  const tagint itag = tag(i);

  F_FLOAT fi[3], fj[3], fk[3];

  //const AtomNeighborsConst d_neighbors_i = k_list.get_neighbors_const(i);
  const int jnum = d_numneigh_short[i];

  // repulsive

  F_FLOAT f_x = 0.0;
  F_FLOAT f_y = 0.0;
  F_FLOAT f_z = 0.0;

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    const int jtype = type(j);
    const tagint jtag = tag(j);

    if (itag > jtag) {
      if ((itag+jtag) % 2 == 0) continue;
    } else if (itag < jtag) {
      if ((itag+jtag) % 2 == 1) continue;
    } else {
      if (x(j,2)  < ztmp) continue;
      if (x(j,2) == ztmp && x(j,1)  < ytmp) continue;
      if (x(j,2) == ztmp && x(j,1) == ytmp && x(j,0) < xtmp) continue;
    }

    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const F_FLOAT cutsq = paramskk(itype,jtype,jtype).cutsq;

    if (rsq > cutsq) continue;

    const F_FLOAT r = sqrt(rsq);
    const F_FLOAT tmp_fce = ters_fc_k(itype,jtype,jtype,r);
    const F_FLOAT tmp_fcd = ters_dfc(itype,jtype,jtype,r);
    const F_FLOAT tmp_exp = exp(-paramskk(itype,jtype,jtype).lam1 * r);
    const F_FLOAT frep = -paramskk(itype,jtype,jtype).biga * tmp_exp *
                          (tmp_fcd - tmp_fce*paramskk(itype,jtype,jtype).lam1) / r;
    const F_FLOAT eng = tmp_fce * paramskk(itype,jtype,jtype).biga * tmp_exp;

    f_x += delx*frep;
    f_y += dely*frep;
    f_z += delz*frep;
    a_f(j,0) -= delx*frep;
    a_f(j,1) -= dely*frep;
    a_f(j,2) -= delz*frep;

    if (EVFLAG) {
      if (eflag) ev.evdwl += eng;
      if (vflag_either || eflag_atom) this->template ev_tally<NEIGHFLAG>(ev,i,j,eng,frep,delx,dely,delz);
    }
  }

  // attractive: bond order

  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    const int jtype = type(j);

    const F_FLOAT delx1 = xtmp - x(j,0);
    const F_FLOAT dely1 = ytmp - x(j,1);
    const F_FLOAT delz1 = ztmp - x(j,2);
    const F_FLOAT rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    const F_FLOAT cutsq1 = paramskk(itype,jtype,jtype).cutsq;

    F_FLOAT bo_ij = 0.0;
    if (rsq1 > cutsq1) continue;
    const F_FLOAT rij = sqrt(rsq1);

    for (int kk = 0; kk < jnum; kk++) {
      if (jj == kk) continue;
      int k = d_neighbors_short(i,kk);
      k &= NEIGHMASK;
      const int ktype = type(k);

      const F_FLOAT delx2 = xtmp - x(k,0);
      const F_FLOAT dely2 = ytmp - x(k,1);
      const F_FLOAT delz2 = ztmp - x(k,2);
      const F_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      const F_FLOAT cutsq2 = paramskk(itype,jtype,ktype).cutsq;

      if (rsq2 > cutsq2) continue;
      const F_FLOAT rik = sqrt(rsq2);
      bo_ij += bondorder(itype,jtype,ktype,rij,delx1,dely1,delz1,rik,delx2,dely2,delz2);
    }

    // attractive: pairwise potential and force

    const F_FLOAT fa = ters_fa_k(itype,jtype,jtype,rij);
    const F_FLOAT dfa = ters_dfa(itype,jtype,jtype,rij);
    const F_FLOAT bij = ters_bij_k(itype,jtype,jtype,bo_ij);
    const F_FLOAT fatt = -0.5*bij * dfa / rij;
    const F_FLOAT prefactor = 0.5*fa * ters_dbij(itype,jtype,jtype,bo_ij);

    f_x += delx1*fatt;
    f_y += dely1*fatt;
    f_z += delz1*fatt;
    F_FLOAT fj_x = -delx1*fatt;
    F_FLOAT fj_y = -dely1*fatt;
    F_FLOAT fj_z = -delz1*fatt;

    if (EVFLAG) {
      const F_FLOAT eng = 0.5*bij * fa;
      if (eflag) ev.evdwl += eng;
      if (vflag_either || eflag_atom)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,eng,fatt,delx1,dely1,delz1);
    }

    // attractive: three-body force

    for (int kk = 0; kk < jnum; kk++) {
      if (jj == kk) continue;
      int k = d_neighbors_short(i,kk);
      k &= NEIGHMASK;
      const int ktype = type(k);

      const F_FLOAT delx2 = xtmp - x(k,0);
      const F_FLOAT dely2 = ytmp - x(k,1);
      const F_FLOAT delz2 = ztmp - x(k,2);
      const F_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      const F_FLOAT cutsq2 = paramskk(itype,jtype,ktype).cutsq;

      if (rsq2 > cutsq2) continue;
      const F_FLOAT rik = sqrt(rsq2);
      ters_dthb(itype,jtype,ktype,prefactor,rij,delx1,dely1,delz1,
                rik,delx2,dely2,delz2,fi,fj,fk);

      f_x += fi[0];
      f_y += fi[1];
      f_z += fi[2];
      fj_x += fj[0];
      fj_y += fj[1];
      fj_z += fj[2];
      a_f(k,0) += fk[0];
      a_f(k,1) += fk[1];
      a_f(k,2) += fk[2];

      if (vflag_atom) {
        F_FLOAT delrij[3], delrik[3];
        delrij[0] = -delx1; delrij[1] = -dely1; delrij[2] = -delz1;
        delrik[0] = -delx2; delrik[1] = -dely2; delrik[2] = -delz2;
        if (vflag_either) this->template v_tally3<NEIGHFLAG>(ev,i,j,k,fj,fk,delrij,delrik);
      }
    }
    a_f(j,0) += fj_x;
    a_f(j,1) += fj_y;
    a_f(j,2) += fj_z;
  }
  a_f(i,0) += f_x;
  a_f(i,1) += f_y;
  a_f(i,2) += f_z;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::operator()(TagPairTersoffComputeHalf<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairTersoffComputeHalf<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::operator()(TagPairTersoffComputeFullA<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  int j,k,jj,kk,jtype,ktype;
  F_FLOAT rsq1, cutsq1, rsq2, cutsq2, rij, rik, bo_ij;
  F_FLOAT fi[3], fj[3], fk[3];
  X_FLOAT delx1, dely1, delz1, delx2, dely2, delz2;

  //const AtomNeighborsConst d_neighbors_i = k_list.get_neighbors_const(i);
  const int jnum = d_numneigh_short[i];

  // repulsive

  F_FLOAT f_x = 0.0;
  F_FLOAT f_y = 0.0;
  F_FLOAT f_z = 0.0;
  for (jj = 0; jj < jnum; jj++) {
    j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    const int jtype = type(j);

    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const F_FLOAT cutsq = paramskk(itype,jtype,jtype).cutsq;

    if (rsq > cutsq) continue;

    const F_FLOAT r = sqrt(rsq);
    const F_FLOAT tmp_fce = ters_fc_k(itype,jtype,jtype,r);
    const F_FLOAT tmp_fcd = ters_dfc(itype,jtype,jtype,r);
    const F_FLOAT tmp_exp = exp(-paramskk(itype,jtype,jtype).lam1 * r);
    const F_FLOAT frep = -paramskk(itype,jtype,jtype).biga * tmp_exp *
                          (tmp_fcd - tmp_fce*paramskk(itype,jtype,jtype).lam1) / r;
    const F_FLOAT eng = tmp_fce * paramskk(itype,jtype,jtype).biga * tmp_exp;

    f_x += delx*frep;
    f_y += dely*frep;
    f_z += delz*frep;

    if (EVFLAG) {
      if (eflag)
        ev.evdwl += 0.5*eng;
      if (vflag_either || eflag_atom)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,eng,frep,delx,dely,delz);
    }
  }

  // attractive: bond order

  for (jj = 0; jj < jnum; jj++) {
    j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    jtype = type(j);

    delx1 = xtmp - x(j,0);
    dely1 = ytmp - x(j,1);
    delz1 = ztmp - x(j,2);
    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    cutsq1 = paramskk(itype,jtype,jtype).cutsq;

    bo_ij = 0.0;
    if (rsq1 > cutsq1) continue;
    rij = sqrt(rsq1);

    for (kk = 0; kk < jnum; kk++) {
      if (jj == kk) continue;
      k = d_neighbors_short(i,kk);
      k &= NEIGHMASK;
      ktype = type(k);

      delx2 = xtmp - x(k,0);
      dely2 = ytmp - x(k,1);
      delz2 = ztmp - x(k,2);
      rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      cutsq2 = paramskk(itype,jtype,ktype).cutsq;

      if (rsq2 > cutsq2) continue;
      rik = sqrt(rsq2);
      bo_ij += bondorder(itype,jtype,ktype,rij,delx1,dely1,delz1,rik,delx2,dely2,delz2);
    }

    // attractive: pairwise potential and force

    const F_FLOAT fa = ters_fa_k(itype,jtype,jtype,rij);
    const F_FLOAT dfa = ters_dfa(itype,jtype,jtype,rij);
    const F_FLOAT bij = ters_bij_k(itype,jtype,jtype,bo_ij);
    const F_FLOAT fatt = -0.5*bij * dfa / rij;
    const F_FLOAT prefactor = 0.5*fa * ters_dbij(itype,jtype,jtype,bo_ij);
    const F_FLOAT eng = 0.5*bij * fa;

    f_x += delx1*fatt;
    f_y += dely1*fatt;
    f_z += delz1*fatt;

    if (EVFLAG) {
      if (eflag) ev.evdwl += 0.5*eng;
      if (vflag_either || eflag_atom)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,eng,fatt,delx1,dely1,delz1);
    }

    // attractive: three-body force

    for (kk = 0; kk < jnum; kk++) {
      if (jj == kk) continue;
      k = d_neighbors_short(i,kk);
      k &= NEIGHMASK;
      ktype = type(k);

      delx2 = xtmp - x(k,0);
      dely2 = ytmp - x(k,1);
      delz2 = ztmp - x(k,2);
      rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      cutsq2 = paramskk(itype,jtype,ktype).cutsq;

      if (rsq2 > cutsq2) continue;
      rik = sqrt(rsq2);
      ters_dthb(itype,jtype,ktype,prefactor,rij,delx1,dely1,delz1,
                rik,delx2,dely2,delz2,fi,fj,fk);

      f_x += fi[0];
      f_y += fi[1];
      f_z += fi[2];

      if (vflag_atom) {
        F_FLOAT delrij[3], delrik[3];
        delrij[0] = -delx1; delrij[1] = -dely1; delrij[2] = -delz1;
        delrik[0] = -delx2; delrik[1] = -dely2; delrik[2] = -delz2;
        if (vflag_either) this->template v_tally3<NEIGHFLAG>(ev,i,j,k,fj,fk,delrij,delrik);
      }
    }
  }
  f(i,0) += f_x;
  f(i,1) += f_y;
  f(i,2) += f_z;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::operator()(TagPairTersoffComputeFullA<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairTersoffComputeFullA<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::operator()(TagPairTersoffComputeFullB<NEIGHFLAG,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type(i);

  int j,k,jj,kk,jtype,ktype,j_jnum;
  F_FLOAT rsq1, cutsq1, rsq2, cutsq2, rij, rik, bo_ij;
  F_FLOAT fj[3], fk[3];
  X_FLOAT delx1, dely1, delz1, delx2, dely2, delz2;

  const int jnum = d_numneigh_short[i];

  F_FLOAT f_x = 0.0;
  F_FLOAT f_y = 0.0;
  F_FLOAT f_z = 0.0;

  // attractive: bond order

  for (jj = 0; jj < jnum; jj++) {
    j = d_neighbors_short(i,jj);
    j &= NEIGHMASK;
    if (j >= nlocal) continue;
    jtype = type(j);

    delx1 = x(j,0) - xtmp;
    dely1 = x(j,1) - ytmp;
    delz1 = x(j,2) - ztmp;
    rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
    cutsq1 = paramskk(jtype,itype,itype).cutsq;

    bo_ij = 0.0;
    if (rsq1 > cutsq1) continue;
    rij = sqrt(rsq1);

    j_jnum = d_numneigh_short[j];

    for (kk = 0; kk < j_jnum; kk++) {
      k = d_neighbors_short(j,kk);
      if (k == i) continue;
      k &= NEIGHMASK;
      ktype = type(k);

      delx2 = x(j,0) - x(k,0);
      dely2 = x(j,1) - x(k,1);
      delz2 = x(j,2) - x(k,2);
      rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      cutsq2 = paramskk(jtype,itype,ktype).cutsq;

      if (rsq2 > cutsq2) continue;
      rik = sqrt(rsq2);
      bo_ij += bondorder(jtype,itype,ktype,rij,delx1,dely1,delz1,rik,delx2,dely2,delz2);

    }

    // attractive: pairwise potential and force

    const F_FLOAT fa = ters_fa_k(jtype,itype,itype,rij);
    const F_FLOAT dfa = ters_dfa(jtype,itype,itype,rij);
    const F_FLOAT bij = ters_bij_k(jtype,itype,itype,bo_ij);
    const F_FLOAT fatt = -0.5*bij * dfa / rij;
    const F_FLOAT prefactor = 0.5*fa * ters_dbij(jtype,itype,itype,bo_ij);
    const F_FLOAT eng = 0.5*bij * fa;

    f_x -= delx1*fatt;
    f_y -= dely1*fatt;
    f_z -= delz1*fatt;

    if (EVFLAG) {
      if (eflag)
        ev.evdwl += 0.5 * eng;
      if (vflag_either || eflag_atom)
        this->template ev_tally<NEIGHFLAG>(ev,i,j,eng,fatt,delx1,dely1,delz1);
    }

    // attractive: three-body force

    for (kk = 0; kk < j_jnum; kk++) {
      k = d_neighbors_short(j,kk);
      if (k == i) continue;
      k &= NEIGHMASK;
      ktype = type(k);

      delx2 = x(j,0) - x(k,0);
      dely2 = x(j,1) - x(k,1);
      delz2 = x(j,2) - x(k,2);
      rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
      cutsq2 = paramskk(jtype,itype,ktype).cutsq;

      if (rsq2 > cutsq2) continue;
      rik = sqrt(rsq2);
      ters_dthbj(jtype,itype,ktype,prefactor,rij,delx1,dely1,delz1,
                rik,delx2,dely2,delz2,fj,fk);
      f_x += fj[0];
      f_y += fj[1];
      f_z += fj[2];

      if (vflag_atom) {
        F_FLOAT delrji[3], delrjk[3];
        delrji[0] = -delx1; delrji[1] = -dely1; delrji[2] = -delz1;
        delrjk[0] = -delx2; delrjk[1] = -dely2; delrjk[2] = -delz2;
        if (vflag_either) v_tally3_atom(ev,i,j,k,fj,fk,delrji,delrjk);
      }

      const F_FLOAT fa_jk = ters_fa_k(jtype,ktype,itype,rik);
      const F_FLOAT prefactor_jk = 0.5*fa_jk * ters_dbij(jtype,ktype,itype,bo_ij);
      ters_dthbk(jtype,ktype,itype,prefactor_jk,rik,delx2,dely2,delz2,
                rij,delx1,dely1,delz1,fk);
      f_x += fk[0];
      f_y += fk[1];
      f_z += fk[2];
    }
  }
  f(i,0) += f_x;
  f(i,1) += f_y;
  f(i,2) += f_z;
}

template<class DeviceType>
template<int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::operator()(TagPairTersoffComputeFullB<NEIGHFLAG,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,EVFLAG>(TagPairTersoffComputeFullB<NEIGHFLAG,EVFLAG>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffKokkos<DeviceType>::ters_fc_k(const int &i, const int &j,
                const int &k, const F_FLOAT &r) const
{
  const F_FLOAT ters_R = paramskk(i,j,k).bigr;
  const F_FLOAT ters_D = paramskk(i,j,k).bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffKokkos<DeviceType>::ters_dfc(const int &i, const int &j,
                const int &k, const F_FLOAT &r) const
{
  const F_FLOAT ters_R = paramskk(i,j,k).bigr;
  const F_FLOAT ters_D = paramskk(i,j,k).bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4/ters_D) * cos(MY_PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffKokkos<DeviceType>::bondorder(const int &i, const int &j, const int &k,
        const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
        const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2) const
{
  F_FLOAT arg, ex_delr;

  const F_FLOAT costheta = (dx1*dx2 + dy1*dy2 + dz1*dz2)/(rij*rik);

  const F_FLOAT param = paramskk(i,j,k).lam3 * (rij-rik);
  if (int(paramskk(i,j,k).powerm) == 3) arg = param*param*param;//pow(paramskk(i,j,k).lam3 * (rij-rik),3.0);
  else arg = param;

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc_k(i,j,k,rik) * ters_gijk(i,j,k,costheta) * ex_delr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffKokkos<DeviceType>::
        ters_gijk(const int &i, const int &j, const int &k, const F_FLOAT &cos) const
{
  const F_FLOAT ters_c = paramskk(i,j,k).c * paramskk(i,j,k).c;
  const F_FLOAT ters_d = paramskk(i,j,k).d * paramskk(i,j,k).d;
  const F_FLOAT hcth = paramskk(i,j,k).h - cos;

  return paramskk(i,j,k).gamma*(1.0 + ters_c/ters_d - ters_c/(ters_d+hcth*hcth));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffKokkos<DeviceType>::
        ters_dgijk(const int &i, const int &j, const int &k, const F_FLOAT &cos) const
{
  const F_FLOAT ters_c = paramskk(i,j,k).c * paramskk(i,j,k).c;
  const F_FLOAT ters_d = paramskk(i,j,k).d * paramskk(i,j,k).d;
  const F_FLOAT hcth = paramskk(i,j,k).h - cos;
  const F_FLOAT numerator = -2.0 * ters_c * hcth;
  const F_FLOAT denominator = 1.0/(ters_d + hcth*hcth);
  return paramskk(i,j,k).gamma * numerator * denominator * denominator;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffKokkos<DeviceType>::ters_fa_k(const int &i, const int &j,
                const int &k, const F_FLOAT &r) const
{
  if (r > paramskk(i,j,k).bigr + paramskk(i,j,k).bigd) return 0.0;
  return -paramskk(i,j,k).bigb * exp(-paramskk(i,j,k).lam2 * r)
          * ters_fc_k(i,j,k,r);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffKokkos<DeviceType>::ters_dfa(const int &i, const int &j,
                const int &k, const F_FLOAT &r) const
{
  if (r > paramskk(i,j,k).bigr + paramskk(i,j,k).bigd) return 0.0;
  return paramskk(i,j,k).bigb * exp(-paramskk(i,j,k).lam2 * r) *
    (paramskk(i,j,k).lam2 * ters_fc_k(i,j,k,r) - ters_dfc(i,j,k,r));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffKokkos<DeviceType>::ters_bij_k(const int &i, const int &j,
                const int &k, const F_FLOAT &bo) const
{
  const F_FLOAT tmp = paramskk(i,j,k).beta * bo;
  if (tmp > paramskk(i,j,k).c1) return 1.0/sqrt(tmp);
  if (tmp > paramskk(i,j,k).c2)
    return (1.0 - pow(tmp,-paramskk(i,j,k).powern) / (2.0*paramskk(i,j,k).powern))/sqrt(tmp);
  if (tmp < paramskk(i,j,k).c4) return 1.0;
  if (tmp < paramskk(i,j,k).c3)
    return 1.0 - pow(tmp,paramskk(i,j,k).powern)/(2.0*paramskk(i,j,k).powern);
  return pow(1.0 + pow(tmp,paramskk(i,j,k).powern), -1.0/(2.0*paramskk(i,j,k).powern));
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairTersoffKokkos<DeviceType>::ters_dbij(const int &i, const int &j,
                const int &k, const F_FLOAT &bo) const
{
  const F_FLOAT tmp = paramskk(i,j,k).beta * bo;
  if (tmp > paramskk(i,j,k).c1) return paramskk(i,j,k).beta * -0.5/sqrt(tmp*tmp);//*pow(tmp,-1.5);
  if (tmp > paramskk(i,j,k).c2)
    return paramskk(i,j,k).beta * (-0.5/sqrt(tmp*tmp) * //*pow(tmp,-1.5) *
           (1.0 - 0.5*(1.0 +  1.0/(2.0*paramskk(i,j,k).powern)) *
           pow(tmp,-paramskk(i,j,k).powern)));
  if (tmp < paramskk(i,j,k).c4) return 0.0;
  if (tmp < paramskk(i,j,k).c3)
    return -0.5*paramskk(i,j,k).beta * pow(tmp,paramskk(i,j,k).powern-1.0);

  const F_FLOAT tmp_n = pow(tmp,paramskk(i,j,k).powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*paramskk(i,j,k).powern)))*tmp_n / bo;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::ters_dthb(
        const int &i, const int &j, const int &k, const F_FLOAT &prefactor,
        const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
        const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2,
        F_FLOAT *fi, F_FLOAT *fj, F_FLOAT *fk) const
{
  // from PairTersoff::attractive
  F_FLOAT rij_hat[3],rik_hat[3];
  F_FLOAT rijinv,rikinv;
  F_FLOAT delrij[3], delrik[3];

  delrij[0] = dx1; delrij[1] = dy1; delrij[2] = dz1;
  delrik[0] = dx2; delrik[1] = dy2; delrik[2] = dz2;

  //rij = sqrt(rsq1);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  //rik = sqrt(rsq2);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  // from PairTersoff::ters_zetaterm_d
  F_FLOAT gijk,dgijk,ex_delr,dex_delr,fc,dfc,cos,tmp;
  F_FLOAT dcosfi[3],dcosfj[3],dcosfk[3];

  fc = ters_fc_k(i,j,k,rik);
  dfc = ters_dfc(i,j,k,rik);
  const F_FLOAT param = paramskk(i,j,k).lam3 * (rij-rik);
  if (int(paramskk(i,j,k).powerm) == 3) tmp = param*param*param;//pow(paramskk(i,j,k).lam3 * (rij-rik),3.0);
  else tmp = param;

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (int(paramskk(i,j,k).powerm) == 3)
    dex_delr = 3.0*param*param*paramskk(i,j,k).lam3*ex_delr;//pow(rij-rik,2.0)*ex_delr;
  else dex_delr = paramskk(i,j,k).lam3 * ex_delr;

  cos = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(i,j,k,cos);
  dgijk = ters_dgijk(i,j,k,cos);

  // from PairTersoff::costheta_d
  vec3_scaleadd(-cos,rij_hat,rik_hat,dcosfj);
  vec3_scale(rijinv,dcosfj,dcosfj);
  vec3_scaleadd(-cos,rik_hat,rij_hat,dcosfk);
  vec3_scale(rikinv,dcosfk,dcosfk);
  vec3_add(dcosfj,dcosfk,dcosfi);
  vec3_scale(-1.0,dcosfi,dcosfi);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,fi);
  vec3_scaleadd(fc*dgijk*ex_delr,dcosfi,fi,fi);
  vec3_scaleadd(fc*gijk*dex_delr,rik_hat,fi,fi);
  vec3_scaleadd(-fc*gijk*dex_delr,rij_hat,fi,fi);
  vec3_scale(prefactor,fi,fi);

  vec3_scale(fc*dgijk*ex_delr,dcosfj,fj);
  vec3_scaleadd(fc*gijk*dex_delr,rij_hat,fj,fj);
  vec3_scale(prefactor,fj,fj);

  vec3_scale(dfc*gijk*ex_delr,rik_hat,fk);
  vec3_scaleadd(fc*dgijk*ex_delr,dcosfk,fk,fk);
  vec3_scaleadd(-fc*gijk*dex_delr,rik_hat,fk,fk);
  vec3_scale(prefactor,fk,fk);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::ters_dthbj(
        const int &i, const int &j, const int &k, const F_FLOAT &prefactor,
        const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
        const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2,
        F_FLOAT *fj, F_FLOAT *fk) const
{
  F_FLOAT rij_hat[3],rik_hat[3];
  F_FLOAT rijinv,rikinv;
  F_FLOAT delrij[3], delrik[3];

  delrij[0] = dx1; delrij[1] = dy1; delrij[2] = dz1;
  delrik[0] = dx2; delrik[1] = dy2; delrik[2] = dz2;

  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  F_FLOAT gijk,dgijk,ex_delr,dex_delr,fc,dfc,cos,tmp;
  F_FLOAT dcosfi[3],dcosfj[3],dcosfk[3];

  fc = ters_fc_k(i,j,k,rik);
  dfc = ters_dfc(i,j,k,rik);
  const F_FLOAT param = paramskk(i,j,k).lam3 * (rij-rik);
  if (int(paramskk(i,j,k).powerm) == 3) tmp = param*param*param;//pow(paramskk(i,j,k).lam3 * (rij-rik),3.0);
  else tmp = param;//paramskk(i,j,k).lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (int(paramskk(i,j,k).powerm) == 3)
    dex_delr = 3.0*param*param*paramskk(i,j,k).lam3*ex_delr;//pow(paramskk(i,j,k).lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else dex_delr = paramskk(i,j,k).lam3 * ex_delr;

  cos = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(i,j,k,cos);
  dgijk = ters_dgijk(i,j,k,cos);

  vec3_scaleadd(-cos,rij_hat,rik_hat,dcosfj);
  vec3_scale(rijinv,dcosfj,dcosfj);
  vec3_scaleadd(-cos,rik_hat,rij_hat,dcosfk);
  vec3_scale(rikinv,dcosfk,dcosfk);
  vec3_add(dcosfj,dcosfk,dcosfi);
  vec3_scale(-1.0,dcosfi,dcosfi);

  vec3_scale(fc*dgijk*ex_delr,dcosfj,fj);
  vec3_scaleadd(fc*gijk*dex_delr,rij_hat,fj,fj);
  vec3_scale(prefactor,fj,fj);

  vec3_scale(dfc*gijk*ex_delr,rik_hat,fk);
  vec3_scaleadd(fc*dgijk*ex_delr,dcosfk,fk,fk);
  vec3_scaleadd(-fc*gijk*dex_delr,rik_hat,fk,fk);
  vec3_scale(prefactor,fk,fk);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::ters_dthbk(
        const int &i, const int &j, const int &k, const F_FLOAT &prefactor,
        const F_FLOAT &rij, const F_FLOAT &dx1, const F_FLOAT &dy1, const F_FLOAT &dz1,
        const F_FLOAT &rik, const F_FLOAT &dx2, const F_FLOAT &dy2, const F_FLOAT &dz2,
        F_FLOAT *fk) const
{
  F_FLOAT rij_hat[3],rik_hat[3];
  F_FLOAT rijinv,rikinv;
  F_FLOAT delrij[3], delrik[3];

  delrij[0] = dx1; delrij[1] = dy1; delrij[2] = dz1;
  delrik[0] = dx2; delrik[1] = dy2; delrik[2] = dz2;

  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  F_FLOAT gijk,dgijk,ex_delr,dex_delr,fc,dfc,cos,tmp;
  F_FLOAT dcosfi[3],dcosfj[3],dcosfk[3];

  fc = ters_fc_k(i,j,k,rik);
  dfc = ters_dfc(i,j,k,rik);
  const F_FLOAT param = paramskk(i,j,k).lam3 * (rij-rik);
  if (int(paramskk(i,j,k).powerm) == 3) tmp = param*param*param;//pow(paramskk(i,j,k).lam3 * (rij-rik),3.0);
  else tmp = param;//paramskk(i,j,k).lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (int(paramskk(i,j,k).powerm) == 3)
    dex_delr = 3.0*param*param*paramskk(i,j,k).lam3*ex_delr;//pow(paramskk(i,j,k).lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else dex_delr = paramskk(i,j,k).lam3 * ex_delr;

  cos = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(i,j,k,cos);
  dgijk = ters_dgijk(i,j,k,cos);

  vec3_scaleadd(-cos,rij_hat,rik_hat,dcosfj);
  vec3_scale(rijinv,dcosfj,dcosfj);
  vec3_scaleadd(-cos,rik_hat,rij_hat,dcosfk);
  vec3_scale(rikinv,dcosfk,dcosfk);
  vec3_add(dcosfj,dcosfk,dcosfi);
  vec3_scale(-1.0,dcosfi,dcosfi);

  vec3_scale(dfc*gijk*ex_delr,rik_hat,fk);
  vec3_scaleadd(fc*dgijk*ex_delr,dcosfk,fk,fk);
  vec3_scaleadd(-fc*gijk*dex_delr,rik_hat,fk,fk);
  vec3_scale(prefactor,fk,fk);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
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
    if (NEIGHFLAG != FULL) a_eatom[j] += epairhalf;
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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::v_tally3(EV_FLOAT &ev, const int &i, const int &j, const int &k,
        F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drij, F_FLOAT *drik) const
{
  // The vatom array is duplicated for OpenMP, atomic for CUDA, and neither for Serial

  auto v_vatom = ScatterViewHelper<NeedDup<NEIGHFLAG,DeviceType>::value,decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup<NEIGHFLAG,DeviceType>::value>();

  F_FLOAT v[6];

  v[0] = THIRD * (drij[0]*fj[0] + drik[0]*fk[0]);
  v[1] = THIRD * (drij[1]*fj[1] + drik[1]*fk[1]);
  v[2] = THIRD * (drij[2]*fj[2] + drik[2]*fk[2]);
  v[3] = THIRD * (drij[0]*fj[1] + drik[0]*fk[1]);
  v[4] = THIRD * (drij[0]*fj[2] + drik[0]*fk[2]);
  v[5] = THIRD * (drij[1]*fj[2] + drik[1]*fk[2]);

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    a_vatom(i,0) += v[0]; a_vatom(i,1) += v[1]; a_vatom(i,2) += v[2];
    a_vatom(i,3) += v[3]; a_vatom(i,4) += v[4]; a_vatom(i,5) += v[5];
    if (NEIGHFLAG != FULL) {
      a_vatom(j,0) += v[0]; a_vatom(j,1) += v[1]; a_vatom(j,2) += v[2];
      a_vatom(j,3) += v[3]; a_vatom(j,4) += v[4]; a_vatom(j,5) += v[5];
      a_vatom(k,0) += v[0]; a_vatom(k,1) += v[1]; a_vatom(k,2) += v[2];
      a_vatom(k,3) += v[3]; a_vatom(k,4) += v[4]; a_vatom(k,5) += v[5];
    }
  }

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairTersoffKokkos<DeviceType>::v_tally3_atom(EV_FLOAT &ev, const int &i, const int &j, const int &k,
        F_FLOAT *fj, F_FLOAT *fk, F_FLOAT *drji, F_FLOAT *drjk) const
{
  F_FLOAT v[6];

  v[0] = THIRD * (drji[0]*fj[0] + drjk[0]*fk[0]);
  v[1] = THIRD * (drji[1]*fj[1] + drjk[1]*fk[1]);
  v[2] = THIRD * (drji[2]*fj[2] + drjk[2]*fk[2]);
  v[3] = THIRD * (drji[0]*fj[1] + drjk[0]*fk[1]);
  v[4] = THIRD * (drji[0]*fj[2] + drjk[0]*fk[2]);
  v[5] = THIRD * (drji[1]*fj[2] + drjk[1]*fk[2]);

  if (vflag_global) {
    ev.v[0] += v[0];
    ev.v[1] += v[1];
    ev.v[2] += v[2];
    ev.v[3] += v[3];
    ev.v[4] += v[4];
    ev.v[5] += v[5];
  }

  if (vflag_atom) {
    d_vatom(i,0) += v[0]; d_vatom(i,1) += v[1]; d_vatom(i,2) += v[2];
    d_vatom(i,3) += v[3]; d_vatom(i,4) += v[4]; d_vatom(i,5) += v[5];
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairTersoffKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairTersoffKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairTersoffKokkos<LMPHostType>;
#endif
}

