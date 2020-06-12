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

/* ----------------------------------------------------------------------------------------
   Contributing authors:
   Stan Moore (Sandia)

   Please cite the related publications:
   J.D. Moore, B.C. Barnes, S. Izvekov, M. Lisal, M.S. Sellers, D.E. Taylor & J.K. Brennan
   "A coarse-grain force field for RDX: Density dependent and energy conserving"
   The Journal of Chemical Physics, 2016, 144, 104501.
------------------------------------------------------------------------------------------- */

#include "pair_multi_lucy_rx_kokkos.h"
#include <cmath>
#include <cstring>
#include "math_const.h"
#include "atom_kokkos.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"
#include "neigh_request.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

enum{NONE,RLINEAR,RSQ};

#define MAXLINE 1024

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define oneFluidParameter (-1)
#define isOneFluid(_site) ( (_site) == oneFluidParameter )

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairMultiLucyRXKokkos<Space>::PairMultiLucyRXKokkos(LAMMPS *lmp) : PairMultiLucyRX(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  update_table = 1;
  k_table = new TableDual();
  h_table = new TableHost(k_table);
  d_table = new TableDevice(k_table);

  k_error_flag = DAT::tdual_int_scalar("pair:error_flag");
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
PairMultiLucyRXKokkos<Space>::~PairMultiLucyRXKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);

  memoryKK->destroy_kokkos(k_cutsq,cutsq);

  delete k_table;
  delete h_table;
  delete d_table;
  tabindex = NULL;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairMultiLucyRXKokkos<Space>::init_style()
{
  PairMultiLucyRX::init_style();

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = (Space == Host) &&
    !(Space == Device);
  neighbor->requests[irequest]->
    kokkos_device = (Space == Device);

  if (neighflag == FULL) {
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
  } else if (neighflag == HALF || neighflag == HALFTHREAD) {
    neighbor->requests[irequest]->full = 0;
    neighbor->requests[irequest]->half = 1;
  } else {
    error->all(FLERR,"Cannot use chosen neighbor list style with multi/lucy/rx/kk");
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairMultiLucyRXKokkos<Space>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  if (update_table)
    create_kokkos_tables();

  if (tabstyle == LOOKUP)
    compute_style<LOOKUP>(eflag_in,vflag_in);
  else if(tabstyle == LINEAR)
    compute_style<LINEAR>(eflag_in,vflag_in);

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
template<int TABSTYLE>
void PairMultiLucyRXKokkos<Space>::compute_style(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  ev_init(eflag,vflag,0);

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

  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  rho = DualViewHelper<Space>::view(atomKK->k_rho);
  uCG = DualViewHelper<Space>::view(atomKK->k_uCG);
  uCGnew = DualViewHelper<Space>::view(atomKK->k_uCGnew);
  dvector = DualViewHelper<Space>::view(atomKK->k_dvector);

  atomKK->sync(Space,X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK | DPDRHO_MASK | UCG_MASK | UCGNEW_MASK | DVECTOR_MASK);
  DualViewHelper<Space>::sync(k_cutsq);

  nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int newton_pair = force->newton_pair;

  {
    const int ntotal = nlocal + nghost;
    if (ntotal > d_mixWtSite1.extent(0)) {
      d_mixWtSite1old = typename AT::t_float_1d("PairMultiLucyRX::mixWtSite1old",ntotal);
      d_mixWtSite2old = typename AT::t_float_1d("PairMultiLucyRX::mixWtSite2old",ntotal);
      d_mixWtSite1 = typename AT::t_float_1d("PairMultiLucyRX::mixWtSite1",ntotal);
      d_mixWtSite2 = typename AT::t_float_1d("PairMultiLucyRX::mixWtSite2",ntotal);
    }

    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXgetMixingWeights>(0,ntotal),*this);
  }

  const int inum = list->inum;
  NeighListKokkos<Space>* k_list = static_cast<NeighListKokkos<Space>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  computeLocalDensity();

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (neighflag == HALF) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALF,1,1,TABSTYLE> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALF,1,0,TABSTYLE> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALF,0,1,TABSTYLE> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALF,0,0,TABSTYLE> >(0,inum),*this);
    }
  } else if (neighflag == HALFTHREAD) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALFTHREAD,1,1,TABSTYLE> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALFTHREAD,1,0,TABSTYLE> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALFTHREAD,0,1,TABSTYLE> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<HALFTHREAD,0,0,TABSTYLE> >(0,inum),*this);
    }
  } else if (neighflag == FULL) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<FULL,1,1,TABSTYLE> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<FULL,1,0,TABSTYLE> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<FULL,0,1,TABSTYLE> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXCompute<FULL,0,0,TABSTYLE> >(0,inum),*this);
    }
  }

  if (evflag) atomKK->modified(Space,F_MASK | ENERGY_MASK | VIRIAL_MASK | UCG_MASK | UCGNEW_MASK);
  else atomKK->modified(Space,F_MASK | UCG_MASK | UCGNEW_MASK);

  DualViewHelper<Space>::modify(k_error_flag);
  k_error_flag.sync_host();
  if (k_error_flag.h_view() == 1)
    error->one(FLERR,"Density < table inner cutoff");
  else if (k_error_flag.h_view() == 2)
    error->one(FLERR,"Density > table outer cutoff");
  else if (k_error_flag.h_view() == 3)
    error->one(FLERR,"Only LOOKUP and LINEAR table styles have been implemented for pair multi/lucy/rx");

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute<Space>(this);

  if (eflag_atom) {
    DualViewHelper<Space>::modify(k_eatom);
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    DualViewHelper<Space>::modify(k_vatom);
    k_vatom.sync_host();
  }
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<Space>::operator()(TagPairMultiLucyRXgetMixingWeights, const int &i) const {
  getMixingWeights(i, d_mixWtSite1old[i], d_mixWtSite2old[i], d_mixWtSite1[i], d_mixWtSite2[i]);
}

template<ExecutionSpace Space>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<Space>::operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>, const int &ii, EV_FLOAT& ev) const {

  // The f array is atomic for Half/Thread neighbor style
  Kokkos::View<typename AT::t_float_1d_3::data_type, typename AT::t_float_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;

  int i,jj,jnum,itype,jtype,itable;
  KK_FLOAT xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwlOld,fpair;
  KK_FLOAT rsq;

  KK_FLOAT mixWtSite1old_i,mixWtSite1old_j;
  KK_FLOAT mixWtSite2old_i,mixWtSite2old_j;
  KK_FLOAT mixWtSite1_i;

  KK_FLOAT pi = MathConst::MY_PI;
  KK_FLOAT A_i, A_j;
  KK_FLOAT fraction_i,fraction_j;
  int jtable;

  int tlm1 = tablength - 1;

  i = d_ilist[ii];
  xtmp = x(i,0);
  ytmp = x(i,1);
  ztmp = x(i,2);
  itype = type[i];
  jnum = d_numneigh[i];

  KK_FLOAT fx_i = 0.0;
  KK_FLOAT fy_i = 0.0;
  KK_FLOAT fz_i = 0.0;

  mixWtSite1old_i = d_mixWtSite1old[i];
  mixWtSite2old_i = d_mixWtSite2old[i];
  mixWtSite1_i = d_mixWtSite1[i];

  for (jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;

    delx = xtmp - x(j,0);
    dely = ytmp - x(j,1);
    delz = ztmp - x(j,2);
    rsq = delx*delx + dely*dely + delz*delz;
    jtype = type[j];

    if (rsq < d_cutsq(itype,jtype)) { // optimize
      fpair = 0.0;

      mixWtSite1old_j = d_mixWtSite1old[j];
      mixWtSite2old_j = d_mixWtSite2old[j];

      //tb = &tables[tabindex[itype][jtype]];
      const int tidx = d_table_const.tabindex(itype,jtype);

      //if (rho[i]*rho[i] < tb->innersq || rho[j]*rho[j] < tb->innersq){
      if (rho[i]*rho[i] < d_table_const.innersq(tidx) || rho[j]*rho[j] < d_table_const.innersq(tidx)){
        DualViewHelper<Space>::view(k_error_flag)() = 1;
      }

      if (TABSTYLE == LOOKUP) {
        //itable = static_cast<int> (((rho[i]*rho[i]) - tb->innersq) * tb->invdelta);
        itable = static_cast<int> (((rho[i]*rho[i]) - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
        //jtable = static_cast<int> (((rho[j]*rho[j]) - tb->innersq) * tb->invdelta);
        jtable = static_cast<int> (((rho[j]*rho[j]) - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
        if (itable >= tlm1 || jtable >= tlm1){
          DualViewHelper<Space>::view(k_error_flag)() = 2;
        }
        //A_i = tb->f[itable];
        A_i = d_table_const.f(tidx,itable);
        //A_j = tb->f[jtable];
        A_j = d_table_const.f(tidx,jtable);

        const KK_FLOAT rfactor = 1.0-sqrt(rsq/d_cutsq(itype,jtype));
        fpair = 0.5*(A_i + A_j)*(4.0-3.0*rfactor)*rfactor*rfactor*rfactor;
        fpair /= sqrt(rsq);

      } else if (TABSTYLE == LINEAR) {

        //itable = static_cast<int> ((rho[i]*rho[i] - tb->innersq) * tb->invdelta);
        itable = static_cast<int> ((rho[i]*rho[i] - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
        //jtable = static_cast<int> (((rho[j]*rho[j]) - tb->innersq) * tb->invdelta);
        jtable = static_cast<int> ((rho[j]*rho[j] - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
        if (itable >= tlm1 || jtable >= tlm1){
          DualViewHelper<Space>::view(k_error_flag)() = 2;
        }
        if(itable<0) itable=0;
        if(itable>=tlm1) itable=tlm1;
        if(jtable<0) jtable=0;
        if(jtable>=tlm1)jtable=tlm1;

        //fraction_i = (((rho[i]*rho[i]) - tb->rsq[itable]) * tb->invdelta);
        fraction_i = (((rho[i]*rho[i]) - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx));
        //fraction_j = (((rho[j]*rho[j]) - tb->rsq[jtable]) * tb->invdelta);
        fraction_j = (((rho[j]*rho[j]) - d_table_const.rsq(tidx,jtable)) * d_table_const.invdelta(tidx));
        if(itable==0) fraction_i=0.0;
        if(itable==tlm1) fraction_i=0.0;
        if(jtable==0) fraction_j=0.0;
        if(jtable==tlm1) fraction_j=0.0;

        //A_i = tb->f[itable] + fraction_i*tb->df[itable];
        A_i = d_table_const.f(tidx,itable) + fraction_i*d_table_const.df(tidx,itable);
        //A_j = tb->f[jtable] + fraction_j*tb->df[jtable];
        A_j = d_table_const.f(tidx,jtable) + fraction_j*d_table_const.df(tidx,jtable);

        const KK_FLOAT rfactor = 1.0-sqrt(rsq/d_cutsq(itype,jtype));
        fpair = 0.5*(A_i + A_j)*(4.0-3.0*rfactor)*rfactor*rfactor*rfactor;
        fpair /= sqrt(rsq);

      } else DualViewHelper<Space>::view(k_error_flag)() = 3;

      if (isite1 == isite2) fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpair;
      else fpair = (sqrt(mixWtSite1old_i*mixWtSite2old_j) + sqrt(mixWtSite2old_i*mixWtSite1old_j))*fpair;

      fx_i += delx*fpair;
      fy_i += dely*fpair;
      fz_i += delz*fpair;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        a_f(j,0) -= delx*fpair;
        a_f(j,1) -= dely*fpair;
        a_f(j,2) -= delz*fpair;
      }
      //if (evflag) ev_tally(i,j,nlocal,newton_pair,0.0,0.0,fpair,delx,dely,delz);
      if (EVFLAG) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,0.0,fpair,delx,dely,delz);
    }
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;

  //tb = &tables[tabindex[itype][itype]];
  const int tidx = d_table_const.tabindex(itype,itype);
  //itable = static_cast<int> (((rho[i]*rho[i]) - tb->innersq) * tb->invdelta);
  itable = static_cast<int> (((rho[i]*rho[i]) - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
  //if (TABSTYLE == LOOKUP) evdwl = tb->e[itable];
  if (TABSTYLE == LOOKUP) {
    evdwl = d_table_const.e(tidx,itable);
  } else if (TABSTYLE == LINEAR) {
    if (itable >= tlm1){
      DualViewHelper<Space>::view(k_error_flag)() = 2;
    }
    if(itable==0) fraction_i=0.0;
    //else fraction_i = (((rho[i]*rho[i]) - tb->rsq[itable]) * tb->invdelta);
    else fraction_i = (((rho[i]*rho[i]) - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx));
    //evdwl = tb->e[itable] + fraction_i*tb->de[itable];
    evdwl = d_table_const.e(tidx,itable) + fraction_i*d_table_const.de(tidx,itable);
  } else DualViewHelper<Space>::view(k_error_flag)() = 3;

  evdwl *=(pi*d_cutsq(itype,itype)*d_cutsq(itype,itype))/84.0;
  evdwlOld = mixWtSite1old_i*evdwl;
  evdwl = mixWtSite1_i*evdwl;

  uCG[i] += evdwlOld;
  uCGnew[i] += evdwl;

  evdwl = evdwlOld;

  //if (evflag) ev_tally(0,0,nlocal,newton_pair,evdwl,0.0,0.0,0.0,0.0,0.0);
  if (EVFLAG)
    ev.evdwl += ((/*FIXME??? (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && */ NEWTON_PAIR)?1.0:0.5)*evdwl;
}

template<ExecutionSpace Space>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<Space>::operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairMultiLucyRXKokkos<Space>::computeLocalDensity()
{
  x = DualViewHelper<Space>::view(atomKK->k_x);
  type = DualViewHelper<Space>::view(atomKK->k_type);
  rho = DualViewHelper<Space>::view(atomKK->k_rho);
  h_rho = atomKK->k_rho.h_view;
  nlocal = atom->nlocal;

  atomKK->sync(Space,X_MASK | TYPE_MASK | DPDRHO_MASK);

  const int inum = list->inum;
  NeighListKokkos<Space>* k_list = static_cast<NeighListKokkos<Space>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  const KK_FLOAT pi = MathConst::MY_PI;

  const bool newton_pair = force->newton_pair;
  const bool one_type = (atom->ntypes == 1);

  // Special cut-off values for when there's only one type.
  cutsq_type11 = cutsq[1][1];
  rcut_type11 = sqrt(cutsq_type11);
  factor_type11 = 84.0/(5.0*pi*rcut_type11*rcut_type11*rcut_type11);

  // zero out density
  int m = nlocal;
  if (newton_pair) m += atom->nghost;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXZero>(0,m),*this);

  // rho = density at each atom
  // loop over neighbors of my atoms

  if (neighflag == HALF) {
    if (newton_pair)
      if (one_type)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALF,1,true> >(0,inum),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALF,1,false> >(0,inum),*this);
    else
      if (one_type)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALF,0,true> >(0,inum),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALF,0,false> >(0,inum),*this);
  } else if (neighflag == HALFTHREAD) {
    if (newton_pair)
      if (one_type)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALFTHREAD,1,true> >(0,inum),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALFTHREAD,1,false> >(0,inum),*this);
    else
      if (one_type)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALFTHREAD,0,true> >(0,inum),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<HALFTHREAD,0,false> >(0,inum),*this);
  } else if (neighflag == FULL) {
    if (newton_pair)
      if (one_type)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<FULL,1,true> >(0,inum),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<FULL,1,false> >(0,inum),*this);
    else
      if (one_type)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<FULL,0,true> >(0,inum),*this);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXComputeLocalDensity<FULL,0,false> >(0,inum),*this);
  }

  atomKK->modified(Space,DPDRHO_MASK);

  // communicate and sum densities (on the host)

  if (newton_pair)
    comm->reverse_comm_pair(this);

  comm->forward_comm_pair(this);
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<Space>::operator()(TagPairMultiLucyRXZero, const int &i) const {
  rho[i] = 0.0;
}

template<ExecutionSpace Space>
template<int NEIGHFLAG, int NEWTON_PAIR, bool ONE_TYPE>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<Space>::operator()(TagPairMultiLucyRXComputeLocalDensity<NEIGHFLAG,NEWTON_PAIR,ONE_TYPE>, const int &ii) const {


  // The rho array is atomic for Half/Thread neighbor style
  Kokkos::View<typename AT::t_float_1d::data_type, typename AT::t_float_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_rho = rho;

  const int i = d_ilist[ii];

  const KK_FLOAT xtmp = x(i,0);
  const KK_FLOAT ytmp = x(i,1);
  const KK_FLOAT ztmp = x(i,2);

  KK_FLOAT rho_i_contrib = 0.0;

  const int itype = type[i];
  const int jnum = d_numneigh[i];

  const KK_FLOAT pi = MathConst::MY_PI;

  for (int jj = 0; jj < jnum; jj++){
    const int j = (d_neighbors(i,jj) & NEIGHMASK);
    const int jtype = type[j];

    const KK_FLOAT delx = xtmp - x(j,0);
    const KK_FLOAT dely = ytmp - x(j,1);
    const KK_FLOAT delz = ztmp - x(j,2);
    const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;

    if (ONE_TYPE) {
      if (rsq < cutsq_type11) {
        const KK_FLOAT rcut = rcut_type11;
        const KK_FLOAT r_over_rcut = sqrt(rsq) / rcut;
        const KK_FLOAT tmpFactor = 1.0 - r_over_rcut;
        const KK_FLOAT tmpFactor4 = tmpFactor*tmpFactor*tmpFactor*tmpFactor;
        const KK_FLOAT factor = factor_type11*(1.0 + 1.5*r_over_rcut)*tmpFactor4;
        rho_i_contrib += factor;
        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
          a_rho[j] += factor;
      }
    } else if (rsq < d_cutsq(itype,jtype)) {
      const KK_FLOAT rcut = sqrt(d_cutsq(itype,jtype));
      const KK_FLOAT tmpFactor = 1.0-sqrt(rsq)/rcut;
      const KK_FLOAT tmpFactor4 = tmpFactor*tmpFactor*tmpFactor*tmpFactor;
      const KK_FLOAT factor = (84.0/(5.0*pi*rcut*rcut*rcut))*(1.0+3.0*sqrt(rsq)/(2.0*rcut))*tmpFactor4;
      rho_i_contrib += factor;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
        a_rho[j] += factor;
    }
  }

  a_rho[i] += rho_i_contrib;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<Space>::getMixingWeights(int id, SPACE_FLOAT &mixWtSite1old, SPACE_FLOAT &mixWtSite2old, SPACE_FLOAT &mixWtSite1, SPACE_FLOAT &mixWtSite2) const
{
  KK_FLOAT fractionOFAold, fractionOFA;
  KK_FLOAT fractionOld1, fraction1;
  KK_FLOAT fractionOld2, fraction2;
  KK_FLOAT nMoleculesOFAold, nMoleculesOFA;
  KK_FLOAT nMoleculesOld1, nMolecules1;
  KK_FLOAT nMoleculesOld2, nMolecules2;
  KK_FLOAT nTotal, nTotalOld;


  nTotal = 0.0;
  nTotalOld = 0.0;
  for (int ispecies = 0; ispecies < nspecies; ispecies++){
    nTotal += dvector(ispecies,id);
    nTotalOld += dvector(ispecies+nspecies,id);
  }

  if (isOneFluid(isite1) == false){
    nMoleculesOld1 = dvector(isite1+nspecies,id);
    nMolecules1 = dvector(isite1,id);
    fractionOld1 = nMoleculesOld1/nTotalOld;
    fraction1 = nMolecules1/nTotal;
  }
  if (isOneFluid(isite2) == false){
    nMoleculesOld2 = dvector(isite2+nspecies,id);
    nMolecules2 = dvector(isite2,id);
    fractionOld2 = nMoleculesOld2/nTotalOld;
    fraction2 = nMolecules2/nTotal;
  }

  if (isOneFluid(isite1) || isOneFluid(isite2)){
    nMoleculesOFAold  = 0.0;
    nMoleculesOFA  = 0.0;
    fractionOFAold  = 0.0;
    fractionOFA  = 0.0;

    for (int ispecies = 0; ispecies < nspecies; ispecies++){
      if (isite1 == ispecies || isite2 == ispecies) continue;
      nMoleculesOFAold += dvector(ispecies+nspecies,id);
      nMoleculesOFA += dvector(ispecies,id);
      fractionOFAold += dvector(ispecies+nspecies,id) / nTotalOld;
      fractionOFA += dvector(ispecies,id) / nTotal;
    }
    if (isOneFluid(isite1)){
      nMoleculesOld1 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules1 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld1 = fractionOFAold;
      fraction1 = fractionOFA;
    }
    if (isOneFluid(isite2)){
      nMoleculesOld2 = 1.0-(nTotalOld-nMoleculesOFAold);
      nMolecules2 = 1.0-(nTotal-nMoleculesOFA);
      fractionOld2 = fractionOFAold;
      fraction2 = fractionOFA;
    }
  }

  if(fractionalWeighting){
    mixWtSite1old = fractionOld1;
    mixWtSite1 = fraction1;
    mixWtSite2old = fractionOld2;
    mixWtSite2 = fraction2;
  } else {
    mixWtSite1old = nMoleculesOld1;
    mixWtSite1 = nMolecules1;
    mixWtSite2old = nMoleculesOld2;
    mixWtSite2 = nMolecules2;
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
int PairMultiLucyRXKokkos<Space>::pack_forward_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist, int iswap_in, DAT::tdual_float_1d &buf,
                               int pbc_flag, int *pbc)
{
  atomKK->sync(Space,DPDRHO_MASK);

  d_sendlist = DualViewHelper<Space>::view(k_sendlist);
  iswap = iswap_in;
  v_buf = DualViewHelper<Space>::view(buf);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXPackForwardComm>(0,n),*this);
  return n;
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<Space>::operator()(TagPairMultiLucyRXPackForwardComm, const int &i) const {
  int j = d_sendlist(iswap, i);
  v_buf[i] = rho[j];
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairMultiLucyRXKokkos<Space>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_float_1d &buf)
{
  first = first_in;
  v_buf = DualViewHelper<Space>::view(buf);
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairMultiLucyRXUnpackForwardComm>(0,n),*this);

  atomKK->modified(Space,DPDRHO_MASK);
}

template<ExecutionSpace Space>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<Space>::operator()(TagPairMultiLucyRXUnpackForwardComm, const int &i) const {
  rho[i + first] = v_buf[i];
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
int PairMultiLucyRXKokkos<Space>::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;

  atomKK->sync(Host,DPDRHO_MASK);

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = h_rho[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairMultiLucyRXKokkos<Space>::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) h_rho[i] = buf[m++];

  atomKK->modified(Host,DPDRHO_MASK);
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
int PairMultiLucyRXKokkos<Space>::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  atomKK->sync(Host,DPDRHO_MASK);

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = h_rho[i];
  return m;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairMultiLucyRXKokkos<Space>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    h_rho[j] += buf[m++];
  }

  atomKK->modified(Host,DPDRHO_MASK);
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<Space>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are atomic for Half/Thread neighbor style
  Kokkos::View<typename AT::t_float_1d::data_type, typename AT::t_float_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = DualViewHelper<Space>::view(k_eatom);
  Kokkos::View<typename AT::t_float_1d_6::data_type, typename AT::t_float_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = DualViewHelper<Space>::view(k_vatom);

  if (EFLAG) {
    if (eflag_atom) {
      const KK_FLOAT epairhalf = 0.5 * epair;
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) v_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) v_eatom[j] += epairhalf;
      } else {
        v_eatom[i] += epairhalf;
      }
    }
  }

  if (VFLAG) {
    const KK_FLOAT v0 = delx*delx*fpair;
    const KK_FLOAT v1 = dely*dely*fpair;
    const KK_FLOAT v2 = delz*delz*fpair;
    const KK_FLOAT v3 = delx*dely*fpair;
    const KK_FLOAT v4 = delx*delz*fpair;
    const KK_FLOAT v5 = dely*delz*fpair;

    if (vflag_global) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          ev.v[0] += 0.5*v0;
          ev.v[1] += 0.5*v1;
          ev.v[2] += 0.5*v2;
          ev.v[3] += 0.5*v3;
          ev.v[4] += 0.5*v4;
          ev.v[5] += 0.5*v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
        ev.v[0] += 0.5*v0;
        ev.v[1] += 0.5*v1;
        ev.v[2] += 0.5*v2;
        ev.v[3] += 0.5*v3;
        ev.v[4] += 0.5*v4;
        ev.v[5] += 0.5*v5;
        }
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
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          v_vatom(i,0) += 0.5*v0;
          v_vatom(i,1) += 0.5*v1;
          v_vatom(i,2) += 0.5*v2;
          v_vatom(i,3) += 0.5*v3;
          v_vatom(i,4) += 0.5*v4;
          v_vatom(i,5) += 0.5*v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
        v_vatom(j,0) += 0.5*v0;
        v_vatom(j,1) += 0.5*v1;
        v_vatom(j,2) += 0.5*v2;
        v_vatom(j,3) += 0.5*v3;
        v_vatom(j,4) += 0.5*v4;
        v_vatom(j,5) += 0.5*v5;
        }
      } else {
        v_vatom(i,0) += 0.5*v0;
        v_vatom(i,1) += 0.5*v1;
        v_vatom(i,2) += 0.5*v2;
        v_vatom(i,3) += 0.5*v3;
        v_vatom(i,4) += 0.5*v4;
        v_vatom(i,5) += 0.5*v5;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairMultiLucyRXKokkos<Space>::create_kokkos_tables()
{
  const int tlm1 = tablength-1;

  memoryKK->create_kokkos(k_table->k_innersq,ntables,"Table::innersq");
  memoryKK->create_kokkos(k_table->k_invdelta,ntables,"Table::invdelta");

  if(tabstyle == LOOKUP) {
    memoryKK->create_kokkos(k_table->k_e,ntables,tlm1,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,tlm1,"Table::f");
  }

  if(tabstyle == LINEAR) {
    memoryKK->create_kokkos(k_table->k_rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(k_table->k_e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(k_table->k_f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(k_table->k_de,ntables,tlm1,"Table::de");
    memoryKK->create_kokkos(k_table->k_df,ntables,tlm1,"Table::df");
  }

  for(int i=0; i < ntables; i++) {
    Table* tb = &tables[i];

    h_table->innersq[i] = tb->innersq;
    h_table->invdelta[i] = tb->invdelta;

    for(int j = 0; j<h_table->rsq.extent(1); j++)
      h_table->rsq(i,j) = tb->rsq[j];
    for(int j = 0; j<h_table->e.extent(1); j++)
      h_table->e(i,j) = tb->e[j];
    for(int j = 0; j<h_table->de.extent(1); j++)
      h_table->de(i,j) = tb->de[j];
    for(int j = 0; j<h_table->f.extent(1); j++)
      h_table->f(i,j) = tb->f[j];
    for(int j = 0; j<h_table->df.extent(1); j++)
      h_table->df(i,j) = tb->df[j];
  }


  k_table->k_innersq.modify_host();
  k_table->k_innersq.sync_device();
  k_table->k_invdelta.modify_host();
  k_table->k_invdelta.sync_device();
  k_table->k_rsq.modify_host();
  k_table->k_rsq.sync_device();
  k_table->k_e.modify_host();
  k_table->k_e.sync_device();
  k_table->k_de.modify_host();
  k_table->k_de.sync_device();
  k_table->k_f.modify_host();
  k_table->k_f.sync_device();
  k_table->k_df.modify_host();
  k_table->k_df.sync_device();
  k_table->k_tabindex.modify_host();
  k_table->k_tabindex.sync_device();

  d_table_const.innersq = d_table->innersq;
  d_table_const.invdelta = d_table->invdelta;
  d_table_const.rsq = d_table->rsq;
  d_table_const.e = d_table->e;
  d_table_const.de = d_table->de;
  d_table_const.f = d_table->f;
  d_table_const.df = d_table->df;

  update_table = 0;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairMultiLucyRXKokkos<Space>::allocate()
{
  allocated = 1;
  const int nt = atom->ntypes + 1;

  memory->create(setflag,nt,nt,"pair:setflag");

  memoryKK->create_kokkos(k_cutsq,cutsq,nt,nt,"pair:cutsq");
  d_cutsq = DualViewHelper<Space>::view(k_cutsq);
  k_cutsq.modify_host();

  memoryKK->create_kokkos(k_table->k_tabindex,tabindex,nt,nt,"pair:tabindex");
  d_table_const.tabindex = d_table->tabindex;

  memset(&setflag[0][0],0,nt*nt*sizeof(int));
  memset(&cutsq[0][0],0,nt*nt*sizeof(KK_FLOAT));
  memset(&tabindex[0][0],0,nt*nt*sizeof(int));
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void PairMultiLucyRXKokkos<Space>::settings(int narg, char **arg)
{
  if (narg < 2) error->all(FLERR,"Illegal pair_style command");

  // new settings

  if (strcmp(arg[0],"lookup") == 0) tabstyle = LOOKUP;
  else if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
  else error->all(FLERR,"Unknown table style in pair_style command");

  tablength = force->inumeric(FLERR,arg[1]);
  if (tablength < 2) error->all(FLERR,"Illegal number of pair table entries");

  // optional keywords

  int iarg = 2;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"fractional") == 0)   fractionalWeighting = true;
    else if (strcmp(arg[iarg],"molecular") == 0)   fractionalWeighting = false;
    else error->all(FLERR,"Illegal pair_style command");
    iarg++;
  }

  // delete old tables, since cannot just change settings

  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);

    d_table_const.tabindex = d_table->tabindex = typename AT::t_int_2d();
    h_table->tabindex = HAT::t_int_2d();
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairMultiLucyRXKokkos<Device>;
template class PairMultiLucyRXKokkos<Host>;
}
