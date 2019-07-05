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
#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "math_const.h"
#include "atom_kokkos.h"
#include "force.h"
#include "comm.h"
#include "neigh_list.h"
#include "memory_kokkos.h"
#include "error.h"
#include "citeme.h"
#include "modify.h"
#include "fix.h"
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

template<class DeviceType>
PairMultiLucyRXKokkos<DeviceType>::PairMultiLucyRXKokkos(LAMMPS *lmp) : PairMultiLucyRX(lmp)
{
  respa_enable = 0;

  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  update_table = 1;
  h_table = new TableHost();
  d_table = new TableDevice();

  k_error_flag = DAT::tdual_int_scalar("pair:error_flag");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairMultiLucyRXKokkos<DeviceType>::~PairMultiLucyRXKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);

  memoryKK->destroy_kokkos(k_cutsq,cutsq);

  delete h_table;
  delete d_table;
  tabindex = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::init_style()
{
  PairMultiLucyRX::init_style();

  // irequest = neigh request made by parent class

  neighflag = lmp->kokkos->neighflag;
  int irequest = neighbor->nrequest - 1;

  neighbor->requests[irequest]->
    kokkos_host = Kokkos::Impl::is_same<DeviceType,LMPHostType>::value &&
    !Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;
  neighbor->requests[irequest]->
    kokkos_device = Kokkos::Impl::is_same<DeviceType,LMPDeviceType>::value;

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

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

template<class DeviceType>
template<int TABSTYLE>
void PairMultiLucyRXKokkos<DeviceType>::compute_style(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  rho = atomKK->k_rho.view<DeviceType>();
  uCG = atomKK->k_uCG.view<DeviceType>();
  uCGnew = atomKK->k_uCGnew.view<DeviceType>();
  dvector = atomKK->k_dvector.view<DeviceType>();

  atomKK->sync(execution_space,X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK | DPDRHO_MASK | UCG_MASK | UCGNEW_MASK | DVECTOR_MASK);
  k_cutsq.template sync<DeviceType>();

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
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
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

  if (evflag) atomKK->modified(execution_space,F_MASK | ENERGY_MASK | VIRIAL_MASK | UCG_MASK | UCGNEW_MASK);
  else atomKK->modified(execution_space,F_MASK | UCG_MASK | UCGNEW_MASK);

  k_error_flag.template modify<DeviceType>();
  k_error_flag.template sync<LMPHostType>();
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

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXgetMixingWeights, const int &i) const {
  getMixingWeights(i, d_mixWtSite1old[i], d_mixWtSite2old[i], d_mixWtSite1[i], d_mixWtSite2[i]);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>, const int &ii, EV_FLOAT& ev) const {

  // The f array is atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;

  int i,jj,jnum,itype,jtype,itable;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwlOld,fpair;
  double rsq;

  double mixWtSite1old_i,mixWtSite1old_j;
  double mixWtSite2old_i,mixWtSite2old_j;
  double mixWtSite1_i;

  double pi = MathConst::MY_PI;
  double A_i, A_j;
  double fraction_i,fraction_j;
  int jtable;

  int tlm1 = tablength - 1;

  i = d_ilist[ii];
  xtmp = x(i,0);
  ytmp = x(i,1);
  ztmp = x(i,2);
  itype = type[i];
  jnum = d_numneigh[i];

  double fx_i = 0.0;
  double fy_i = 0.0;
  double fz_i = 0.0;

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
        k_error_flag.template view<DeviceType>()() = 1;
      }

      if (TABSTYLE == LOOKUP) {
        //itable = static_cast<int> (((rho[i]*rho[i]) - tb->innersq) * tb->invdelta);
        itable = static_cast<int> (((rho[i]*rho[i]) - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
        //jtable = static_cast<int> (((rho[j]*rho[j]) - tb->innersq) * tb->invdelta);
        jtable = static_cast<int> (((rho[j]*rho[j]) - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
        if (itable >= tlm1 || jtable >= tlm1){
          k_error_flag.template view<DeviceType>()() = 2;
        }
        //A_i = tb->f[itable];
        A_i = d_table_const.f(tidx,itable);
        //A_j = tb->f[jtable];
        A_j = d_table_const.f(tidx,jtable);

        const double rfactor = 1.0-sqrt(rsq/d_cutsq(itype,jtype));
        fpair = 0.5*(A_i + A_j)*(4.0-3.0*rfactor)*rfactor*rfactor*rfactor;
        fpair /= sqrt(rsq);

      } else if (TABSTYLE == LINEAR) {

        //itable = static_cast<int> ((rho[i]*rho[i] - tb->innersq) * tb->invdelta);
        itable = static_cast<int> ((rho[i]*rho[i] - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
        //jtable = static_cast<int> (((rho[j]*rho[j]) - tb->innersq) * tb->invdelta);
        jtable = static_cast<int> ((rho[j]*rho[j] - d_table_const.innersq(tidx)) * d_table_const.invdelta(tidx));
        if (itable >= tlm1 || jtable >= tlm1){
          k_error_flag.template view<DeviceType>()() = 2;
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

        const double rfactor = 1.0-sqrt(rsq/d_cutsq(itype,jtype));
        fpair = 0.5*(A_i + A_j)*(4.0-3.0*rfactor)*rfactor*rfactor*rfactor;
        fpair /= sqrt(rsq);

      } else k_error_flag.template view<DeviceType>()() = 3;

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
      k_error_flag.template view<DeviceType>()() = 2;
    }
    if(itable==0) fraction_i=0.0;
    //else fraction_i = (((rho[i]*rho[i]) - tb->rsq[itable]) * tb->invdelta);
    else fraction_i = (((rho[i]*rho[i]) - d_table_const.rsq(tidx,itable)) * d_table_const.invdelta(tidx));
    //evdwl = tb->e[itable] + fraction_i*tb->de[itable];
    evdwl = d_table_const.e(tidx,itable) + fraction_i*d_table_const.de(tidx,itable);
  } else k_error_flag.template view<DeviceType>()() = 3;

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

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int TABSTYLE>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>(TagPairMultiLucyRXCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG,TABSTYLE>(), ii, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::computeLocalDensity()
{
  x = atomKK->k_x.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  rho = atomKK->k_rho.view<DeviceType>();
  h_rho = atomKK->k_rho.h_view;
  nlocal = atom->nlocal;

  atomKK->sync(execution_space,X_MASK | TYPE_MASK | DPDRHO_MASK);

  const int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  const double pi = MathConst::MY_PI;

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

  atomKK->modified(execution_space,DPDRHO_MASK);

  // communicate and sum densities (on the host)

  if (newton_pair)
    comm->reverse_comm_pair(this);

  comm->forward_comm_pair(this);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXZero, const int &i) const {
  rho[i] = 0.0;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, bool ONE_TYPE>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXComputeLocalDensity<NEIGHFLAG,NEWTON_PAIR,ONE_TYPE>, const int &ii) const {


  // The rho array is atomic for Half/Thread neighbor style
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_rho = rho;

  const int i = d_ilist[ii];

  const double xtmp = x(i,0);
  const double ytmp = x(i,1);
  const double ztmp = x(i,2);

  double rho_i_contrib = 0.0;

  const int itype = type[i];
  const int jnum = d_numneigh[i];

  const double pi = MathConst::MY_PI;

  for (int jj = 0; jj < jnum; jj++){
    const int j = (d_neighbors(i,jj) & NEIGHMASK);
    const int jtype = type[j];

    const double delx = xtmp - x(j,0);
    const double dely = ytmp - x(j,1);
    const double delz = ztmp - x(j,2);
    const double rsq = delx*delx + dely*dely + delz*delz;

    if (ONE_TYPE) {
      if (rsq < cutsq_type11) {
        const double rcut = rcut_type11;
        const double r_over_rcut = sqrt(rsq) / rcut;
        const double tmpFactor = 1.0 - r_over_rcut;
        const double tmpFactor4 = tmpFactor*tmpFactor*tmpFactor*tmpFactor;
        const double factor = factor_type11*(1.0 + 1.5*r_over_rcut)*tmpFactor4;
        rho_i_contrib += factor;
        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
          a_rho[j] += factor;
      }
    } else if (rsq < d_cutsq(itype,jtype)) {
      const double rcut = sqrt(d_cutsq(itype,jtype));
      const double tmpFactor = 1.0-sqrt(rsq)/rcut;
      const double tmpFactor4 = tmpFactor*tmpFactor*tmpFactor*tmpFactor;
      const double factor = (84.0/(5.0*pi*rcut*rcut*rcut))*(1.0+3.0*sqrt(rsq)/(2.0*rcut))*tmpFactor4;
      rho_i_contrib += factor;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
        a_rho[j] += factor;
    }
  }

  a_rho[i] += rho_i_contrib;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::getMixingWeights(int id, double &mixWtSite1old, double &mixWtSite2old, double &mixWtSite1, double &mixWtSite2) const
{
  double fractionOFAold, fractionOFA;
  double fractionOld1, fraction1;
  double fractionOld2, fraction2;
  double nMoleculesOFAold, nMoleculesOFA;
  double nMoleculesOld1, nMolecules1;
  double nMoleculesOld2, nMolecules2;
  double nTotal, nTotalOld;


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

template<class DeviceType>
int PairMultiLucyRXKokkos<DeviceType>::pack_forward_comm_kokkos(int n, DAT::tdual_int_2d k_sendlist, int iswap_in, DAT::tdual_xfloat_1d &buf,
                               int pbc_flag, int *pbc)
{
  atomKK->sync(execution_space,DPDRHO_MASK);

  d_sendlist = k_sendlist.view<DeviceType>();
  iswap = iswap_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagPairMultiLucyRXPackForwardComm>(0,n),*this);
  return n;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXPackForwardComm, const int &i) const {
  int j = d_sendlist(iswap, i);
  v_buf[i] = rho[j];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::unpack_forward_comm_kokkos(int n, int first_in, DAT::tdual_xfloat_1d &buf)
{
  first = first_in;
  v_buf = buf.view<DeviceType>();
  Kokkos::parallel_for(Kokkos::RangePolicy<LMPDeviceType, TagPairMultiLucyRXUnpackForwardComm>(0,n),*this);

  atomKK->modified(execution_space,DPDRHO_MASK);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::operator()(TagPairMultiLucyRXUnpackForwardComm, const int &i) const {
  rho[i + first] = v_buf[i];
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMultiLucyRXKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
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

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) h_rho[i] = buf[m++];

  atomKK->modified(Host,DPDRHO_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
int PairMultiLucyRXKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  atomKK->sync(Host,DPDRHO_MASK);

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = h_rho[i];
  return m;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
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

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairMultiLucyRXKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are atomic for Half/Thread neighbor style
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();

  if (EFLAG) {
    if (eflag_atom) {
      const E_FLOAT epairhalf = 0.5 * epair;
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) v_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) v_eatom[j] += epairhalf;
      } else {
        v_eatom[i] += epairhalf;
      }
    }
  }

  if (VFLAG) {
    const E_FLOAT v0 = delx*delx*fpair;
    const E_FLOAT v1 = dely*dely*fpair;
    const E_FLOAT v2 = delz*delz*fpair;
    const E_FLOAT v3 = delx*dely*fpair;
    const E_FLOAT v4 = delx*delz*fpair;
    const E_FLOAT v5 = dely*delz*fpair;

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

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::create_kokkos_tables()
{
  const int tlm1 = tablength-1;

  memoryKK->create_kokkos(d_table->innersq,h_table->innersq,ntables,"Table::innersq");
  memoryKK->create_kokkos(d_table->invdelta,h_table->invdelta,ntables,"Table::invdelta");

  if(tabstyle == LOOKUP) {
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,tlm1,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,tlm1,"Table::f");
  }

  if(tabstyle == LINEAR) {
    memoryKK->create_kokkos(d_table->rsq,h_table->rsq,ntables,tablength,"Table::rsq");
    memoryKK->create_kokkos(d_table->e,h_table->e,ntables,tablength,"Table::e");
    memoryKK->create_kokkos(d_table->f,h_table->f,ntables,tablength,"Table::f");
    memoryKK->create_kokkos(d_table->de,h_table->de,ntables,tlm1,"Table::de");
    memoryKK->create_kokkos(d_table->df,h_table->df,ntables,tlm1,"Table::df");
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


  Kokkos::deep_copy(d_table->innersq,h_table->innersq);
  Kokkos::deep_copy(d_table->invdelta,h_table->invdelta);
  Kokkos::deep_copy(d_table->rsq,h_table->rsq);
  Kokkos::deep_copy(d_table->e,h_table->e);
  Kokkos::deep_copy(d_table->de,h_table->de);
  Kokkos::deep_copy(d_table->f,h_table->f);
  Kokkos::deep_copy(d_table->df,h_table->df);
  Kokkos::deep_copy(d_table->tabindex,h_table->tabindex);

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

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::allocate()
{
  allocated = 1;
  const int nt = atom->ntypes + 1;

  memory->create(setflag,nt,nt,"pair:setflag");

  memoryKK->create_kokkos(k_cutsq,cutsq,nt,nt,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_cutsq.template modify<LMPHostType>();

  memoryKK->create_kokkos(d_table->tabindex,h_table->tabindex,tabindex,nt,nt,"pair:tabindex");
  d_table_const.tabindex = d_table->tabindex;

  memset(&setflag[0][0],0,nt*nt*sizeof(int));
  memset(&cutsq[0][0],0,nt*nt*sizeof(double));
  memset(&tabindex[0][0],0,nt*nt*sizeof(int));
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairMultiLucyRXKokkos<DeviceType>::settings(int narg, char **arg)
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

    d_table_const.tabindex = d_table->tabindex = typename ArrayTypes<DeviceType>::t_int_2d();
    h_table->tabindex = typename ArrayTypes<LMPHostType>::t_int_2d();
  }
  allocated = 0;

  ntables = 0;
  tables = NULL;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairMultiLucyRXKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class PairMultiLucyRXKokkos<LMPHostType>;
#endif
}
