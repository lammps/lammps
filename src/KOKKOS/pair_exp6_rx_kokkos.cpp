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
   Contributing author: Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_exp6_rx_kokkos.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "fix.h"
#include <float.h>
#include "atom_masks.h"
#include "neigh_request.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

#define MAXLINE 1024
#define DELTA 4

#ifdef DBL_EPSILON
  #define MY_EPSILON (10.0*DBL_EPSILON)
#else
  #define MY_EPSILON (10.0*2.220446049250313e-16)
#endif

#define oneFluidApproxParameter (-1)
#define isOneFluidApprox(_site) ( (_site) == oneFluidApproxParameter )

#define exp6PotentialType (1)
#define isExp6PotentialType(_type) ( (_type) == exp6PotentialType )

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairExp6rxKokkos<DeviceType>::PairExp6rxKokkos(LAMMPS *lmp) : PairExp6rx(lmp)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  k_error_flag = DAT::tdual_int_scalar("pair:error_flag");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairExp6rxKokkos<DeviceType>::~PairExp6rxKokkos()
{
  if (copymode) return;

  memory->destroy_kokkos(k_eatom,eatom);
  memory->destroy_kokkos(k_vatom,vatom);

  memory->destroy_kokkos(k_cutsq,cutsq);

  for (int i=0; i < nparams; ++i) {
    delete[] params[i].name;
    delete[] params[i].potential;
  }
  memory->destroy_kokkos(k_params,params);

  memory->destroy_kokkos(k_mol2param,mol2param);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::init_style()
{
  PairExp6rx::init_style();

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
    error->all(FLERR,"Cannot use chosen neighbor list style with reax/c/kk");
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memory->destroy_kokkos(k_eatom,eatom);
    memory->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.d_view;
  }
  if (vflag_atom) {
    memory->destroy_kokkos(k_vatom,vatom);
    memory->create_kokkos(k_vatom,vatom,maxvatom,6,"pair:vatom");
    d_vatom = k_vatom.d_view;
  }

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  uCG = atomKK->k_uCG.view<DeviceType>();
  uCGnew = atomKK->k_uCGnew.view<DeviceType>();
  dvector = atomKK->k_dvector.view<DeviceType>();
  nlocal = atom->nlocal;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  newton_pair = force->newton_pair;

  atomKK->sync(execution_space,X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK | UCG_MASK | UCGNEW_MASK | DVECTOR_MASK);
  if (evflag) atomKK->modified(execution_space,F_MASK | ENERGY_MASK | VIRIAL_MASK | UCG_MASK | UCGNEW_MASK);
  else atomKK->modified(execution_space,F_MASK | UCG_MASK | UCGNEW_MASK);
  k_cutsq.template sync<DeviceType>();

  // Initialize the Exp6 parameter data for both the local
  // and ghost atoms. Make the parameter data persistent
  // and exchange like any other atom property later.

  {
     const int np_total = nlocal + atom->nghost;

     PairExp6ParamData.epsilon1     = typename AT::t_float_1d("PairExp6ParamData.epsilon1"    ,np_total);
     PairExp6ParamData.alpha1       = typename AT::t_float_1d("PairExp6ParamData.alpha1"      ,np_total);
     PairExp6ParamData.rm1          = typename AT::t_float_1d("PairExp6ParamData.rm1"         ,np_total);
     PairExp6ParamData.mixWtSite1    = typename AT::t_float_1d("PairExp6ParamData.mixWtSite1"   ,np_total);
     PairExp6ParamData.epsilon2     = typename AT::t_float_1d("PairExp6ParamData.epsilon2"    ,np_total);
     PairExp6ParamData.alpha2       = typename AT::t_float_1d("PairExp6ParamData.alpha2"      ,np_total);
     PairExp6ParamData.rm2          = typename AT::t_float_1d("PairExp6ParamData.rm2"         ,np_total);
     PairExp6ParamData.mixWtSite2    = typename AT::t_float_1d("PairExp6ParamData.mixWtSite2"   ,np_total);
     PairExp6ParamData.epsilonOld1  = typename AT::t_float_1d("PairExp6ParamData.epsilonOld1" ,np_total);
     PairExp6ParamData.alphaOld1    = typename AT::t_float_1d("PairExp6ParamData.alphaOld1"   ,np_total);
     PairExp6ParamData.rmOld1       = typename AT::t_float_1d("PairExp6ParamData.rmOld1"      ,np_total);
     PairExp6ParamData.mixWtSite1old = typename AT::t_float_1d("PairExp6ParamData.mixWtSite1old",np_total);
     PairExp6ParamData.epsilonOld2  = typename AT::t_float_1d("PairExp6ParamData.epsilonOld2" ,np_total);
     PairExp6ParamData.alphaOld2    = typename AT::t_float_1d("PairExp6ParamData.alphaOld2"   ,np_total);
     PairExp6ParamData.rmOld2       = typename AT::t_float_1d("PairExp6ParamData.rmOld2"      ,np_total);
     PairExp6ParamData.mixWtSite2old = typename AT::t_float_1d("PairExp6ParamData.mixWtSite2old",np_total);

     Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxgetMixingWeights>(0,np_total),*this);
  }

  k_error_flag.template modify<DeviceType>();
  k_error_flag.template sync<LMPHostType>();
  if (k_error_flag.h_view() == 1)
    error->all(FLERR,"The number of molecules in CG particle is less than 10*DBL_EPSILON.");
  else if (k_error_flag.h_view() == 2)
    error->all(FLERR,"Computed fraction less than -10*DBL_EPSILON");

  int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (neighflag == HALF) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<HALF,1,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<HALF,1,0> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<HALF,0,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<HALF,0,0> >(0,inum),*this);
    }
  } else if (neighflag == HALFTHREAD) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<HALFTHREAD,1,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<HALFTHREAD,1,0> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<HALFTHREAD,0,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<HALFTHREAD,0,0> >(0,inum),*this);
    }
  }

  k_error_flag.template modify<DeviceType>();
  k_error_flag.template sync<LMPHostType>();
  if (k_error_flag.h_view())
    error->all(FLERR,"alpha_ij is 6.0 in pair exp6");

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

  copymode = 0;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::operator()(TagPairExp6rxgetMixingWeights, const int &i) const {
  getMixingWeights (i, PairExp6ParamData.epsilon1[i],
                    PairExp6ParamData.alpha1[i],
                    PairExp6ParamData.rm1[i],
                    PairExp6ParamData.mixWtSite1[i],
                    PairExp6ParamData.epsilon2[i],
                    PairExp6ParamData.alpha2[i],
                    PairExp6ParamData.rm2[i],
                    PairExp6ParamData.mixWtSite2[i],
                    PairExp6ParamData.epsilonOld1[i],
                    PairExp6ParamData.alphaOld1[i],
                    PairExp6ParamData.rmOld1[i],
                    PairExp6ParamData.mixWtSite1old[i],
                    PairExp6ParamData.epsilonOld2[i],
                    PairExp6ParamData.alphaOld2[i],
                    PairExp6ParamData.rmOld2[i],
                    PairExp6ParamData.mixWtSite2old[i]);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::operator()(TagPairExp6rxCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  // These arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_uCG = uCG;
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_uCGnew = uCGnew;

  int i,jj,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwlOld,fpair;
  double rsq,r2inv,r6inv,forceExp6,factor_lj;
  double rCut,rCutInv,rCut2inv,rCut6inv,rCutExp,urc,durc;
  double rm2ij,rm6ij;
  double r,rexp;

  double alphaOld12_ij, rmOld12_ij, epsilonOld12_ij;
  double alphaOld21_ij, rmOld21_ij, epsilonOld21_ij;
  double alpha12_ij, rm12_ij, epsilon12_ij;
  double alpha21_ij, rm21_ij, epsilon21_ij;
  double rminv, buck1, buck2;
  double epsilonOld1_i,alphaOld1_i,rmOld1_i;
  double epsilonOld1_j,alphaOld1_j,rmOld1_j;
  double epsilonOld2_i,alphaOld2_i,rmOld2_i;
  double epsilonOld2_j,alphaOld2_j,rmOld2_j;
  double epsilon1_i,alpha1_i,rm1_i;
  double epsilon1_j,alpha1_j,rm1_j;
  double epsilon2_i,alpha2_i,rm2_i;
  double epsilon2_j,alpha2_j,rm2_j;
  double evdwlOldEXP6_12, evdwlOldEXP6_21, fpairOldEXP6_12, fpairOldEXP6_21;
  double evdwlEXP6_12, evdwlEXP6_21;
  double mixWtSite1old_i, mixWtSite1old_j;
  double mixWtSite2old_i, mixWtSite2old_j;
  double mixWtSite1_i, mixWtSite1_j;
  double mixWtSite2_i, mixWtSite2_j;

  const int nRep = 12;
  const double shift = 1.05;
  double rin1, aRep, uin1, win1, uin1rep, rin1exp, rin6, rin6inv;

  evdwlOld = 0.0;
  evdwl = 0.0;

  i = d_ilist[ii];
  xtmp = x(i,0);
  ytmp = x(i,1);
  ztmp = x(i,2);
  itype = type[i];
  jnum = d_numneigh[i];

  double fx_i = 0.0;
  double fy_i = 0.0;
  double fz_i = 0.0;
  double uCG_i = 0.0;
  double uCGnew_i = 0.0;

  {
     epsilon1_i     = PairExp6ParamData.epsilon1[i];
     alpha1_i       = PairExp6ParamData.alpha1[i];
     rm1_i          = PairExp6ParamData.rm1[i];
     mixWtSite1_i    = PairExp6ParamData.mixWtSite1[i];
     epsilon2_i     = PairExp6ParamData.epsilon2[i];
     alpha2_i       = PairExp6ParamData.alpha2[i];
     rm2_i          = PairExp6ParamData.rm2[i];
     mixWtSite2_i    = PairExp6ParamData.mixWtSite2[i];
     epsilonOld1_i  = PairExp6ParamData.epsilonOld1[i];
     alphaOld1_i    = PairExp6ParamData.alphaOld1[i];
     rmOld1_i       = PairExp6ParamData.rmOld1[i];
     mixWtSite1old_i = PairExp6ParamData.mixWtSite1old[i];
     epsilonOld2_i  = PairExp6ParamData.epsilonOld2[i];
     alphaOld2_i    = PairExp6ParamData.alphaOld2[i];
     rmOld2_i       = PairExp6ParamData.rmOld2[i];
     mixWtSite2old_i = PairExp6ParamData.mixWtSite2old[i];
  }

  for (jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    factor_lj = special_lj[sbmask(j)];
    j &= NEIGHMASK;

    delx = xtmp - x(j,0);
    dely = ytmp - x(j,1);
    delz = ztmp - x(j,2);

    rsq = delx*delx + dely*dely + delz*delz;
    jtype = type[j];

    if (rsq < d_cutsq(itype,jtype)) { // optimize
      r2inv = 1.0/rsq;
      r6inv = r2inv*r2inv*r2inv;

      r = sqrt(rsq);
      rCut2inv = 1.0/d_cutsq(itype,jtype);
      rCut6inv = rCut2inv*rCut2inv*rCut2inv;
      rCut = sqrt(d_cutsq(itype,jtype));
      rCutInv = 1.0/rCut;

      //
      // A. Compute the exp-6 potential
      //

      // A1.  Get alpha, epsilon and rm for particle j

      {
         epsilon1_j     = PairExp6ParamData.epsilon1[j];
         alpha1_j       = PairExp6ParamData.alpha1[j];
         rm1_j          = PairExp6ParamData.rm1[j];
         mixWtSite1_j    = PairExp6ParamData.mixWtSite1[j];
         epsilon2_j     = PairExp6ParamData.epsilon2[j];
         alpha2_j       = PairExp6ParamData.alpha2[j];
         rm2_j          = PairExp6ParamData.rm2[j];
         mixWtSite2_j    = PairExp6ParamData.mixWtSite2[j];
         epsilonOld1_j  = PairExp6ParamData.epsilonOld1[j];
         alphaOld1_j    = PairExp6ParamData.alphaOld1[j];
         rmOld1_j       = PairExp6ParamData.rmOld1[j];
         mixWtSite1old_j = PairExp6ParamData.mixWtSite1old[j];
         epsilonOld2_j  = PairExp6ParamData.epsilonOld2[j];
         alphaOld2_j    = PairExp6ParamData.alphaOld2[j];
         rmOld2_j       = PairExp6ParamData.rmOld2[j];
         mixWtSite2old_j = PairExp6ParamData.mixWtSite2old[j];
      }

      // A2.  Apply Lorentz-Berthelot mixing rules for the i-j pair
      alphaOld12_ij = sqrt(alphaOld1_i*alphaOld2_j);
      rmOld12_ij = 0.5*(rmOld1_i + rmOld2_j);
      epsilonOld12_ij = sqrt(epsilonOld1_i*epsilonOld2_j);
      alphaOld21_ij = sqrt(alphaOld2_i*alphaOld1_j);
      rmOld21_ij = 0.5*(rmOld2_i + rmOld1_j);
      epsilonOld21_ij = sqrt(epsilonOld2_i*epsilonOld1_j);

      alpha12_ij = sqrt(alpha1_i*alpha2_j);
      rm12_ij = 0.5*(rm1_i + rm2_j);
      epsilon12_ij = sqrt(epsilon1_i*epsilon2_j);
      alpha21_ij = sqrt(alpha2_i*alpha1_j);
      rm21_ij = 0.5*(rm2_i + rm1_j);
      epsilon21_ij = sqrt(epsilon2_i*epsilon1_j);

      if(rmOld12_ij!=0.0 && rmOld21_ij!=0.0){
        if(alphaOld21_ij == 6.0 || alphaOld12_ij == 6.0)
          k_error_flag.d_view() = 1;

        // A3.  Compute some convenient quantities for evaluating the force
        rminv = 1.0/rmOld12_ij;
        buck1 = epsilonOld12_ij / (alphaOld12_ij - 6.0);
        rexp = expValue(alphaOld12_ij*(1.0-r*rminv));
        rm2ij = rmOld12_ij*rmOld12_ij;
        rm6ij = rm2ij*rm2ij*rm2ij;

        // Compute the shifted potential
        rCutExp = expValue(alphaOld12_ij*(1.0-rCut*rminv));
        buck2 = 6.0*alphaOld12_ij;
        urc = buck1*(6.0*rCutExp - alphaOld12_ij*rm6ij*rCut6inv);
        durc = -buck1*buck2*(rCutExp* rminv - rCutInv*rm6ij*rCut6inv);
        rin1 = shift*rmOld12_ij*func_rin(alphaOld12_ij);
        if(r < rin1){
          rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
          rin6inv = 1.0/rin6;

          rin1exp = expValue(alphaOld12_ij*(1.0-rin1*rminv));

          uin1 = buck1*(6.0*rin1exp - alphaOld12_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

          win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

          aRep = win1*powint(rin1,nRep)/nRep;

          uin1rep = aRep/powint(rin1,nRep);

          forceExp6 = double(nRep)*aRep/powint(r,nRep);
          fpairOldEXP6_12 = factor_lj*forceExp6*r2inv;

          evdwlOldEXP6_12 = uin1 - uin1rep + aRep/powint(r,nRep);
        } else {
          forceExp6 = buck1*buck2*(r*rexp*rminv - rm6ij*r6inv) + r*durc;
          fpairOldEXP6_12 = factor_lj*forceExp6*r2inv;

          evdwlOldEXP6_12 = buck1*(6.0*rexp - alphaOld12_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
        }

        // A3.  Compute some convenient quantities for evaluating the force
        rminv = 1.0/rmOld21_ij;
        buck1 = epsilonOld21_ij / (alphaOld21_ij - 6.0);
        buck2 = 6.0*alphaOld21_ij;
        rexp = expValue(alphaOld21_ij*(1.0-r*rminv));
        rm2ij = rmOld21_ij*rmOld21_ij;
        rm6ij = rm2ij*rm2ij*rm2ij;

        // Compute the shifted potential
        rCutExp = expValue(alphaOld21_ij*(1.0-rCut*rminv));
        buck2 = 6.0*alphaOld21_ij;
        urc = buck1*(6.0*rCutExp - alphaOld21_ij*rm6ij*rCut6inv);
        durc = -buck1*buck2*(rCutExp* rminv - rCutInv*rm6ij*rCut6inv);
        rin1 = shift*rmOld21_ij*func_rin(alphaOld21_ij);

        if(r < rin1){
          rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
          rin6inv = 1.0/rin6;

          rin1exp = expValue(alphaOld21_ij*(1.0-rin1*rminv));

          uin1 = buck1*(6.0*rin1exp - alphaOld21_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

          win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

          aRep = win1*powint(rin1,nRep)/nRep;

          uin1rep = aRep/powint(rin1,nRep);

          forceExp6 = double(nRep)*aRep/powint(r,nRep);
          fpairOldEXP6_21 = factor_lj*forceExp6*r2inv;

          evdwlOldEXP6_21 = uin1 - uin1rep + aRep/powint(r,nRep);
        } else {
          forceExp6 = buck1*buck2*(r*rexp*rminv - rm6ij*r6inv) + r*durc;
          fpairOldEXP6_21 = factor_lj*forceExp6*r2inv;

          evdwlOldEXP6_21 = buck1*(6.0*rexp - alphaOld21_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
        }

        if (isite1 == isite2)
          evdwlOld = sqrt(mixWtSite1old_i*mixWtSite2old_j)*evdwlOldEXP6_12;
        else
          evdwlOld = sqrt(mixWtSite1old_i*mixWtSite2old_j)*evdwlOldEXP6_12 + sqrt(mixWtSite2old_i*mixWtSite1old_j)*evdwlOldEXP6_21;

        evdwlOld *= factor_lj;

        uCG_i += 0.5*evdwlOld;
        if (NEWTON_PAIR || j < nlocal)
          a_uCG[j] += 0.5*evdwlOld;
      }

      if(rm12_ij!=0.0 && rm21_ij!=0.0){
        if(alpha21_ij == 6.0 || alpha12_ij == 6.0)
          k_error_flag.d_view() = 1;

        // A3.  Compute some convenient quantities for evaluating the force
        rminv = 1.0/rm12_ij;
        buck1 = epsilon12_ij / (alpha12_ij - 6.0);
        buck2 = 6.0*alpha12_ij;
        rexp = expValue(alpha12_ij*(1.0-r*rminv));
        rm2ij = rm12_ij*rm12_ij;
        rm6ij = rm2ij*rm2ij*rm2ij;

        // Compute the shifted potential
        rCutExp = expValue(alpha12_ij*(1.0-rCut*rminv));
        urc = buck1*(6.0*rCutExp - alpha12_ij*rm6ij*rCut6inv);
        durc = -buck1*buck2*(rCutExp*rminv - rCutInv*rm6ij*rCut6inv);
        rin1 = shift*rm12_ij*func_rin(alpha12_ij);

        if(r < rin1){
          rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
          rin6inv = 1.0/rin6;

          rin1exp = expValue(alpha12_ij*(1.0-rin1*rminv));

          uin1 = buck1*(6.0*rin1exp - alpha12_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

          win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

          aRep = win1*powint(rin1,nRep)/nRep;

          uin1rep = aRep/powint(rin1,nRep);

          evdwlEXP6_12 = uin1 - uin1rep + aRep/powint(r,nRep);
        } else {
          evdwlEXP6_12 = buck1*(6.0*rexp - alpha12_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
        }

        rminv = 1.0/rm21_ij;
        buck1 = epsilon21_ij / (alpha21_ij - 6.0);
        buck2 = 6.0*alpha21_ij;
        rexp = expValue(alpha21_ij*(1.0-r*rminv));
        rm2ij = rm21_ij*rm21_ij;
        rm6ij = rm2ij*rm2ij*rm2ij;

        // Compute the shifted potential
        rCutExp = expValue(alpha21_ij*(1.0-rCut*rminv));
        urc = buck1*(6.0*rCutExp - alpha21_ij*rm6ij*rCut6inv);
        durc = -buck1*buck2*(rCutExp*rminv - rCutInv*rm6ij*rCut6inv);
        rin1 = shift*rm21_ij*func_rin(alpha21_ij);

        if(r < rin1){
          rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
          rin6inv = 1.0/rin6;

          rin1exp = expValue(alpha21_ij*(1.0-rin1*rminv));

          uin1 = buck1*(6.0*rin1exp - alpha21_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

          win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

          aRep = win1*powint(rin1,nRep)/nRep;

          uin1rep = aRep/powint(rin1,nRep);

          evdwlEXP6_21 = uin1 - uin1rep + aRep/powint(r,nRep);
        } else {
          evdwlEXP6_21 = buck1*(6.0*rexp - alpha21_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
        }

        //
        // Apply Mixing Rule to get the overall force for the CG pair
        //
        if (isite1 == isite2) fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12;
        else fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12 + sqrt(mixWtSite2old_i*mixWtSite1old_j)*fpairOldEXP6_21;

        fx_i += delx*fpair;
        fy_i += dely*fpair;
        fz_i += delz*fpair;
        if (NEWTON_PAIR || j < nlocal) {
          a_f(j,0) -= delx*fpair;
          a_f(j,1) -= dely*fpair;
          a_f(j,2) -= delz*fpair;
        }

        if (isite1 == isite2) evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12;
        else evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12 + sqrt(mixWtSite2_i*mixWtSite1_j)*evdwlEXP6_21;
        evdwl *= factor_lj;

        uCGnew_i   += 0.5*evdwl;
        if (NEWTON_PAIR || j < nlocal)
          a_uCGnew[j] += 0.5*evdwl;
        evdwl = evdwlOld;
        if (EVFLAG)
          ev.evdwl += ((NEWTON_PAIR||(j<nlocal))?1.0:0.5)*evdwl;
        //if (vflag_either || eflag_atom) 
        if (EVFLAG) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,evdwl,fpair,delx,dely,delz);
      }
    }
  }

  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
  a_uCG[i] += uCG_i;
  a_uCGnew[i] += uCGnew_i;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::operator()(TagPairExp6rxCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairExp6rxCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_cutsq.template modify<LMPHostType>();

  memory->create(cut,n+1,n+1,"pair:cut_lj");
}


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairExp6rx::coeff(narg,arg);
  
  if (scalingFlag == POLYNOMIAL)
    for (int i = 0; i < 6; i++) {
      s_coeffAlpha[i] = coeffAlpha[i];
      s_coeffEps[i] = coeffEps[i];
      s_coeffRm[i] = coeffRm[i];
    }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::read_file(char *file)
{
  int params_per_line = 5;
  char **words = new char*[params_per_line+1];

  memory->destroy_kokkos(k_params,params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  fp = NULL;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open exp6/rx potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each set of params from potential file
  // one set of params can span multiple lines

  int n,nwords,ispecies;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in exp6/rx potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    for (ispecies = 0; ispecies < nspecies; ispecies++)
      if (strcmp(words[0],&atom->dname[ispecies][0]) == 0) break;
    if (ispecies == nspecies) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      k_params.template modify<LMPHostType>();
      maxparam += DELTA;
      memory->grow_kokkos(k_params,params,maxparam,
                          "pair:params");
    }

    params[nparams].ispecies = ispecies;

    n = strlen(&atom->dname[ispecies][0]) + 1;
    params[nparams].name = new char[n];
    strcpy(params[nparams].name,&atom->dname[ispecies][0]);

    n = strlen(words[1]) + 1;
    params[nparams].potential = new char[n];
    strcpy(params[nparams].potential,words[1]);
    if (strcmp(params[nparams].potential,"exp6") == 0){
      params[nparams].alpha = atof(words[2]);
      params[nparams].epsilon = atof(words[3]);
      params[nparams].rm = atof(words[4]);
      if (params[nparams].epsilon <= 0.0 || params[nparams].rm <= 0.0 ||
          params[nparams].alpha < 0.0)
        error->all(FLERR,"Illegal exp6/rx parameters.  Rm and Epsilon must be greater than zero.  Alpha cannot be negative.");
    } else {
      error->all(FLERR,"Illegal exp6/rx parameters.  Interaction potential does not exist.");
    }
    nparams++;
  }

  delete [] words;

  k_params.template modify<LMPHostType>();
  k_params.template sync<DeviceType>();
  d_params = k_params.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::setup()
{
  int i,j,n;

  // set mol2param for all combinations
  // must be a single exact match to lines read from file

  memory->destroy_kokkos(k_mol2param,mol2param);
  memory->create_kokkos(k_mol2param,mol2param,nspecies,"pair:mol2param");

  for (i = 0; i < nspecies; i++) {
    n = -1;
    for (j = 0; j < nparams; j++) {
      if (i == params[j].ispecies) {
        if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
        n = j;
      }
    }
    mol2param[i] = n;
  }

  k_mol2param.template modify<LMPHostType>();
  k_mol2param.template sync<DeviceType>();
  d_mol2param = k_mol2param.template view<DeviceType>();

  neighflag = lmp->kokkos->neighflag;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::getMixingWeights(int id,double &epsilon1,double &alpha1,double &rm1, double &mixWtSite1,double &epsilon2,double &alpha2,double &rm2,double &mixWtSite2,double &epsilon1_old,double &alpha1_old,double &rm1_old, double &mixWtSite1old,double &epsilon2_old,double &alpha2_old,double &rm2_old,double &mixWtSite2old) const
{
  int iparam, jparam;
  double rmi, rmj, rmij, rm3ij;
  double epsiloni, epsilonj, epsilonij;
  double alphai, alphaj, alphaij;
  double epsilon_old, rm3_old, alpha_old;
  double epsilon, rm3, alpha;
  double xMolei, xMolej, xMolei_old, xMolej_old;

  double fractionOFAold, fractionOFA;
  double fractionOld1, fraction1;
  double fractionOld2, fraction2;
  double nMoleculesOFAold, nMoleculesOFA;
  double nMoleculesOld1, nMolecules1;
  double nMoleculesOld2, nMolecules2;
  double nTotal, nTotalold;

  rm3 = 0.0;
  epsilon = 0.0;
  alpha = 0.0;
  epsilon_old = 0.0;
  rm3_old = 0.0;
  alpha_old = 0.0;
  fractionOFA = 0.0;
  fractionOFAold = 0.0;
  nMoleculesOFA = 0.0;
  nMoleculesOFAold = 0.0;
  nTotal = 0.0;
  nTotalold = 0.0;

  // Compute the total number of molecules in the old and new CG particle as well as the total number of molecules in the fluid portion of the old and new CG particle
  for (int ispecies = 0; ispecies < nspecies; ispecies++){
    nTotal += dvector(ispecies,id);
    nTotalold += dvector(ispecies+nspecies,id);

    iparam = d_mol2param[ispecies];

    if (iparam < 0 || d_params[iparam].potentialType != exp6PotentialType ) continue;
    if (isOneFluidApprox(isite1) || isOneFluidApprox(isite2)) {
      if (isite1 == d_params[iparam].ispecies || isite2 == d_params[iparam].ispecies) continue;
      nMoleculesOFAold += dvector(ispecies+nspecies,id);
      nMoleculesOFA += dvector(ispecies,id);
    }
  }
  if(nTotal < MY_EPSILON || nTotalold < MY_EPSILON)
    k_error_flag.d_view() = 1;

  // Compute the mole fraction of molecules within the fluid portion of the particle (One Fluid Approximation)
  fractionOFAold = nMoleculesOFAold / nTotalold;
  fractionOFA = nMoleculesOFA / nTotal;

  for (int ispecies = 0; ispecies < nspecies; ispecies++) {
    iparam = d_mol2param[ispecies];
    if (iparam < 0 || d_params[iparam].potentialType != exp6PotentialType ) continue;

    // If Site1 matches a pure species, then grab the parameters
    if (isite1 == d_params[iparam].ispecies){
      rm1_old = d_params[iparam].rm;
      rm1 = d_params[iparam].rm;
      epsilon1_old = d_params[iparam].epsilon;
      epsilon1 = d_params[iparam].epsilon;
      alpha1_old = d_params[iparam].alpha;
      alpha1 = d_params[iparam].alpha;

      // Compute the mole fraction of Site1
      nMoleculesOld1 = dvector(ispecies+nspecies,id);
      nMolecules1 = dvector(ispecies,id);
      fractionOld1 = nMoleculesOld1/nTotalold;
      fraction1 = nMolecules1/nTotal;
    }

    // If Site2 matches a pure species, then grab the parameters
    if (isite2 == d_params[iparam].ispecies){
      rm2_old = d_params[iparam].rm;
      rm2 = d_params[iparam].rm;
      epsilon2_old = d_params[iparam].epsilon;
      epsilon2 = d_params[iparam].epsilon;
      alpha2_old = d_params[iparam].alpha;
      alpha2 = d_params[iparam].alpha;

      // Compute the mole fraction of Site2
      nMoleculesOld2 = dvector(ispecies+nspecies,id);
      nMolecules2 = dvector(ispecies,id);
      fractionOld2 = dvector(ispecies+nspecies,id)/nTotalold;
    }

    // If Site1 or Site2 matches is a fluid, then compute the paramters
    if (isOneFluidApprox(isite1) || isOneFluidApprox(isite2)) {
      if (isite1 == d_params[iparam].ispecies || isite2 == d_params[iparam].ispecies) continue;
      rmi = d_params[iparam].rm;
      epsiloni = d_params[iparam].epsilon;
      alphai = d_params[iparam].alpha;
      if(nMoleculesOFA<MY_EPSILON) xMolei = 0.0;
      else xMolei = dvector(ispecies,id)/nMoleculesOFA;
      if(nMoleculesOFAold<MY_EPSILON) xMolei_old = 0.0;
      else xMolei_old = dvector(ispecies+nspecies,id)/nMoleculesOFAold;

      for (int jspecies = 0; jspecies < nspecies; jspecies++) {
        jparam = d_mol2param[jspecies];
        if (jparam < 0 || d_params[jparam].potentialType != exp6PotentialType ) continue;
        if (isite1 == d_params[jparam].ispecies || isite2 == d_params[jparam].ispecies) continue;
        rmj = d_params[jparam].rm;
        epsilonj = d_params[jparam].epsilon;
        alphaj = d_params[jparam].alpha;
        if(nMoleculesOFA<MY_EPSILON) xMolej = 0.0;
        else xMolej = dvector(jspecies,id)/nMoleculesOFA;
        if(nMoleculesOFAold<MY_EPSILON) xMolej_old = 0.0;
        else xMolej_old = dvector(jspecies+nspecies,id)/nMoleculesOFAold;

        rmij = (rmi+rmj)/2.0;
        rm3ij = rmij*rmij*rmij;
        epsilonij = sqrt(epsiloni*epsilonj);
        alphaij = sqrt(alphai*alphaj);

        if(fractionOFAold > 0.0){
          rm3_old += xMolei_old*xMolej_old*rm3ij;
          epsilon_old += xMolei_old*xMolej_old*rm3ij*epsilonij;
          alpha_old += xMolei_old*xMolej_old*rm3ij*epsilonij*alphaij;
        }
        if(fractionOFA > 0.0){
          rm3 += xMolei*xMolej*rm3ij;
          epsilon += xMolei*xMolej*rm3ij*epsilonij;
          alpha += xMolei*xMolej*rm3ij*epsilonij*alphaij;
        }
      }
    }
  }

  if (isOneFluidApprox(isite1)){
    rm1 = cbrt(rm3);
    if(rm1 < MY_EPSILON) {
      rm1 = 0.0;
      epsilon1 = 0.0;
      alpha1 = 0.0;
    } else {
      epsilon1 = epsilon / rm3;
      alpha1 = alpha / epsilon1 / rm3;
    }
    nMolecules1 = 1.0-(nTotal-nMoleculesOFA);
    fraction1 = fractionOFA;

    rm1_old = cbrt(rm3_old);
    if(rm1_old < MY_EPSILON) {
      rm1_old = 0.0;
      epsilon1_old = 0.0;
      alpha1_old = 0.0;
    } else {
      epsilon1_old = epsilon_old / rm3_old;
      alpha1_old = alpha_old / epsilon1_old / rm3_old;
    }
    nMoleculesOld1 = 1.0-(nTotalold-nMoleculesOFAold);
    fractionOld1 = fractionOFAold;

    if(scalingFlag == EXPONENT){
      exponentScaling(nMoleculesOFA,epsilon1,rm1);
      exponentScaling(nMoleculesOFAold,epsilon1_old,rm1_old);
    } else if(scalingFlag == POLYNOMIAL){
      polynomialScaling(nMoleculesOFA,alpha1,epsilon1,rm1);
      polynomialScaling(nMoleculesOFAold,alpha1_old,epsilon1_old,rm1_old);
    }
  }

  if (isOneFluidApprox(isite2)){
    rm2 = cbrt(rm3);
    if(rm2 < MY_EPSILON) {
      rm2 = 0.0;
      epsilon2 = 0.0;
      alpha2 = 0.0;
    } else {
      epsilon2 = epsilon / rm3;
      alpha2 = alpha / epsilon2 / rm3;
    }
    nMolecules2 = 1.0-(nTotal-nMoleculesOFA);
    fraction2 = fractionOFA;

    rm2_old = cbrt(rm3_old);
    if(rm2_old < MY_EPSILON) {
      rm2_old = 0.0;
      epsilon2_old = 0.0;
      alpha2_old = 0.0;
    } else {
      epsilon2_old = epsilon_old / rm3_old;
      alpha2_old = alpha_old / epsilon2_old / rm3_old;
    }
    nMoleculesOld2 = 1.0-(nTotalold-nMoleculesOFAold);
    fractionOld2 = fractionOFAold;

    if(scalingFlag == EXPONENT){
      exponentScaling(nMoleculesOFA,epsilon2,rm2);
      exponentScaling(nMoleculesOFAold,epsilon2_old,rm2_old);
    } else if(scalingFlag == POLYNOMIAL){
      polynomialScaling(nMoleculesOFA,alpha2,epsilon2,rm2);
      polynomialScaling(nMoleculesOFAold,alpha2_old,epsilon2_old,rm2_old);
    }
  }

  // Check that no fractions are less than zero
  if(fraction1 < 0.0 || nMolecules1 < 0.0){
    if(fraction1 < -MY_EPSILON || nMolecules1 < -MY_EPSILON){
      k_error_flag.d_view() = 2;
    }
    nMolecules1 = 0.0;
    fraction1 = 0.0;
  }
  if(fraction2 < 0.0 || nMolecules2 < 0.0){
    if(fraction2 < -MY_EPSILON || nMolecules2 < -MY_EPSILON){
      k_error_flag.d_view() = 2;
    }
    nMolecules2 = 0.0;
    fraction2 = 0.0;
  }
  if(fractionOld1 < 0.0 || nMoleculesOld1 < 0.0){
    if(fractionOld1 < -MY_EPSILON || nMoleculesOld1 < -MY_EPSILON){
      k_error_flag.d_view() = 2;
    }
    nMoleculesOld1 = 0.0;
    fractionOld1 = 0.0;
  }
  if(fractionOld2 < 0.0 || nMoleculesOld2 < 0.0){
    if(fractionOld2 < -MY_EPSILON || nMoleculesOld2 < -MY_EPSILON){
      k_error_flag.d_view() = 2;
    }
    nMoleculesOld2 = 0.0;
    fractionOld2 = 0.0;
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
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::exponentScaling(double phi, double &epsilon, double &rm) const
{
  double powfuch;

  if(exponentEpsilon < 0.0){
    powfuch = pow(phi,-exponentEpsilon);
    if(powfuch<MY_EPSILON) epsilon = 0.0;
    else epsilon *= 1.0/powfuch;
  } else {
    epsilon *= pow(phi,exponentEpsilon);
  }

  if(exponentR < 0.0){
    powfuch = pow(phi,-exponentR);
    if(powfuch<MY_EPSILON) rm = 0.0;
    else rm *= 1.0/powfuch;
  } else {
    rm *= pow(phi,exponentR);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::polynomialScaling(double phi, double &alpha, double &epsilon, double &rm) const
{
    double phi2 = phi*phi;
    double phi3 = phi2*phi;
    double phi4 = phi2*phi2;
    double phi5 = phi2*phi3;

    alpha = (s_coeffAlpha[0]*phi5 + s_coeffAlpha[1]*phi4 + s_coeffAlpha[2]*phi3 + s_coeffAlpha[3]*phi2 + s_coeffAlpha[4]*phi + s_coeffAlpha[5]);
    epsilon *= (s_coeffEps[0]*phi5 + s_coeffEps[1]*phi4 + s_coeffEps[2]*phi3 + s_coeffEps[3]*phi2 + s_coeffEps[4]*phi + s_coeffEps[5]);
    rm *= (s_coeffEps[0]*phi5 + s_coeffEps[1]*phi4 + s_coeffEps[2]*phi3 + s_coeffEps[3]*phi2 + s_coeffEps[4]*phi + s_coeffEps[5]);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairExp6rxKokkos<DeviceType>::func_rin(const double &alpha) const
{
  double function;

  const double a = 3.7682065;
  const double b = -1.4308614;

  function = a+b*sqrt(alpha);
  function = expValue(function);

  return function;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
double PairExp6rxKokkos<DeviceType>::expValue(double value) const
{
  double returnValue;
  if(value < DBL_MIN_EXP) returnValue = 0.0;
  else returnValue = exp(value);

  return returnValue;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
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
KOKKOS_INLINE_FUNCTION
int PairExp6rxKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairExp6rxKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class PairExp6rxKokkos<LMPHostType>;
#endif
}