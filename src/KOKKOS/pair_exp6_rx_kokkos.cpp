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

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include "pair_exp6_rx_kokkos.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_special_kokkos.h"
#include "memory_kokkos.h"
#include "error.h"
#include "fix.h"
#include "atom_masks.h"
#include "neigh_request.h"
#include "atom_kokkos.h"
#include "kokkos.h"

#include <cmath>
#include <cfloat>
#include <cstring>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace MathSpecialKokkos;

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
  kokkosable = 1;
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

  memoryKK->destroy_kokkos(k_eatom,eatom);
  memoryKK->destroy_kokkos(k_vatom,vatom);

  memoryKK->destroy_kokkos(k_cutsq,cutsq);

  for (int i=0; i < nparams; ++i) {
    delete[] params[i].name;
    delete[] params[i].potential;
  }
  memoryKK->destroy_kokkos(k_params,params);

  memoryKK->destroy_kokkos(k_mol2param,mol2param);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::init_style()
{
  PairExp6rx::init_style();

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
  if (neighflag == FULL) request->enable_full();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

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
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
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

     if (np_total > (int)PairExp6ParamData.epsilon1.extent(0)) {
       PairExp6ParamData.epsilon1      = typename AT::t_float_1d("PairExp6ParamData.epsilon1"     ,np_total);
       PairExp6ParamData.alpha1        = typename AT::t_float_1d("PairExp6ParamData.alpha1"       ,np_total);
       PairExp6ParamData.rm1           = typename AT::t_float_1d("PairExp6ParamData.rm1"          ,np_total);
       PairExp6ParamData.mixWtSite1    = typename AT::t_float_1d("PairExp6ParamData.mixWtSite1"   ,np_total);
       PairExp6ParamData.epsilon2      = typename AT::t_float_1d("PairExp6ParamData.epsilon2"     ,np_total);
       PairExp6ParamData.alpha2        = typename AT::t_float_1d("PairExp6ParamData.alpha2"       ,np_total);
       PairExp6ParamData.rm2           = typename AT::t_float_1d("PairExp6ParamData.rm2"          ,np_total);
       PairExp6ParamData.mixWtSite2    = typename AT::t_float_1d("PairExp6ParamData.mixWtSite2"   ,np_total);
       PairExp6ParamData.epsilonOld1   = typename AT::t_float_1d("PairExp6ParamData.epsilonOld1"  ,np_total);
       PairExp6ParamData.alphaOld1     = typename AT::t_float_1d("PairExp6ParamData.alphaOld1"    ,np_total);
       PairExp6ParamData.rmOld1        = typename AT::t_float_1d("PairExp6ParamData.rmOld1"       ,np_total);
       PairExp6ParamData.mixWtSite1old = typename AT::t_float_1d("PairExp6ParamData.mixWtSite1old",np_total);
       PairExp6ParamData.epsilonOld2   = typename AT::t_float_1d("PairExp6ParamData.epsilonOld2"  ,np_total);
       PairExp6ParamData.alphaOld2     = typename AT::t_float_1d("PairExp6ParamData.alphaOld2"    ,np_total);
       PairExp6ParamData.rmOld2        = typename AT::t_float_1d("PairExp6ParamData.rmOld2"       ,np_total);
       PairExp6ParamData.mixWtSite2old = typename AT::t_float_1d("PairExp6ParamData.mixWtSite2old",np_total);

       PairExp6ParamDataVect.epsilon          = typename AT::t_float_1d("PairExp6ParamDataVect.epsilon"         ,np_total);
       PairExp6ParamDataVect.rm3              = typename AT::t_float_1d("PairExp6ParamDataVect.rm3"             ,np_total);
       PairExp6ParamDataVect.alpha            = typename AT::t_float_1d("PairExp6ParamDataVect.alpha"           ,np_total);
       PairExp6ParamDataVect.xMolei           = typename AT::t_float_1d("PairExp6ParamDataVect.xMolei"          ,np_total);
       PairExp6ParamDataVect.epsilon_old      = typename AT::t_float_1d("PairExp6ParamDataVect.epsilon_old"     ,np_total);
       PairExp6ParamDataVect.rm3_old          = typename AT::t_float_1d("PairExp6ParamDataVect.rm3_old"         ,np_total);
       PairExp6ParamDataVect.alpha_old        = typename AT::t_float_1d("PairExp6ParamDataVect.alpha_old"       ,np_total);
       PairExp6ParamDataVect.xMolei_old       = typename AT::t_float_1d("PairExp6ParamDataVect.xMolei_old"      ,np_total);
       PairExp6ParamDataVect.fractionOFA      = typename AT::t_float_1d("PairExp6ParamDataVect.fractionOFA"     ,np_total);
       PairExp6ParamDataVect.fraction1        = typename AT::t_float_1d("PairExp6ParamDataVect.fraction1"       ,np_total);
       PairExp6ParamDataVect.fraction2        = typename AT::t_float_1d("PairExp6ParamDataVect.fraction2"       ,np_total);
       PairExp6ParamDataVect.nMoleculesOFA    = typename AT::t_float_1d("PairExp6ParamDataVect.nMoleculesOFA"   ,np_total);
       PairExp6ParamDataVect.nMolecules1      = typename AT::t_float_1d("PairExp6ParamDataVect.nMolecules1"     ,np_total);
       PairExp6ParamDataVect.nMolecules2      = typename AT::t_float_1d("PairExp6ParamDataVect.nMolecules2"     ,np_total);
       PairExp6ParamDataVect.nTotal           = typename AT::t_float_1d("PairExp6ParamDataVect.nTotal"          ,np_total);
       PairExp6ParamDataVect.fractionOFAold   = typename AT::t_float_1d("PairExp6ParamDataVect.fractionOFAold"  ,np_total);
       PairExp6ParamDataVect.fractionOld1     = typename AT::t_float_1d("PairExp6ParamDataVect.fractionOld1"    ,np_total);
       PairExp6ParamDataVect.fractionOld2     = typename AT::t_float_1d("PairExp6ParamDataVect.fractionOld2"    ,np_total);
       PairExp6ParamDataVect.nMoleculesOFAold = typename AT::t_float_1d("PairExp6ParamDataVect.nMoleculesOFAold",np_total);
       PairExp6ParamDataVect.nMoleculesOld1   = typename AT::t_float_1d("PairExp6ParamDataVect.nMoleculesOld1"  ,np_total);
       PairExp6ParamDataVect.nMoleculesOld2   = typename AT::t_float_1d("PairExp6ParamDataVect.nMoleculesOld2"  ,np_total);
       PairExp6ParamDataVect.nTotalold        = typename AT::t_float_1d("PairExp6ParamDataVect.nTotalold"       ,np_total);
     } else
       Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxZeroMixingWeights>(0,np_total),*this);

#ifdef LMP_KOKKOS_GPU
     Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxgetMixingWeights>(0,np_total),*this);
#else
     int errorFlag = 0;
     getMixingWeightsVect (np_total, errorFlag, PairExp6ParamData.epsilon1,
                                                PairExp6ParamData.alpha1,
                                                PairExp6ParamData.rm1,
                                                PairExp6ParamData.mixWtSite1,
                                                PairExp6ParamData.epsilon2,
                                                PairExp6ParamData.alpha2,
                                                PairExp6ParamData.rm2,
                                                PairExp6ParamData.mixWtSite2,
                                                PairExp6ParamData.epsilonOld1,
                                                PairExp6ParamData.alphaOld1,
                                                PairExp6ParamData.rmOld1,
                                                PairExp6ParamData.mixWtSite1old,
                                                PairExp6ParamData.epsilonOld2,
                                                PairExp6ParamData.alphaOld2,
                                                PairExp6ParamData.rmOld2,
                                                PairExp6ParamData.mixWtSite2old);
     if (errorFlag == 1)
       error->all(FLERR,"The number of molecules in CG particle is less than 10*DBL_EPSILON.");
     else if (errorFlag == 2)
       error->all(FLERR,"Computed fraction less than -10*DBL_EPSILON");
#endif
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

#ifdef LMP_KOKKOS_GPU  // Use atomics

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
  } else if (neighflag == FULL) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<FULL,1,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<FULL,1,0> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<FULL,0,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCompute<FULL,0,0> >(0,inum),*this);
    }
  }

#else // No atomics

  nthreads = lmp->kokkos->nthreads;
  int nmax = f.extent(0);
  if (nmax > (int)t_f.extent(1)) {
    t_f = t_f_array_thread("pair_exp6_rx:t_f",nthreads,nmax);
    t_uCG = t_efloat_1d_thread("pair_exp6_rx:t_uCG",nthreads,nmax);
    t_uCGnew = t_efloat_1d_thread("pair_exp6_rx:t_UCGnew",nthreads,nmax);
  }

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxZeroDupViews>(0,nmax),*this);

  if (neighflag == HALF) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<HALF,1,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<HALF,1,0> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<HALF,0,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<HALF,0,0> >(0,inum),*this);
    }
  } else if (neighflag == HALFTHREAD) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<HALFTHREAD,1,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<HALFTHREAD,1,0> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<HALFTHREAD,0,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<HALFTHREAD,0,0> >(0,inum),*this);
    }
  } else if (neighflag == FULL) {
    if (newton_pair) {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<FULL,1,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<FULL,1,0> >(0,inum),*this);
    } else {
      if (evflag) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<FULL,0,1> >(0,inum),*this,ev);
      else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxComputeNoAtomics<FULL,0,0> >(0,inum),*this);
    }
  }

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairExp6rxCollapseDupViews>(0,nmax),*this);

#endif

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
void PairExp6rxKokkos<DeviceType>::operator()(TagPairExp6rxZeroMixingWeights, const int &i) const {
  PairExp6ParamData.epsilon1[i] = 0.0;
  PairExp6ParamData.alpha1[i] = 0.0;
  PairExp6ParamData.rm1[i] = 0.0;
  PairExp6ParamData.mixWtSite1[i] = 0.0;
  PairExp6ParamData.epsilon2[i] = 0.0;
  PairExp6ParamData.alpha2[i] = 0.0;
  PairExp6ParamData.rm2[i] = 0.0;
  PairExp6ParamData.mixWtSite2[i] = 0.0;
  PairExp6ParamData.epsilonOld1[i] = 0.0;
  PairExp6ParamData.alphaOld1[i] = 0.0;
  PairExp6ParamData.rmOld1[i] = 0.0;
  PairExp6ParamData.mixWtSite1old[i] = 0.0;
  PairExp6ParamData.epsilonOld2[i] = 0.0;
  PairExp6ParamData.alphaOld2[i] = 0.0;
  PairExp6ParamData.rmOld2[i] = 0.0;
  PairExp6ParamData.mixWtSite2old[i] = 0.0;
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

  {
    const bool one_type = (ntypes == 1);
    if (isite1 == isite2)
      if (one_type)
        this->vectorized_operator<NEIGHFLAG,NEWTON_PAIR,EVFLAG,true, true, true>(ii, ev);
      else
        this->vectorized_operator<NEIGHFLAG,NEWTON_PAIR,EVFLAG,true, true,false>(ii, ev);
    else
      if (one_type)
        this->vectorized_operator<NEIGHFLAG,NEWTON_PAIR,EVFLAG,false,true, true>(ii, ev);
      else
        this->vectorized_operator<NEIGHFLAG,NEWTON_PAIR,EVFLAG,false,true,false>(ii, ev);
    return;
  }

  // These arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_uCG = uCG;
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_uCGnew = uCGnew;

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

      evdwlOldEXP6_12 = 0.0;
      evdwlOldEXP6_21 = 0.0;
      evdwlEXP6_12 = 0.0;
      evdwlEXP6_21 = 0.0;
      fpairOldEXP6_12 = 0.0;
      fpairOldEXP6_21 = 0.0;

      if (rmOld12_ij!=0.0 && rmOld21_ij!=0.0) {
        if (alphaOld21_ij == 6.0 || alphaOld12_ij == 6.0)
          k_error_flag.template view<DeviceType>()() = 1;

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
        if (r < rin1) {
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

        if (r < rin1) {
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
        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
          a_uCG[j] += 0.5*evdwlOld;
      }

      if (rm12_ij!=0.0 && rm21_ij!=0.0) {
        if (alpha21_ij == 6.0 || alpha12_ij == 6.0)
          k_error_flag.template view<DeviceType>()() = 1;

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

        if (r < rin1) {
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

        if (r < rin1) {
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
      }

      //
      // Apply Mixing Rule to get the overall force for the CG pair
      //
      if (isite1 == isite2) fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12;
      else fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12 + sqrt(mixWtSite2old_i*mixWtSite1old_j)*fpairOldEXP6_21;

      fx_i += delx*fpair;
      fy_i += dely*fpair;
      fz_i += delz*fpair;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        a_f(j,0) -= delx*fpair;
        a_f(j,1) -= dely*fpair;
        a_f(j,2) -= delz*fpair;
      }

      if (isite1 == isite2) evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12;
      else evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12 + sqrt(mixWtSite2_i*mixWtSite1_j)*evdwlEXP6_21;
      evdwl *= factor_lj;

      uCGnew_i   += 0.5*evdwl;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
        a_uCGnew[j] += 0.5*evdwl;
      evdwl = evdwlOld;
      if (EVFLAG)
        ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*evdwl;
      //if (vflag_either || eflag_atom)
      if (EVFLAG) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,evdwl,fpair,delx,dely,delz);
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

// Experimental thread-safety using duplicated data instead of atomics

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::operator()(TagPairExp6rxComputeNoAtomics<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii, EV_FLOAT& ev) const {

  {
    const bool one_type = (ntypes == 1);
    if (isite1 == isite2)
      if (one_type)
        this->vectorized_operator<NEIGHFLAG,NEWTON_PAIR,EVFLAG,true, false, true>(ii, ev);
      else
        this->vectorized_operator<NEIGHFLAG,NEWTON_PAIR,EVFLAG,true, false,false>(ii, ev);
    else
      if (one_type)
        this->vectorized_operator<NEIGHFLAG,NEWTON_PAIR,EVFLAG,false,false, true>(ii, ev);
      else
        this->vectorized_operator<NEIGHFLAG,NEWTON_PAIR,EVFLAG,false,false,false>(ii, ev);
    return;
  }

  int tid = 0;
#ifndef LMP_KOKKOS_GPU
  typedef Kokkos::Experimental::UniqueToken<
    DeviceType, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;
  unique_token_type unique_token;
  tid = unique_token.acquire();
#endif

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

      evdwlOldEXP6_12 = 0.0;
      evdwlOldEXP6_21 = 0.0;
      evdwlEXP6_12 = 0.0;
      evdwlEXP6_21 = 0.0;
      fpairOldEXP6_12 = 0.0;
      fpairOldEXP6_21 = 0.0;

      if (rmOld12_ij!=0.0 && rmOld21_ij!=0.0) {
        if (alphaOld21_ij == 6.0 || alphaOld12_ij == 6.0)
          k_error_flag.template view<DeviceType>()() = 1;

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
        if (r < rin1) {
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

        if (r < rin1) {
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
        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
          t_uCG(tid,j) += 0.5*evdwlOld;
      }

      if (rm12_ij!=0.0 && rm21_ij!=0.0) {
        if (alpha21_ij == 6.0 || alpha12_ij == 6.0)
          k_error_flag.template view<DeviceType>()() = 1;

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

        if (r < rin1) {
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

        if (r < rin1) {
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
      }

      //
      // Apply Mixing Rule to get the overall force for the CG pair
      //
      if (isite1 == isite2) fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12;
      else fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12 + sqrt(mixWtSite2old_i*mixWtSite1old_j)*fpairOldEXP6_21;

      fx_i += delx*fpair;
      fy_i += dely*fpair;
      fz_i += delz*fpair;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        t_f(tid,j,0) -= delx*fpair;
        t_f(tid,j,1) -= dely*fpair;
        t_f(tid,j,2) -= delz*fpair;
      }

      if (isite1 == isite2) evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12;
      else evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12 + sqrt(mixWtSite2_i*mixWtSite1_j)*evdwlEXP6_21;
      evdwl *= factor_lj;

      uCGnew_i += 0.5*evdwl;
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal))
        t_uCGnew(tid,j) += 0.5*evdwl;
      evdwl = evdwlOld;
      if (EVFLAG)
        ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*evdwl;
      //if (vflag_either || eflag_atom)
      if (EVFLAG) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,evdwl,fpair,delx,dely,delz);
    }
  }

  t_f(tid,i,0) += fx_i;
  t_f(tid,i,1) += fy_i;
  t_f(tid,i,2) += fz_i;
  t_uCG(tid,i) += uCG_i;
  t_uCGnew(tid,i) += uCGnew_i;

#ifndef LMP_KOKKOS_GPU
  unique_token.release(tid);
#endif
}

// Experimental thread-safe approach using duplicated data instead of atomics and
// temporary local short vector arrays for the inner j-loop to increase vectorization.

template<int n>
  KOKKOS_INLINE_FUNCTION
double __powint(const double& x, const int)
{
   static_assert(n == 12, "__powint<> only supports specific integer powers.");

   if (n == 12)
   {
     // Do x^12 here ... x^12 = (x^3)^4
     double x3 = x*x*x;
     return x3*x3*x3*x3;
   }
}

template<class DeviceType>
  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, bool Site1EqSite2, bool UseAtomics, bool OneType>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::vectorized_operator(const int &ii, EV_FLOAT& ev) const
{
  // These arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_uCG = uCG;
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_uCGnew = uCGnew;

  int tid = 0;
#ifndef LMP_KOKKOS_GPU
  typedef Kokkos::Experimental::UniqueToken<
    DeviceType, Kokkos::Experimental::UniqueTokenScope::Global> unique_token_type;
  unique_token_type unique_token;
  tid = unique_token.acquire();
#endif

  const int nRep = 12;
  const double shift = 1.05;

  const int i = d_ilist[ii];
  const double xtmp = x(i,0);
  const double ytmp = x(i,1);
  const double ztmp = x(i,2);
  const int itype = type[i];
  const int jnum = d_numneigh[i];

  double fx_i = 0.0;
  double fy_i = 0.0;
  double fz_i = 0.0;
  double uCG_i = 0.0;
  double uCGnew_i = 0.0;

  // Constant values for this atom.
  const double epsilon1_i      = PairExp6ParamData.epsilon1[i];
  const double alpha1_i        = PairExp6ParamData.alpha1[i];
  const double rm1_i           = PairExp6ParamData.rm1[i];
  const double mixWtSite1_i    = PairExp6ParamData.mixWtSite1[i];
  const double epsilon2_i      = PairExp6ParamData.epsilon2[i];
  const double alpha2_i        = PairExp6ParamData.alpha2[i];
  const double rm2_i           = PairExp6ParamData.rm2[i];
  const double mixWtSite2_i    = PairExp6ParamData.mixWtSite2[i];
  const double epsilonOld1_i   = PairExp6ParamData.epsilonOld1[i];
  const double alphaOld1_i     = PairExp6ParamData.alphaOld1[i];
  const double rmOld1_i        = PairExp6ParamData.rmOld1[i];
  const double mixWtSite1old_i = PairExp6ParamData.mixWtSite1old[i];
  const double epsilonOld2_i   = PairExp6ParamData.epsilonOld2[i];
  const double alphaOld2_i     = PairExp6ParamData.alphaOld2[i];
  const double rmOld2_i        = PairExp6ParamData.rmOld2[i];
  const double mixWtSite2old_i = PairExp6ParamData.mixWtSite2old[i];

  const double cutsq_type11 = d_cutsq(1,1);
  const double rCut2inv_type11 = 1.0/ cutsq_type11;
  const double rCut6inv_type11 = rCut2inv_type11*rCut2inv_type11*rCut2inv_type11;
  const double rCut_type11 = sqrt( cutsq_type11 );
  const double rCutInv_type11 = 1.0/rCut_type11;

  // Do error testing locally.
  bool hasError = false;

  // Process this many neighbors concurrently -- if possible.
  const int batchSize = 8;

  int neigh_j[batchSize];
  double evdwlOld_j[batchSize];
  double uCGnew_j[batchSize];
  double fpair_j[batchSize];
  double delx_j[batchSize];
  double dely_j[batchSize];
  double delz_j[batchSize];
  double cutsq_j[batchSize];

  for (int jptr = 0; jptr < jnum; )
  {
    // The core computation here is very expensive so let's only bother with
    // those that pass rsq < cutsq.

    for (int j = 0; j < batchSize; ++j)
    {
      evdwlOld_j[j] = 0.0;
      uCGnew_j[j] = 0.0;
      fpair_j[j] = 0.0;
      //delx_j[j] = 0.0;
      //dely_j[j] = 0.0;
      //delz_j[j] = 0.0;
      //cutsq_j[j] = 0.0;
    }

    int niters = 0;

    for (; (jptr < jnum) && (niters < batchSize); ++jptr)
    {
      const int j = d_neighbors(i,jptr) & NEIGHMASK;

      const double delx = xtmp - x(j,0);
      const double dely = ytmp - x(j,1);
      const double delz = ztmp - x(j,2);

      const double rsq = delx*delx + dely*dely + delz*delz;
      const int jtype = type[j];

      const double cutsq_ij = (OneType) ? cutsq_type11 : d_cutsq(itype,jtype);

      if (rsq < cutsq_ij)
      {
        delx_j [niters] = delx;
        dely_j [niters] = dely;
        delz_j [niters] = delz;
        if (!OneType)
          cutsq_j[niters] = cutsq_ij;

        neigh_j[niters] = d_neighbors(i,jptr);

        ++niters;
      }
    }

    // reduction here.
    #ifdef KOKKOS_ENABLE_PRAGMA_SIMD
    #pragma simd reduction(+: fx_i, fy_i, fz_i, uCG_i, uCGnew_i) reduction(|: hasError)
    #endif
    for (int jlane = 0; jlane < niters; jlane++)
    {
      int j = neigh_j[jlane];
      const double factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = delx_j[jlane];
      const double dely = dely_j[jlane];
      const double delz = delz_j[jlane];

      const double rsq = delx*delx + dely*dely + delz*delz;
      // const int jtype = type[j];

      // if (rsq < d_cutsq(itype,jtype)) // optimize
      {
        const double r2inv = 1.0/rsq;
        const double r6inv = r2inv*r2inv*r2inv;

        const double r = sqrt(rsq);
        const double rCut2inv = (OneType) ? rCut2inv_type11 : (1.0/ cutsq_j[jlane]);
        const double rCut6inv = (OneType) ? rCut6inv_type11 : (rCut2inv*rCut2inv*rCut2inv);
        const double rCut =     (OneType) ? rCut_type11     : (sqrt( cutsq_j[jlane] ));
        const double rCutInv =  (OneType) ? rCutInv_type11  : (1.0/rCut);

        //
        // A. Compute the exp-6 potential
        //

        // A1.  Get alpha, epsilon and rm for particle j

        const double epsilon1_j      = PairExp6ParamData.epsilon1[j];
        const double alpha1_j        = PairExp6ParamData.alpha1[j];
        const double rm1_j           = PairExp6ParamData.rm1[j];
        const double mixWtSite1_j    = PairExp6ParamData.mixWtSite1[j];
        const double epsilon2_j      = PairExp6ParamData.epsilon2[j];
        const double alpha2_j        = PairExp6ParamData.alpha2[j];
        const double rm2_j           = PairExp6ParamData.rm2[j];
        const double mixWtSite2_j    = PairExp6ParamData.mixWtSite2[j];
        const double epsilonOld1_j   = PairExp6ParamData.epsilonOld1[j];
        const double alphaOld1_j     = PairExp6ParamData.alphaOld1[j];
        const double rmOld1_j        = PairExp6ParamData.rmOld1[j];
        const double mixWtSite1old_j = PairExp6ParamData.mixWtSite1old[j];
        const double epsilonOld2_j   = PairExp6ParamData.epsilonOld2[j];
        const double alphaOld2_j     = PairExp6ParamData.alphaOld2[j];
        const double rmOld2_j        = PairExp6ParamData.rmOld2[j];
        const double mixWtSite2old_j = PairExp6ParamData.mixWtSite2old[j];

        // A2.  Apply Lorentz-Berthelot mixing rules for the i-j pair
        const double alphaOld12_ij = sqrt(alphaOld1_i*alphaOld2_j);
        const double rmOld12_ij = 0.5*(rmOld1_i + rmOld2_j);
        const double epsilonOld12_ij = sqrt(epsilonOld1_i*epsilonOld2_j);
        const double alphaOld21_ij = sqrt(alphaOld2_i*alphaOld1_j);
        const double rmOld21_ij = 0.5*(rmOld2_i + rmOld1_j);
        const double epsilonOld21_ij = sqrt(epsilonOld2_i*epsilonOld1_j);

        const double alpha12_ij = sqrt(alpha1_i*alpha2_j);
        const double rm12_ij = 0.5*(rm1_i + rm2_j);
        const double epsilon12_ij = sqrt(epsilon1_i*epsilon2_j);
        const double alpha21_ij = sqrt(alpha2_i*alpha1_j);
        const double rm21_ij = 0.5*(rm2_i + rm1_j);
        const double epsilon21_ij = sqrt(epsilon2_i*epsilon1_j);

        double evdwlOldEXP6_12 = 0.0;
        double evdwlOldEXP6_21 = 0.0;
        double evdwlEXP6_12 = 0.0;
        double evdwlEXP6_21 = 0.0;
        double fpairOldEXP6_12 = 0.0;
        double fpairOldEXP6_21 = 0.0;

        if (rmOld12_ij!=0.0 && rmOld21_ij!=0.0)
        {
          hasError |= (alphaOld21_ij == 6.0 || alphaOld12_ij == 6.0);

          // A3.  Compute some convenient quantities for evaluating the force
          double rminv = 1.0/rmOld12_ij;
          double buck1 = epsilonOld12_ij / (alphaOld12_ij - 6.0);
          double rexp = expValue(alphaOld12_ij*(1.0-r*rminv));
          double rm2ij = rmOld12_ij*rmOld12_ij;
          double rm6ij = rm2ij*rm2ij*rm2ij;

          // Compute the shifted potential
          double rCutExp = expValue(alphaOld12_ij*(1.0-rCut*rminv));
          double buck2 = 6.0*alphaOld12_ij;
          double urc = buck1*(6.0*rCutExp - alphaOld12_ij*rm6ij*rCut6inv);
          double durc = -buck1*buck2*(rCutExp* rminv - rCutInv*rm6ij*rCut6inv);
          double rin1 = shift*rmOld12_ij*func_rin(alphaOld12_ij);

          if (r < rin1) {
            const double rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
            const double rin6inv = 1.0/rin6;

            const double rin1exp = expValue(alphaOld12_ij*(1.0-rin1*rminv));

            const double uin1 = buck1*(6.0*rin1exp - alphaOld12_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

            const double win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

            const double aRep = win1*__powint<12>(rin1,nRep)/nRep;

            const double uin1rep = aRep/__powint<12>(rin1,nRep);

            const double forceExp6 = double(nRep)*aRep/__powint<12>(r,nRep);
            fpairOldEXP6_12 = factor_lj*forceExp6*r2inv;

            evdwlOldEXP6_12 = uin1 - uin1rep + aRep/__powint<12>(r,nRep);
          } else {
            const double forceExp6 = buck1*buck2*(r*rexp*rminv - rm6ij*r6inv) + r*durc;
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

          if (r < rin1) {
            const double rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
            const double rin6inv = 1.0/rin6;

            const double rin1exp = expValue(alphaOld21_ij*(1.0-rin1*rminv));

            const double uin1 = buck1*(6.0*rin1exp - alphaOld21_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

            const double win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

            const double aRep = win1*__powint<12>(rin1,nRep)/nRep;

            const double uin1rep = aRep/__powint<12>(rin1,nRep);

            const double forceExp6 = double(nRep)*aRep/__powint<12>(r,nRep);
            fpairOldEXP6_21 = factor_lj*forceExp6*r2inv;

            evdwlOldEXP6_21 = uin1 - uin1rep + aRep/__powint<12>(r,nRep);
          } else {
            const double forceExp6 = buck1*buck2*(r*rexp*rminv - rm6ij*r6inv) + r*durc;
            fpairOldEXP6_21 = factor_lj*forceExp6*r2inv;

            evdwlOldEXP6_21 = buck1*(6.0*rexp - alphaOld21_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
          }

          double evdwlOld;
          if (Site1EqSite2)
            evdwlOld = sqrt(mixWtSite1old_i*mixWtSite2old_j)*evdwlOldEXP6_12;
          else
            evdwlOld = sqrt(mixWtSite1old_i*mixWtSite2old_j)*evdwlOldEXP6_12 + sqrt(mixWtSite2old_i*mixWtSite1old_j)*evdwlOldEXP6_21;

          evdwlOld *= factor_lj;

          uCG_i += 0.5*evdwlOld;

          evdwlOld_j[jlane] = evdwlOld;
        }

        if (rm12_ij!=0.0 && rm21_ij!=0.0)
        {
          hasError |= (alpha21_ij == 6.0 || alpha12_ij == 6.0);

          // A3.  Compute some convenient quantities for evaluating the force
          double rminv = 1.0/rm12_ij;
          double buck1 = epsilon12_ij / (alpha12_ij - 6.0);
          double buck2 = 6.0*alpha12_ij;
          double rexp = expValue(alpha12_ij*(1.0-r*rminv));
          double rm2ij = rm12_ij*rm12_ij;
          double rm6ij = rm2ij*rm2ij*rm2ij;

          // Compute the shifted potential
          double rCutExp = expValue(alpha12_ij*(1.0-rCut*rminv));
          double urc = buck1*(6.0*rCutExp - alpha12_ij*rm6ij*rCut6inv);
          double durc = -buck1*buck2*(rCutExp*rminv - rCutInv*rm6ij*rCut6inv);
          double rin1 = shift*rm12_ij*func_rin(alpha12_ij);

          if (r < rin1) {
            const double rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
            const double rin6inv = 1.0/rin6;

            const double rin1exp = expValue(alpha12_ij*(1.0-rin1*rminv));

            const double uin1 = buck1*(6.0*rin1exp - alpha12_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

            const double win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

            const double aRep = win1*__powint<12>(rin1,nRep)/nRep;

            const double uin1rep = aRep/__powint<12>(rin1,nRep);

            evdwlEXP6_12 = uin1 - uin1rep + aRep/__powint<12>(r,nRep);
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

          if (r < rin1) {
            const double rin6 = rin1*rin1*rin1*rin1*rin1*rin1;
            const double rin6inv = 1.0/rin6;

            const double rin1exp = expValue(alpha21_ij*(1.0-rin1*rminv));

            const double uin1 = buck1*(6.0*rin1exp - alpha21_ij*rm6ij*rin6inv) - urc - durc*(rin1-rCut);

            const double win1 = buck1*buck2*(rin1*rin1exp*rminv - rm6ij*rin6inv) + rin1*durc;

            const double aRep = win1*__powint<12>(rin1,nRep)/nRep;

            const double uin1rep = aRep/__powint<12>(rin1,nRep);

            evdwlEXP6_21 = uin1 - uin1rep + aRep/__powint<12>(r,nRep);
          } else {
            evdwlEXP6_21 = buck1*(6.0*rexp - alpha21_ij*rm6ij*r6inv) - urc - durc*(r-rCut);
          }
        }

        //
        // Apply Mixing Rule to get the overall force for the CG pair
        //
        double fpair;
        if (Site1EqSite2)
          fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12;
        else
          fpair = sqrt(mixWtSite1old_i*mixWtSite2old_j)*fpairOldEXP6_12 + sqrt(mixWtSite2old_i*mixWtSite1old_j)*fpairOldEXP6_21;

        double evdwl;
        if (Site1EqSite2)
          evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12;
        else
          evdwl = sqrt(mixWtSite1_i*mixWtSite2_j)*evdwlEXP6_12 + sqrt(mixWtSite2_i*mixWtSite1_j)*evdwlEXP6_21;

        evdwl *= factor_lj;

        fpair_j[jlane] = fpair;

        fx_i += delx*fpair;
        fy_i += dely*fpair;
        fz_i += delz*fpair;

        uCGnew_i += 0.5*evdwl;
        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD))
          uCGnew_j[jlane] = 0.5*evdwl;

      } // if rsq < cutsq

    } // end jlane loop.

    for (int jlane = 0; jlane < niters; jlane++)
    {
      const int j = neigh_j[jlane] & NEIGHMASK;

      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        if (UseAtomics) {
          a_uCG(j) += 0.5*evdwlOld_j[jlane];
        } else {
          t_uCG(tid,j) += 0.5*evdwlOld_j[jlane];
        }
      }
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        if (UseAtomics) {
          a_uCGnew(j) += uCGnew_j[jlane];
        } else {
          t_uCGnew(tid,j) += uCGnew_j[jlane];
        }
      }
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
        if (UseAtomics) {
          a_f(j,0) -= delx_j[jlane]*fpair_j[jlane];
          a_f(j,1) -= dely_j[jlane]*fpair_j[jlane];
          a_f(j,2) -= delz_j[jlane]*fpair_j[jlane];
        } else {
          t_f(tid,j,0) -= delx_j[jlane]*fpair_j[jlane];
          t_f(tid,j,1) -= dely_j[jlane]*fpair_j[jlane];
          t_f(tid,j,2) -= delz_j[jlane]*fpair_j[jlane];
        }
      }

      double evdwl = evdwlOld_j[jlane];
      if (EVFLAG)
        ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR||(j<nlocal)))?1.0:0.5)*evdwl;
      //if (vflag_either || eflag_atom)
      if (EVFLAG) this->template ev_tally<NEIGHFLAG,NEWTON_PAIR>(ev,i,j,evdwl,fpair_j[jlane],delx_j[jlane],dely_j[jlane],delz_j[jlane]);
    }
  }

  if (hasError)
    k_error_flag.template view<DeviceType>()() = 1;

  if (UseAtomics)
  {
    a_f(i,0) += fx_i;
    a_f(i,1) += fy_i;
    a_f(i,2) += fz_i;
    a_uCG(i) += uCG_i;
    a_uCGnew(i) += uCGnew_i;
  }
  else
  {
    t_f(tid,i,0) += fx_i;
    t_f(tid,i,1) += fy_i;
    t_f(tid,i,2) += fz_i;
    t_uCG(tid,i) += uCG_i;
    t_uCGnew(tid,i) += uCGnew_i;
  }

#ifndef LMP_KOKKOS_GPU
  unique_token.release(tid);
#endif
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::operator()(TagPairExp6rxComputeNoAtomics<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(TagPairExp6rxComputeNoAtomics<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ii, ev);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::operator()(TagPairExp6rxCollapseDupViews, const int &i) const {
  for (int n = 0; n < nthreads; n++) {
    f(i,0) += t_f(n,i,0);
    f(i,1) += t_f(n,i,1);
    f(i,2) += t_f(n,i,2);
    uCG(i) += t_uCG(n,i);
    uCGnew(i) += t_uCGnew(n,i);
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::operator()(TagPairExp6rxZeroDupViews, const int &i) const {
  for (int n = 0; n < nthreads; n++) {
    t_f(n,i,0) = 0.0;
    t_f(n,i,1) = 0.0;
    t_f(n,i,2) = 0.0;
    t_uCG(n,i) = 0.0;
    t_uCGnew(n,i) = 0.0;
  }
}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::allocate()
{
  allocated = 1;
  ntypes = atom->ntypes;

  memory->create(setflag,ntypes+1,ntypes+1,"pair:setflag");
  for (int i = 1; i <= ntypes; i++)
    for (int j = i; j <= ntypes; j++)
      setflag[i][j] = 0;

  memoryKK->create_kokkos(k_cutsq,cutsq,ntypes+1,ntypes+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_cutsq.template modify<LMPHostType>();

  memory->create(cut,ntypes+1,ntypes+1,"pair:cut_lj");
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

  k_params.template modify<LMPHostType>();
  k_params.template sync<DeviceType>();
  d_params = k_params.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::read_file(char *file)
{
  int params_per_line = 5;
  char **words = new char*[params_per_line+1];

  memoryKK->destroy_kokkos(k_params,params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  fp = nullptr;
  if (comm->me == 0) {
    fp = utils::open_potential(file,lmp,nullptr);
    if (fp == nullptr) {
      char str[128];
      snprintf(str,128,"Cannot open exp6/rx potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each set of params from potential file
  // one set of params can span multiple lines

  int n,nwords,ispecies;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (true) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == nullptr) {
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
    nwords = utils::count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == nullptr) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = utils::count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in exp6/rx potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(nullptr," \t\n\r\f"))) continue;

    for (ispecies = 0; ispecies < nspecies; ispecies++)
      if (strcmp(words[0],&atom->dvname[ispecies][0]) == 0) break;
    if (ispecies == nspecies) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      k_params.template modify<LMPHostType>();
      maxparam += DELTA;
      memoryKK->grow_kokkos(k_params,params,maxparam,
                          "pair:params");
    }

    params[nparams].ispecies = ispecies;
    params[nparams].name = utils::strdup(&atom->dvname[ispecies][0]);
    params[nparams].potential = utils::strdup(words[1]);

    if (strcmp(params[nparams].potential,"exp6") == 0) {
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
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairExp6rxKokkos<DeviceType>::setup()
{
  int i,j,n;

  // set mol2param for all combinations
  // must be a single exact match to lines read from file

  memoryKK->destroy_kokkos(k_mol2param,mol2param);
  memoryKK->create_kokkos(k_mol2param,mol2param,nspecies,"pair:mol2param");

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
  for (int ispecies = 0; ispecies < nspecies; ispecies++) {
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
  if (nTotal < MY_EPSILON || nTotalold < MY_EPSILON)
    k_error_flag.template view<DeviceType>()() = 1;

  // Compute the mole fraction of molecules within the fluid portion of the particle (One Fluid Approximation)
  fractionOFAold = nMoleculesOFAold / nTotalold;
  fractionOFA = nMoleculesOFA / nTotal;

  for (int ispecies = 0; ispecies < nspecies; ispecies++) {
    iparam = d_mol2param[ispecies];
    if (iparam < 0 || d_params[iparam].potentialType != exp6PotentialType ) continue;

    // If Site1 matches a pure species, then grab the parameters
    if (isite1 == d_params[iparam].ispecies) {
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
    if (isite2 == d_params[iparam].ispecies) {
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
      fraction2 = nMolecules2/nTotal;
    }

    // If Site1 or Site2 matches is a fluid, then compute the parameters
    if (isOneFluidApprox(isite1) || isOneFluidApprox(isite2)) {
      if (isite1 == d_params[iparam].ispecies || isite2 == d_params[iparam].ispecies) continue;
      rmi = d_params[iparam].rm;
      epsiloni = d_params[iparam].epsilon;
      alphai = d_params[iparam].alpha;
      if (nMoleculesOFA<MY_EPSILON) xMolei = 0.0;
      else xMolei = dvector(ispecies,id)/nMoleculesOFA;
      if (nMoleculesOFAold<MY_EPSILON) xMolei_old = 0.0;
      else xMolei_old = dvector(ispecies+nspecies,id)/nMoleculesOFAold;

      for (int jspecies = 0; jspecies < nspecies; jspecies++) {
        jparam = d_mol2param[jspecies];
        if (jparam < 0 || d_params[jparam].potentialType != exp6PotentialType ) continue;
        if (isite1 == d_params[jparam].ispecies || isite2 == d_params[jparam].ispecies) continue;
        rmj = d_params[jparam].rm;
        epsilonj = d_params[jparam].epsilon;
        alphaj = d_params[jparam].alpha;
        if (nMoleculesOFA<MY_EPSILON) xMolej = 0.0;
        else xMolej = dvector(jspecies,id)/nMoleculesOFA;
        if (nMoleculesOFAold<MY_EPSILON) xMolej_old = 0.0;
        else xMolej_old = dvector(jspecies+nspecies,id)/nMoleculesOFAold;

        rmij = (rmi+rmj)/2.0;
        rm3ij = rmij*rmij*rmij;
        epsilonij = sqrt(epsiloni*epsilonj);
        alphaij = sqrt(alphai*alphaj);

        if (fractionOFAold > 0.0) {
          rm3_old += xMolei_old*xMolej_old*rm3ij;
          epsilon_old += xMolei_old*xMolej_old*rm3ij*epsilonij;
          alpha_old += xMolei_old*xMolej_old*rm3ij*epsilonij*alphaij;
        }
        if (fractionOFA > 0.0) {
          rm3 += xMolei*xMolej*rm3ij;
          epsilon += xMolei*xMolej*rm3ij*epsilonij;
          alpha += xMolei*xMolej*rm3ij*epsilonij*alphaij;
        }
      }
    }
  }

  if (isOneFluidApprox(isite1)) {
    rm1 = cbrt(rm3);
    if (rm1 < MY_EPSILON) {
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
    if (rm1_old < MY_EPSILON) {
      rm1_old = 0.0;
      epsilon1_old = 0.0;
      alpha1_old = 0.0;
    } else {
      epsilon1_old = epsilon_old / rm3_old;
      alpha1_old = alpha_old / epsilon1_old / rm3_old;
    }
    nMoleculesOld1 = 1.0-(nTotalold-nMoleculesOFAold);
    fractionOld1 = fractionOFAold;

    if (scalingFlag == EXPONENT) {
      exponentScaling(nMoleculesOFA,epsilon1,rm1);
      exponentScaling(nMoleculesOFAold,epsilon1_old,rm1_old);
    } else if (scalingFlag == POLYNOMIAL) {
      polynomialScaling(nMoleculesOFA,alpha1,epsilon1,rm1);
      polynomialScaling(nMoleculesOFAold,alpha1_old,epsilon1_old,rm1_old);
    }
  }

  if (isOneFluidApprox(isite2)) {
    rm2 = cbrt(rm3);
    if (rm2 < MY_EPSILON) {
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
    if (rm2_old < MY_EPSILON) {
      rm2_old = 0.0;
      epsilon2_old = 0.0;
      alpha2_old = 0.0;
    } else {
      epsilon2_old = epsilon_old / rm3_old;
      alpha2_old = alpha_old / epsilon2_old / rm3_old;
    }
    nMoleculesOld2 = 1.0-(nTotalold-nMoleculesOFAold);
    fractionOld2 = fractionOFAold;

    if (scalingFlag == EXPONENT) {
      exponentScaling(nMoleculesOFA,epsilon2,rm2);
      exponentScaling(nMoleculesOFAold,epsilon2_old,rm2_old);
    } else if (scalingFlag == POLYNOMIAL) {
      polynomialScaling(nMoleculesOFA,alpha2,epsilon2,rm2);
      polynomialScaling(nMoleculesOFAold,alpha2_old,epsilon2_old,rm2_old);
    }
  }

  // Check that no fractions are less than zero
  if (fraction1 < 0.0 || nMolecules1 < 0.0) {
    if (fraction1 < -MY_EPSILON || nMolecules1 < -MY_EPSILON) {
      k_error_flag.template view<DeviceType>()() = 2;
    }
    nMolecules1 = 0.0;
    fraction1 = 0.0;
  }
  if (fraction2 < 0.0 || nMolecules2 < 0.0) {
    if (fraction2 < -MY_EPSILON || nMolecules2 < -MY_EPSILON) {
      k_error_flag.template view<DeviceType>()() = 2;
    }
    nMolecules2 = 0.0;
    fraction2 = 0.0;
  }
  if (fractionOld1 < 0.0 || nMoleculesOld1 < 0.0) {
    if (fractionOld1 < -MY_EPSILON || nMoleculesOld1 < -MY_EPSILON) {
      k_error_flag.template view<DeviceType>()() = 2;
    }
    nMoleculesOld1 = 0.0;
    fractionOld1 = 0.0;
  }
  if (fractionOld2 < 0.0 || nMoleculesOld2 < 0.0) {
    if (fractionOld2 < -MY_EPSILON || nMoleculesOld2 < -MY_EPSILON) {
      k_error_flag.template view<DeviceType>()() = 2;
    }
    nMoleculesOld2 = 0.0;
    fractionOld2 = 0.0;
  }

  if (fractionalWeighting) {
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

#ifdef _OPENMP
void partition_range( const int begin, const int end, int &thread_begin, int &thread_end, const int chunkSize = 1)
{
   int threadId = omp_get_thread_num();
   int nThreads = omp_get_num_threads();

   const int len = end - begin;
   const int nBlocks = (len + (chunkSize - 1)) / chunkSize;
   const int nBlocksPerThread = nBlocks / nThreads;
   const int nRemaining = nBlocks - nBlocksPerThread * nThreads;
   int block_lo, block_hi;
   if (threadId < nRemaining)
   {
      block_lo = threadId * nBlocksPerThread + threadId;
      block_hi = block_lo + nBlocksPerThread + 1;
   }
   else
   {
      block_lo = threadId * nBlocksPerThread + nRemaining;
      block_hi = block_lo + nBlocksPerThread;
   }

   thread_begin = std::min(begin + block_lo * chunkSize, end);
   thread_end   = std::min(begin + block_hi * chunkSize, end);
   //printf("tid: %d %d %d %d %d\n", threadId, block_lo, block_hi, thread_begin, thread_end);
}
#endif

/* ---------------------------------------------------------------------- */

#ifndef LMP_KOKKOS_GPU
template<class DeviceType>
  template<class ArrayT>
void PairExp6rxKokkos<DeviceType>::getMixingWeightsVect(const int np_total, int errorFlag,
                          ArrayT &epsilon1, ArrayT &alpha1, ArrayT &rm1,  ArrayT &mixWtSite1, ArrayT &epsilon2, ArrayT &alpha2, ArrayT &rm2, ArrayT &mixWtSite2, ArrayT &epsilon1_old, ArrayT &alpha1_old, ArrayT &rm1_old,  ArrayT &mixWtSite1old, ArrayT &epsilon2_old, ArrayT &alpha2_old, ArrayT &rm2_old, ArrayT &mixWtSite2old) const
{
  ArrayT epsilon          = PairExp6ParamDataVect.epsilon         ;
  ArrayT rm3              = PairExp6ParamDataVect.rm3             ;
  ArrayT alpha            = PairExp6ParamDataVect.alpha           ;
  ArrayT xMolei           = PairExp6ParamDataVect.xMolei          ;

  ArrayT epsilon_old      = PairExp6ParamDataVect.epsilon_old     ;
  ArrayT rm3_old          = PairExp6ParamDataVect.rm3_old         ;
  ArrayT alpha_old        = PairExp6ParamDataVect.alpha_old       ;
  ArrayT xMolei_old       = PairExp6ParamDataVect.xMolei_old      ;

  ArrayT fractionOFA      = PairExp6ParamDataVect.fractionOFA     ;
  ArrayT fraction1        = PairExp6ParamDataVect.fraction1       ;
  ArrayT fraction2        = PairExp6ParamDataVect.fraction2       ;
  ArrayT nMoleculesOFA    = PairExp6ParamDataVect.nMoleculesOFA   ;
  ArrayT nMolecules1      = PairExp6ParamDataVect.nMolecules1     ;
  ArrayT nMolecules2      = PairExp6ParamDataVect.nMolecules2     ;
  ArrayT nTotal           = PairExp6ParamDataVect.nTotal          ;

  ArrayT fractionOFAold   = PairExp6ParamDataVect.fractionOFAold  ;
  ArrayT fractionOld1     = PairExp6ParamDataVect.fractionOld1    ;
  ArrayT fractionOld2     = PairExp6ParamDataVect.fractionOld2    ;
  ArrayT nMoleculesOFAold = PairExp6ParamDataVect.nMoleculesOFAold;
  ArrayT nMoleculesOld1   = PairExp6ParamDataVect.nMoleculesOld1  ;
  ArrayT nMoleculesOld2   = PairExp6ParamDataVect.nMoleculesOld2  ;
  ArrayT nTotalold        = PairExp6ParamDataVect.nTotalold       ;

  int errorFlag1 = 0, errorFlag2 = 0;

#ifdef _OPENMP
  #pragma omp parallel reduction(+: errorFlag1, errorFlag2)
#endif
  {
    int idx_begin = 0, idx_end = np_total;
#ifdef _OPENMP
    partition_range( 0, np_total, idx_begin, idx_end, 16 );
#endif

  // Zero out all of the terms first.
  #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
  #pragma ivdep
  #endif
  for (int id = idx_begin; id < idx_end; ++id)
  {
     rm3[id] = 0.0;
     epsilon[id] = 0.0;
     alpha[id] = 0.0;
     epsilon_old[id] = 0.0;
     rm3_old[id] = 0.0;
     alpha_old[id] = 0.0;
     fractionOFA[id] = 0.0;
     fractionOFAold[id] = 0.0;
     nMoleculesOFA[id] = 0.0;
     nMoleculesOFAold[id] = 0.0;
     nTotal[id] = 0.0;
     nTotalold[id] = 0.0;
  }

  // Compute the total number of molecules in the old and new CG particle as well as the total number of molecules in the fluid portion of the old and new CG particle
  for (int ispecies = 0; ispecies < nspecies; ispecies++)
  {
    #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
    #pragma ivdep
    #endif
    for (int id = idx_begin; id < idx_end; ++id)
    {
      nTotal[id] += dvector(ispecies,id);
      nTotalold[id] += dvector(ispecies+nspecies,id);
    }

    const int iparam = d_mol2param[ispecies];

    if (iparam < 0 || d_params[iparam].potentialType != exp6PotentialType ) continue;
    if (isOneFluidApprox(isite1) || isOneFluidApprox(isite2)) {
      if (isite1 == d_params[iparam].ispecies || isite2 == d_params[iparam].ispecies) continue;

      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      for (int id = idx_begin; id < idx_end; ++id)
      {
        nMoleculesOFAold[id] += dvector(ispecies+nspecies,id);
        nMoleculesOFA[id] += dvector(ispecies,id);
      }
    }
  }

  // Make a reduction.
  #ifdef KOKKOS_ENABLE_PRAGMA_simd
  #pragma omp simd reduction(+:errorFlag1)
  #endif
  for (int id = idx_begin; id < idx_end; ++id)
  {
    if (nTotal[id] < MY_EPSILON || nTotalold[id] < MY_EPSILON)
      errorFlag1 = 1;

    // Compute the mole fraction of molecules within the fluid portion of the particle (One Fluid Approximation)
    fractionOFAold[id] = nMoleculesOFAold[id] / nTotalold[id];
    fractionOFA[id] = nMoleculesOFA[id] / nTotal[id];
  }

  for (int ispecies = 0; ispecies < nspecies; ispecies++) {
    const int iparam = d_mol2param[ispecies];
    if (iparam < 0 || d_params[iparam].potentialType != exp6PotentialType ) continue;

    // If Site1 matches a pure species, then grab the parameters
    if (isite1 == d_params[iparam].ispecies)
    {
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      for (int id = idx_begin; id < idx_end; ++id)
      {
        rm1_old[id] = d_params[iparam].rm;
        rm1[id] = d_params[iparam].rm;
        epsilon1_old[id] = d_params[iparam].epsilon;
        epsilon1[id] = d_params[iparam].epsilon;
        alpha1_old[id] = d_params[iparam].alpha;
        alpha1[id] = d_params[iparam].alpha;

        // Compute the mole fraction of Site1
        nMoleculesOld1[id] = dvector(ispecies+nspecies,id);
        nMolecules1[id] = dvector(ispecies,id);
        fractionOld1[id] = nMoleculesOld1[id]/nTotalold[id];
        fraction1[id] = nMolecules1[id]/nTotal[id];
      }
    }

    // If Site2 matches a pure species, then grab the parameters
    if (isite2 == d_params[iparam].ispecies)
    {
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      for (int id = idx_begin; id < idx_end; ++id)
      {
        rm2_old[id] = d_params[iparam].rm;
        rm2[id] = d_params[iparam].rm;
        epsilon2_old[id] = d_params[iparam].epsilon;
        epsilon2[id] = d_params[iparam].epsilon;
        alpha2_old[id] = d_params[iparam].alpha;
        alpha2[id] = d_params[iparam].alpha;

        // Compute the mole fraction of Site2
        nMoleculesOld2[id] = dvector(ispecies+nspecies,id);
        nMolecules2[id] = dvector(ispecies,id);
        fractionOld2[id] = nMoleculesOld2[id]/nTotalold[id];
        fraction2[id] = nMolecules2[id]/nTotal[id];
      }
    }

    // If Site1 or Site2 matches is a fluid, then compute the parameters
    if (isOneFluidApprox(isite1) || isOneFluidApprox(isite2)) {
      if (isite1 == d_params[iparam].ispecies || isite2 == d_params[iparam].ispecies) continue;

      const double rmi = d_params[iparam].rm;
      const double epsiloni = d_params[iparam].epsilon;
      const double alphai = d_params[iparam].alpha;

      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      for (int id = idx_begin; id < idx_end; ++id)
      {
        if (nMoleculesOFA[id]<MY_EPSILON) xMolei[id] = 0.0;
        else xMolei[id] = dvector(ispecies,id)/nMoleculesOFA[id];
        if (nMoleculesOFAold[id]<MY_EPSILON) xMolei_old[id] = 0.0;
        else xMolei_old[id] = dvector(ispecies+nspecies,id)/nMoleculesOFAold[id];
      }

      for (int jspecies = 0; jspecies < nspecies; jspecies++) {
        const int jparam = d_mol2param[jspecies];
        if (jparam < 0 || d_params[jparam].potentialType != exp6PotentialType ) continue;
        if (isite1 == d_params[jparam].ispecies || isite2 == d_params[jparam].ispecies) continue;

        const double rmj = d_params[jparam].rm;
        const double epsilonj = d_params[jparam].epsilon;
        const double alphaj = d_params[jparam].alpha;

        const double rmij = (rmi+rmj)/2.0;
        const double rm3ij = rmij*rmij*rmij;
        const double epsilonij = sqrt(epsiloni*epsilonj);
        const double alphaij = sqrt(alphai*alphaj);

        #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
        #pragma ivdep
        #endif
        for (int id = idx_begin; id < idx_end; ++id)
        {
          double xMolej, xMolej_old;
          if (nMoleculesOFA[id]<MY_EPSILON) xMolej = 0.0;
          else xMolej = dvector(jspecies,id)/nMoleculesOFA[id];
          if (nMoleculesOFAold[id]<MY_EPSILON) xMolej_old = 0.0;
          else xMolej_old = dvector(jspecies+nspecies,id)/nMoleculesOFAold[id];

          if (fractionOFAold[id] > 0.0) {
            rm3_old[id] += xMolei_old[id]*xMolej_old*rm3ij;
            epsilon_old[id] += xMolei_old[id]*xMolej_old*rm3ij*epsilonij;
            alpha_old[id] += xMolei_old[id]*xMolej_old*rm3ij*epsilonij*alphaij;
          }
          if (fractionOFA[id] > 0.0) {
            rm3[id] += xMolei[id]*xMolej*rm3ij;
            epsilon[id] += xMolei[id]*xMolej*rm3ij*epsilonij;
            alpha[id] += xMolei[id]*xMolej*rm3ij*epsilonij*alphaij;
          }
        }
      }
    }
  }

  if (isOneFluidApprox(isite1))
  {
    #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
    #pragma ivdep
    #endif
    for (int id = idx_begin; id < idx_end; ++id)
    {
      rm1[id] = cbrt(rm3[id]);
      if (rm1[id] < MY_EPSILON) {
        rm1[id] = 0.0;
        epsilon1[id] = 0.0;
        alpha1[id] = 0.0;
      } else {
        epsilon1[id] = epsilon[id] / rm3[id];
        alpha1[id] = alpha[id] / epsilon1[id] / rm3[id];
      }
      nMolecules1[id] = 1.0-(nTotal[id]-nMoleculesOFA[id]);
      fraction1[id] = fractionOFA[id];

      rm1_old[id] = cbrt(rm3_old[id]);
      if (rm1_old[id] < MY_EPSILON) {
        rm1_old[id] = 0.0;
        epsilon1_old[id] = 0.0;
        alpha1_old[id] = 0.0;
      } else {
        epsilon1_old[id] = epsilon_old[id] / rm3_old[id];
        alpha1_old[id] = alpha_old[id] / epsilon1_old[id] / rm3_old[id];
      }
      nMoleculesOld1[id] = 1.0-(nTotalold[id]-nMoleculesOFAold[id]);
      fractionOld1[id] = fractionOFAold[id];
    }

    if (scalingFlag == EXPONENT) {
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      for (int id = idx_begin; id < idx_end; ++id)
      {
        exponentScaling(nMoleculesOFA[id],epsilon1[id],rm1[id]);
        exponentScaling(nMoleculesOFAold[id],epsilon1_old[id],rm1_old[id]);
      }
    }
    else if (scalingFlag == POLYNOMIAL) {
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      for (int id = idx_begin; id < idx_end; ++id)
      {
        polynomialScaling(nMoleculesOFA[id],alpha1[id],epsilon1[id],rm1[id]);
        polynomialScaling(nMoleculesOFAold[id],alpha1_old[id],epsilon1_old[id],rm1_old[id]);
      }
    }
  }

  if (isOneFluidApprox(isite2))
  {
    #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
    #pragma ivdep
    #endif
    for (int id = idx_begin; id < idx_end; ++id)
    {
      rm2[id] = cbrt(rm3[id]);
      if (rm2[id] < MY_EPSILON) {
        rm2[id] = 0.0;
        epsilon2[id] = 0.0;
        alpha2[id] = 0.0;
      } else {
        epsilon2[id] = epsilon[id] / rm3[id];
        alpha2[id] = alpha[id] / epsilon2[id] / rm3[id];
      }
      nMolecules2[id] = 1.0-(nTotal[id]-nMoleculesOFA[id]);
      fraction2[id] = fractionOFA[id];

      rm2_old[id] = cbrt(rm3_old[id]);
      if (rm2_old[id] < MY_EPSILON) {
        rm2_old[id] = 0.0;
        epsilon2_old[id] = 0.0;
        alpha2_old[id] = 0.0;
      } else {
        epsilon2_old[id] = epsilon_old[id] / rm3_old[id];
        alpha2_old[id] = alpha_old[id] / epsilon2_old[id] / rm3_old[id];
      }
      nMoleculesOld2[id] = 1.0-(nTotalold[id]-nMoleculesOFAold[id]);
      fractionOld2[id] = fractionOFAold[id];
    }

    if (scalingFlag == EXPONENT) {
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      for (int id = idx_begin; id < idx_end; ++id)
      {
        exponentScaling(nMoleculesOFA[id],epsilon2[id],rm2[id]);
        exponentScaling(nMoleculesOFAold[id],epsilon2_old[id],rm2_old[id]);
      }
    }
    else if (scalingFlag == POLYNOMIAL) {
      #ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      for (int id = idx_begin; id < idx_end; ++id)
      {
        polynomialScaling(nMoleculesOFA[id],alpha2[id],epsilon2[id],rm2[id]);
        polynomialScaling(nMoleculesOFAold[id],alpha2_old[id],epsilon2_old[id],rm2_old[id]);
      }
    }
  }

  // Check that no fractions are less than zero
  #ifdef KOKKOS_ENABLE_PRAGMA_SIMD
  #pragma omp simd reduction(+:errorFlag2)
  #endif
  for (int id = idx_begin; id < idx_end; ++id)
  {
    if (fraction1[id] < 0.0 || nMolecules1[id] < 0.0) {
      if (fraction1[id] < -MY_EPSILON || nMolecules1[id] < -MY_EPSILON) {
        errorFlag2 = 2;
      }
      nMolecules1[id] = 0.0;
      fraction1[id] = 0.0;
    }
    if (fraction2[id] < 0.0 || nMolecules2[id] < 0.0) {
      if (fraction2[id] < -MY_EPSILON || nMolecules2[id] < -MY_EPSILON) {
        errorFlag2 = 2;
      }
      nMolecules2[id] = 0.0;
      fraction2[id] = 0.0;
    }
    if (fractionOld1[id] < 0.0 || nMoleculesOld1[id] < 0.0) {
      if (fractionOld1[id] < -MY_EPSILON || nMoleculesOld1[id] < -MY_EPSILON) {
        errorFlag2 = 2;
      }
      nMoleculesOld1[id] = 0.0;
      fractionOld1[id] = 0.0;
    }
    if (fractionOld2[id] < 0.0 || nMoleculesOld2[id] < 0.0) {
      if (fractionOld2[id] < -MY_EPSILON || nMoleculesOld2[id] < -MY_EPSILON) {
        errorFlag2 = 2;
      }
      nMoleculesOld2[id] = 0.0;
      fractionOld2[id] = 0.0;
    }

    if (fractionalWeighting) {
      mixWtSite1old[id] = fractionOld1[id];
      mixWtSite1[id] = fraction1[id];
      mixWtSite2old[id] = fractionOld2[id];
      mixWtSite2[id] = fraction2[id];
    } else {
      mixWtSite1old[id] = nMoleculesOld1[id];
      mixWtSite1[id] = nMolecules1[id];
      mixWtSite2old[id] = nMoleculesOld2[id];
      mixWtSite2[id] = nMolecules2[id];
    }
  }

  } // end parallel region

  if (errorFlag1 > 0)
    errorFlag = 1;

  if (errorFlag2 > 0)
    errorFlag = 2;
}
#endif

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairExp6rxKokkos<DeviceType>::exponentScaling(double phi, double &epsilon, double &rm) const
{
  double powfuch;

  if (exponentEpsilon < 0.0) {
    powfuch = pow(phi,-exponentEpsilon);
    if (powfuch<MY_EPSILON) epsilon = 0.0;
    else epsilon *= 1.0/powfuch;
  } else {
    epsilon *= pow(phi,exponentEpsilon);
  }

  if (exponentR < 0.0) {
    powfuch = pow(phi,-exponentR);
    if (powfuch<MY_EPSILON) rm = 0.0;
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
    rm *= (s_coeffRm[0]*phi5 + s_coeffRm[1]*phi4 + s_coeffRm[2]*phi3 + s_coeffRm[3]*phi2 + s_coeffRm[4]*phi + s_coeffRm[5]);
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
  if (value < DBL_MIN_EXP) returnValue = 0.0;
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
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();

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
#ifdef LMP_KOKKOS_GPU
template class PairExp6rxKokkos<LMPHostType>;
#endif
}
