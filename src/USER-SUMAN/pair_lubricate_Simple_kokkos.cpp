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

#include "pair_lubricate_Simple_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain_kokkos.h"
#include "error.h"
#include "fix.h"
#include "fix_wall.h"
#include "force.h"
#include "input.h"
#include "kokkos.h"
#include "math_const.h"
#include "math_special_kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecialKokkos;

// same as fix_wall.cpp

enum { EDGE, CONSTANT, VARIABLE };

//#define _NO_RANDOM

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLubricateSimpleKokkos<DeviceType>::PairLubricateSimpleKokkos(LAMMPS *lmp) : PairLubricateSimple(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  domainKK = (DomainKokkos *) domain;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | V_MASK | F_MASK | OMEGA_MASK | TORQUE_MASK | TYPE_MASK | VIRIAL_MASK | RADIUS_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLubricateSimpleKokkos<DeviceType>::~PairLubricateSimpleKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cut_inner,cut_inner);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLubricateSimpleKokkos<DeviceType>::init_style()
{
  printf("Inside PairLubricateSimpleKokkos::init_style()\n");
  PairLubricateSimple::init_style();

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
  //  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  //if (neighflag == FULL) request->enable_full();
  // if (neighflag != FULL)
  //   error->all(FLERR,"Must use full neighbor list style with lubricate/simple/kk");
  printf("Leaving PairLubricateSimpleKokkos::init_style()  shearing= %i\n",shearing);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLubricateSimpleKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{ 
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);
  
  // reallocate per-atom arrays if necessary

  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_cut_inner.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  vxmu2f = force->vxmu2f;

  // loop over neighbors of my atoms

  int inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  copymode = 1;

  // FLD contribution to force  & torque

  R0  = 6.0*MY_PI*mu;
  RT0 = 8.0*MY_PI*mu;
  RS0 = 20.0/3.0*MY_PI*mu;

  h_rate = domainKK->h_rate;
  h_ratelo = domainKK->h_ratelo;

  xprd = domainKK->xprd;
  yprd = domainKK->yprd;
  zprd = domainKK->zprd;
  
  if(flagfld) {

    if(shearing) {
      double *h_rate = domainKK->h_rate;
      Ef[0][0] = h_rate[0]/xprd;
      Ef[1][1] = h_rate[1]/yprd;
      Ef[2][2] = h_rate[2]/zprd;
      Ef[0][1] = Ef[1][0] = 0.5 * h_rate[5]/yprd;
      Ef[0][2] = Ef[2][0] = 0.5 * h_rate[4]/zprd;
      Ef[1][2] = Ef[2][1] = 0.5 * h_rate[3]/zprd;
      // needs to be a k_Ef
    }
    
    EV_FLOAT ev;

    // printf("pair::compute()  flagfld= %i  shearing= %i  nlocal= %i  inum= %i  vxmu2f= %f  prd= %f %f %f\n",flagfld,shearing,atom->nlocal,inum,vxmu2f,xprd,yprd,zprd);
    // printf("                 h_rate= %f %f %f %f %f %f\n",h_rate[0],h_rate[1],h_rate[2],h_rate[3],h_rate[4],h_rate[5]);
    // printf("                 h_ratelo= %f %f %f\n",h_ratelo[0],h_ratelo[1],h_ratelo[2]);
    // printf("                 Ef= %f %f %f %f %f %f\n",Ef[0][0],Ef[1][1],Ef[2][2],Ef[0][1],Ef[0][2],Ef[1][2]);
    
#if 1
    domainKK->x2lamda(atom->nlocal);
    if(shearing && vflag_either) {
      if (neighflag == HALF) {
	if(newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleComputeFLD<HALF, 1, 1> >(0,nlocal),*this,ev);
	else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleComputeFLD<HALF, 0, 1> >(0,nlocal),*this,ev);
      } else if (neighflag == HALFTHREAD) {
	if(newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleComputeFLD<HALFTHREAD, 1, 1> >(0,nlocal),*this,ev);
	else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleComputeFLD<HALFTHREAD, 0, 1> >(0,nlocal),*this,ev);
      } else if (neighflag == FULL) {
	if(newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleComputeFLD<FULL, 1, 1> >(0,nlocal),*this,ev);
	else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleComputeFLD<FULL, 0, 1> >(0,nlocal),*this,ev);
      }      
    } else { // value of NEWTON_PAIR not used, so just 0
      if (neighflag == HALF) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleComputeFLD<HALF, 0, 0> >(0,nlocal),*this);
      else if (neighflag == HALFTHREAD) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleComputeFLD<HALFTHREAD, 0, 0> >(0,nlocal),*this);
      else if (neighflag == FULL) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleComputeFLD<FULL, 0, 0> >(0,nlocal),*this);
    }
    domainKK->lamda2x(atom->nlocal);

    if (vflag_global) {
      virial[0] += ev.v[0];
      virial[1] += ev.v[1];
      virial[2] += ev.v[2];
      virial[3] += ev.v[3];
      virial[4] += ev.v[4];
      virial[5] += ev.v[5];
    }
#endif
  }

  EV_FLOAT ev;

  if(flagHI) {
  
    if (flagfld) { // FLAGFLD == 1
      if (vflag_either) { // VFLAG == 1
	if (neighflag == HALF) {
	  if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALF,1,1,1> >(0,inum),*this,ev);
	  else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALF,0,1,1> >(0,inum),*this,ev);
	} else if (neighflag == HALFTHREAD) {
	  if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALFTHREAD,1,1,1> >(0,inum),*this,ev);
	  else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALFTHREAD,0,1,1> >(0,inum),*this,ev);
	} else if (neighflag == FULL) {
	  if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<FULL,1,1,1> >(0,inum),*this,ev);
	  else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<FULL,0,1,1> >(0,inum),*this,ev);
	}
      } else {  // VFLAG==0
	if (neighflag == HALF) {
	  if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALF,1,0,1> >(0,inum),*this);
	  else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALF,0,0,1> >(0,inum),*this);
	} else if (neighflag == HALFTHREAD) {
	  if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALFTHREAD,1,0,1> >(0,inum),*this);
	  else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALFTHREAD,0,0,1> >(0,inum),*this);
	} else if (neighflag == FULL) {
	  if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<FULL,1,0,1> >(0,inum),*this);
	  else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<FULL,0,0,1> >(0,inum),*this);
	}
      }
    } else { // FLAGFLD == 0
      if (evflag) { // VFLAG== 1
	if (neighflag == HALF) {
	  if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALF,1,1,0> >(0,inum),*this,ev);
	  else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALF,0,1,0> >(0,inum),*this,ev);
	} else if (neighflag == HALFTHREAD) {
	  if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALFTHREAD,1,1,0> >(0,inum),*this,ev);
	  else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALFTHREAD,0,1,0> >(0,inum),*this,ev);
	} else if (neighflag == FULL) {
	  if (newton_pair) Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<FULL,1,1,0> >(0,inum),*this,ev);
	  else Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<FULL,0,1,0> >(0,inum),*this,ev);
	}
      } else { // VFLAG == 0
	if (neighflag == HALF) {
	  if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALF,1,0,0> >(0,inum),*this);
	  else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALF,0,0,0> >(0,inum),*this);
	} else if (neighflag == HALFTHREAD) {
	  if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALFTHREAD,1,0,0> >(0,inum),*this);
	  else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<HALFTHREAD,0,0,0> >(0,inum),*this);
	} else if (neighflag == FULL) {
	  if (newton_pair) Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<FULL,1,0,0> >(0,inum),*this);
	  else Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairLubricateSimpleCompute<FULL,0,0,0> >(0,inum),*this);
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

  } // if(flagHI)
  
  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }
  
  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}


template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int SHEARING>
KOKKOS_INLINE_FUNCTION
void PairLubricateSimpleKokkos<DeviceType>::operator()(TagPairLubricateSimpleComputeFLD<NEIGHFLAG, NEWTON_PAIR, SHEARING>, const int ii, EV_FLOAT &ev) const {
  
  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;
  
  const int i = d_ilist[ii];
  const LMP_FLOAT radi = radius[i];

  LMP_FLOAT vstream[3];
  vstream[0] = h_rate[0]*x(i,0) + h_rate[5]*x(i,1) + h_rate[4]*x(i,2) + h_ratelo[0];
  vstream[1] = h_rate[1]*x(i,1) + h_rate[3]*x(i,2) + h_ratelo[1];
  vstream[2] = h_rate[2]*x(i,2) + h_ratelo[2];
  
  a_f(i,0) += vxmu2f*R0*radi*(vstream[0]-v(i,0));
  a_f(i,1) += vxmu2f*R0*radi*(vstream[1]-v(i,1));
  a_f(i,2) += vxmu2f*R0*radi*(vstream[2]-v(i,2));
  
  const LMP_FLOAT radi3 = radi*radi*radi;
  
  a_torque(i,0) -= vxmu2f*RT0*radi3*(omega(i,0)+0.5*h_rate[3]/zprd);
  a_torque(i,1) -= vxmu2f*RT0*radi3*(omega(i,1)-0.5*h_rate[4]/zprd);
  a_torque(i,2) -= vxmu2f*RT0*radi3*(omega(i,2)+0.5*h_rate[5]/yprd);
  
  // Ef = (grad(vstream) + (grad(vstream))^T) / 2
  // set Ef from h_rate in strain units

  if(SHEARING){
    double vRS0 = 1.0; //-vxmu2f * RS0*radi3;
    v_tally_tensor<NEIGHFLAG, NEWTON_PAIR>(ev,i,i,
  					   vRS0*Ef[0][0],vRS0*Ef[1][1],vRS0*Ef[2][2],
  					   vRS0*Ef[0][1],vRS0*Ef[0][2],vRS0*Ef[1][2]);
  }

}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int SHEARING>
KOKKOS_INLINE_FUNCTION
void PairLubricateSimpleKokkos<DeviceType>::operator()(TagPairLubricateSimpleComputeFLD<NEIGHFLAG,NEWTON_PAIR,SHEARING>, const int ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,SHEARING>(TagPairLubricateSimpleComputeFLD<NEIGHFLAG,NEWTON_PAIR,SHEARING>(), ii, ev);
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
KOKKOS_INLINE_FUNCTION
void PairLubricateSimpleKokkos<DeviceType>::operator()(TagPairLubricateSimpleCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>, const int ii, EV_FLOAT &ev) const {
  
  // The f and torque arrays are atomic for Half/Thread neighbor style
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_f = f;
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > a_torque = torque;
  
  const int i = d_ilist[ii];
  const X_FLOAT xtmp = x(i,0);
  const X_FLOAT ytmp = x(i,1);
  const X_FLOAT ztmp = x(i,2);
  const int itype = type[i];
  const LMP_FLOAT radi = radius[i];
  const int jnum = d_numneigh[i];

  LMP_FLOAT wsx,wsy,wsz;
  LMP_FLOAT wdx,wdy,wdz;
  
  F_FLOAT fx_i = 0.0;
  F_FLOAT fy_i = 0.0;
  F_FLOAT fz_i = 0.0;
  
  F_FLOAT torquex_i = 0.0;
  F_FLOAT torquey_i = 0.0;
  F_FLOAT torquez_i = 0.0;
    
  for (int jj = 0; jj < jnum; jj++) {
    int j = d_neighbors(i,jj);
    j &= NEIGHMASK;
    
    const X_FLOAT delx = xtmp - x(j,0);
    const X_FLOAT dely = ytmp - x(j,1);
    const X_FLOAT delz = ztmp - x(j,2);
    const X_FLOAT rsq = delx*delx + dely*dely + delz*delz;
    const int jtype = type[j];
    const LMP_FLOAT radj = radius[j];

    double radsum = radi + radj; // Sum of two particle's radii
      
    if(rsq < d_cutsq(itype,jtype) && rsq < 1.45*1.45*radsum*radsum) {
      const LMP_FLOAT r = sqrt(rsq);
      
      // scalar resistances XA and YA
	
      LMP_FLOAT h_sep = r - radsum;
	
      // if less than minimum gap (scaled inner cutoff by particle radii), use minimum gap instead
	
      if (r < d_cut_inner(itype,jtype)*radsum) h_sep = d_cut_inner(itype,jtype)*radsum - radsum;

      const LMP_FLOAT wsx = (omega(i,0) + omega(j,0)) * 0.5;
      const LMP_FLOAT wsy = (omega(i,1) + omega(j,1)) * 0.5;
      const LMP_FLOAT wsz = (omega(i,2) + omega(j,2)) * 0.5;

      const LMP_FLOAT wdx = (omega(i,0) - omega(j,0)) * 0.5;
      const LMP_FLOAT wdy = (omega(i,1) - omega(j,1)) * 0.5;
      const LMP_FLOAT wdz = (omega(i,2) - omega(j,2)) * 0.5;

      const LMP_FLOAT nx = -delx / r;
      const LMP_FLOAT ny = -dely / r;
      const LMP_FLOAT nz = -delz / r;

      const LMP_FLOAT hinv = (radi + radj) / (2.0 * h_sep);

      const LMP_FLOAT lhinv = log(hinv);
      //if(lhinv < 0.0) error->all(FLERR,"Using pair lubricate with cutoff problem: log(1/h) is negative");

      const LMP_FLOAT beta0 = radj / radi;
      const LMP_FLOAT beta1 = 1.0 + beta0;

      const LMP_FLOAT b0p2 = beta0 * beta0;
      const LMP_FLOAT b1p2 = beta1 * beta1;
      const LMP_FLOAT b1p3 = beta1 * b1p2;
      const LMP_FLOAT mupradi = mu * MY_PI * radi;

      // scalar resistances

      LMP_FLOAT XA11, YA11, YB11, YC12;
      if(flaglog) {
	  XA11 = 6.0*mupradi*( hinv*2.0*b0p2 + lhinv*beta0*(1.0+ 7.0*beta0+ b0p2 )/5.0 )/b1p3;
	  YA11 = 1.6*mupradi*lhinv*beta0*(2.0+beta0+2.0*b0p2)/b1p3;
	  YB11 = -0.8*mupradi*radi*lhinv*beta0*(4.0+beta0)/b1p2;
	  YC12 = 0.8*mupradi*radi*radi*lhinv*b0p2/beta1*(1.0-4.0/beta0); //YC12*(1-4/beta0)
      } else {
	XA11 = 12.0*mupradi*hinv*b0p2/b1p3;
      }

      // relative velocity components U^2-U^1

      LMP_FLOAT vr1 = v(j,0) - v(i,0);
      LMP_FLOAT vr2 = v(j,1) - v(i,1);
      LMP_FLOAT vr3 = v(j,2) - v(i,2);

      // normal component (vr.n)n

      LMP_FLOAT vnnr = vr1*nx + vr2*ny + vr3*nz;
      LMP_FLOAT vn1 = vnnr * nx;
      LMP_FLOAT vn2 = vnnr * ny;
      LMP_FLOAT vn3 = vnnr * nz;

      // tangential component vr - (vr.n)n

      LMP_FLOAT vt1 = vr1 - vn1;
      LMP_FLOAT vt2 = vr2 - vn2;
      LMP_FLOAT vt3 = vr3 - vn3;

      // force due to squeeze type motion : f = XA11 nn dot vr
      
      F_FLOAT fx = XA11 * vn1;
      F_FLOAT fy = XA11 * vn2;
      F_FLOAT fz = XA11 * vn3;

      // force due to all shear kind of motions

      if(flaglog) {
	LMP_FLOAT ybfac = (1.0-beta0*(1.0+4.0*beta0)/(4.0+beta0))*YB11;
	//f+=YA11*vt1-(r1+r2)*YA11*ws cross n + ybfac* wd cross n
	fx += YA11 * (vt1-(radi+radj)*(wsy*nz-wsz*ny)) + ybfac * (wdy*nz-wdz*ny);
	fy += YA11 * (vt2-(radi+radj)*(wsz*nx-wsx*nz)) + ybfac * (wdz*nx-wdx*nz);
	fz += YA11 * (vt3-(radi+radj)*(wsx*ny-wsy*nx)) + ybfac * (wdx*ny-wdy*nx);
      }
      
      // scale forces for appropriate units

      fx *= vxmu2f;
      fy *= vxmu2f;
      fz *= vxmu2f;
	
      fx_i += fx;
      fy_i += fy;
      fz_i += fz;

      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
	a_f(j,0) -= fx;
	a_f(j,1) -= fy;
	a_f(j,2) -= fz;
      }

      // torque due to this force

      if (flaglog) {

	LMP_FLOAT wsdotn = wsx*nx + wsy*ny + wsz*nz;
	LMP_FLOAT wddotn = wdx*nx + wdy*ny + wdz*nz;

	// squeeze+shear+pump contributions to torque
	
	//YB11*(vr cross n + (r1+r2)* tang(ws)) +YC12*(1-4/beta)*tang(wd)

	// YB11 & YC12 are not correct if radi != radj, but pair lubricate requires mono-disperse, so ok...
	
	F_FLOAT tx = YB11*(vr2*nz -vr3*ny + (radi+radj)*(wsx-wsdotn*nx));
	F_FLOAT ty = YB11*(vr3*nx -vr1*nz + (radi+radj)*(wsy-wsdotn*ny));
	F_FLOAT tz = YB11*(vr1*ny -vr2*nx + (radi+radj)*(wsz-wsdotn*nz));

	F_FLOAT ttx = YC12*(wdx-wddotn*nx);
	F_FLOAT tty = YC12*(wdy-wddotn*ny);
	F_FLOAT ttz = YC12*(wdz-wddotn*nz);

	// if(i == 1) printf("i= %i  x= %f %f %f  x= %f %f %f  rsq= %f  tq= %f %f %f\n",i,
	// 		  x(i,0),x(i,1),x(i,2),
	// 		  x(j,0),x(j,1),x(j,2), rsq, tx,ty,tz);
	// if(j == 1) printf("  j= %i  x= %f %f %f  x= %f %f %f  rsq= %f tq= %f %f %f\n",i,
	// 		  x(i,0),x(i,1),x(i,2),
	// 		  x(j,0),x(j,1),x(j,2), rsq, tx,ty,tz);
	
	// torque is same on both particles ?
	
	torquex_i += vxmu2f * (tx + ttx);
	torquey_i += vxmu2f * (ty + tty);
	torquez_i += vxmu2f * (tz + ttz);
	  
	if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || j < nlocal)) {
	  a_torque(j,0) += vxmu2f * (tx - ttx);
	  a_torque(j,1) += vxmu2f * (ty - tty);
	  a_torque(j,2) += vxmu2f * (tz - ttz);
	}
      } // if(flaglog)
      
      if (VFLAG) ev_tally_xyz<NEIGHFLAG, NEWTON_PAIR>(ev, i, j, fx, fy, fz, delx, dely, delz);
      
    } // if(rsq)
  } // for(jj)
  
  a_f(i,0) += fx_i;
  a_f(i,1) += fy_i;
  a_f(i,2) += fz_i;
  a_torque(i,0) += torquex_i;
  a_torque(i,1) += torquey_i;
  a_torque(i,2) += torquez_i;
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int VFLAG, int FLAGFLD>
KOKKOS_INLINE_FUNCTION
void PairLubricateSimpleKokkos<DeviceType>::operator()(TagPairLubricateSimpleCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>, const int ii) const {
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>(TagPairLubricateSimpleCompute<NEIGHFLAG,NEWTON_PAIR,VFLAG,FLAGFLD>(), ii, ev);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairLubricateSimpleKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT & ev, int i, int j,
                                                        F_FLOAT fx, F_FLOAT fy, F_FLOAT fz,
                                                        X_FLOAT delx, X_FLOAT dely, X_FLOAT delz) const
{
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();

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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairLubricateSimpleKokkos<DeviceType>::v_tally_tensor(EV_FLOAT & ev, int i, int j,
							   F_FLOAT vxx, F_FLOAT vyy, F_FLOAT vzz,
							   F_FLOAT vxy, F_FLOAT vxz, F_FLOAT vyz) const
{
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value> > v_vatom = k_vatom.view<DeviceType>();

  const F_FLOAT v0 = vxx;
  const F_FLOAT v1 = vyy;
  const F_FLOAT v2 = vzz;
  const F_FLOAT v3 = vxy;
  const F_FLOAT v4 = vxz;
  const F_FLOAT v5 = vyz;
  
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

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLubricateSimpleKokkos<DeviceType>::allocate()
{
  printf("Inside PairLubricateSimpleKokkos::allocate()\n");
  
  PairLubricateSimple::allocate();

  int n = atom->ntypes;
  
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  
  memory->destroy(cut_inner);
  memoryKK->create_kokkos(k_cut_inner,cut_inner,n+1,n+1,"pair:cut_inner");
  d_cut_inner = k_cut_inner.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLubricateSimpleKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg != 5 && narg != 7) error->all(FLERR, "Illegal pair_style command");

  PairLubricateSimple::settings(narg,arg);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLubricateSimpleKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLubricateSimple::init_one(i,j);
  double cutinnerm = cut_inner[i][j];

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  
  k_cut_inner.h_view(i,j) = k_cut_inner.h_view(j,i) = cutinnerm;
  k_cut_inner.template modify<LMPHostType>();

  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLubricateSimpleKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairLubricateSimple::coeff(narg,arg);
}

namespace LAMMPS_NS {
template class PairLubricateSimpleKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLubricateSimpleKokkos<LMPHostType>;
#endif
}
