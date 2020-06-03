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
   Contributing author: Ray Shan (Materials Design)
------------------------------------------------------------------------- */

#include "angle_class2_kokkos.h"
#include <cmath>
#include <cstdlib>
#include "atom_kokkos.h"
#include "neighbor_kokkos.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
AngleClass2Kokkos<Space>::AngleClass2Kokkos(LAMMPS *lmp) : AngleClass2(lmp)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = Space;
  neighborKK = (NeighborKokkos *) neighbor;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
AngleClass2Kokkos<Space>::~AngleClass2Kokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void AngleClass2Kokkos<Space>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"angle:eatom");
    d_eatom = DualViewHelper<Space>::view(k_eatom);
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"angle:vatom");
    d_vatom = DualViewHelper<Space>::view(k_vatom);
  }

  //atomKK->sync(execution_space,datamask_read);
  //if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  //else atomKK->modified(execution_space,F_MASK);

  DualViewHelper<Space>::sync(k_theta0);
  DualViewHelper<Space>::sync(k_k2);
  DualViewHelper<Space>::sync(k_k3);
  DualViewHelper<Space>::sync(k_k4);
  DualViewHelper<Space>::sync(k_bb_k);
  DualViewHelper<Space>::sync(k_bb_r1);
  DualViewHelper<Space>::sync(k_bb_r2);
  DualViewHelper<Space>::sync(k_ba_k1);
  DualViewHelper<Space>::sync(k_ba_k2);
  DualViewHelper<Space>::sync(k_ba_r1);
  DualViewHelper<Space>::sync(k_ba_r2);
  DualViewHelper<Space>::sync(k_setflag);
  DualViewHelper<Space>::sync(k_setflag_a);
  DualViewHelper<Space>::sync(k_setflag_bb);
  DualViewHelper<Space>::sync(k_setflag_ba);

  x = DualViewHelper<Space>::view(atomKK->k_x);
  f = DualViewHelper<Space>::view(atomKK->k_f);
  DualViewHelper<Space>::sync(neighborKK->k_anglelist);
  anglelist = DualViewHelper<Space>::view(neighborKK->k_anglelist);
  int nanglelist = neighborKK->nanglelist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleClass2Compute<1,1> >(0,nanglelist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleClass2Compute<0,1> >(0,nanglelist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleClass2Compute<1,0> >(0,nanglelist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleClass2Compute<0,0> >(0,nanglelist),*this);
    }
  }

  if (eflag_global) energy += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    DualViewHelper<Space>::modify(k_eatom);
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    DualViewHelper<Space>::modify(k_vatom);
    k_vatom.sync_host();
  }

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void AngleClass2Kokkos<Space>::operator()(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<typename AT::t_float_1d_3::data_type, typename AT::t_float_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = anglelist(n,0);
  const int i2 = anglelist(n,1);
  const int i3 = anglelist(n,2);
  const int type = anglelist(n,3);

  // 1st bond

  const KK_FLOAT delx1 = x(i1,0) - x(i2,0);
  const KK_FLOAT dely1 = x(i1,1) - x(i2,1);
  const KK_FLOAT delz1 = x(i1,2) - x(i2,2);

  const KK_FLOAT rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
  const KK_FLOAT r1 = sqrt(rsq1);

  // 2nd bond

  const KK_FLOAT delx2 = x(i3,0) - x(i2,0);
  const KK_FLOAT dely2 = x(i3,1) - x(i2,1);
  const KK_FLOAT delz2 = x(i3,2) - x(i2,2);

  const KK_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
  const KK_FLOAT r2 = sqrt(rsq2);

  // angle (cos and sin)

  KK_FLOAT c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  KK_FLOAT s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;
  s = 1.0/s;

  // force & energy for angle term

  const KK_FLOAT dtheta = acos(c) - d_theta0[type];
  const KK_FLOAT dtheta2 = dtheta*dtheta;
  const KK_FLOAT dtheta3 = dtheta2*dtheta;
  const KK_FLOAT dtheta4 = dtheta3*dtheta;

  const KK_FLOAT de_angle = 2.0*d_k2[type]*dtheta + 3.0*d_k3[type]*dtheta2 + 4.0*d_k4[type]*dtheta3;

  const KK_FLOAT a = -de_angle*s;
  const KK_FLOAT a11 = a*c / rsq1;
  const KK_FLOAT a12 = -a / (r1*r2);
  const KK_FLOAT a22 = a*c / rsq2;

  KK_FLOAT f1[3],f3[3];
  f1[0] = a11*delx1 + a12*delx2;
  f1[1] = a11*dely1 + a12*dely2;
  f1[2] = a11*delz1 + a12*delz2;
  f3[0] = a22*delx2 + a12*delx1;
  f3[1] = a22*dely2 + a12*dely1;
  f3[2] = a22*delz2 + a12*delz1;

  KK_FLOAT eangle = 0.0;
  if (eflag) eangle = d_k2[type]*dtheta2 + d_k3[type]*dtheta3 + d_k4[type]*dtheta4;

  // force & energy for bond-bond term

  const KK_FLOAT dr1 = r1 - d_bb_r1[type];
  const KK_FLOAT dr2 = r2 - d_bb_r2[type];
  const KK_FLOAT tk1 = d_bb_k[type] * dr1;
  const KK_FLOAT tk2 = d_bb_k[type] * dr2;

  f1[0] -= delx1*tk2/r1;
  f1[1] -= dely1*tk2/r1;
  f1[2] -= delz1*tk2/r1;

  f3[0] -= delx2*tk1/r2;
  f3[1] -= dely2*tk1/r2;
  f3[2] -= delz2*tk1/r2;

  if (eflag) eangle += d_bb_k[type]*dr1*dr2;

  // force & energy for bond-angle term

  const KK_FLOAT aa1 = s * dr1 * d_ba_k1[type];
  const KK_FLOAT aa2 = s * dr2 * d_ba_k2[type];

  KK_FLOAT aa11 = aa1 * c / rsq1;
  KK_FLOAT aa12 = -aa1 / (r1 * r2);
  KK_FLOAT aa21 = aa2 * c / rsq1;
  KK_FLOAT aa22 = -aa2 / (r1 * r2);

  const KK_FLOAT vx11 = (aa11 * delx1) + (aa12 * delx2);
  const KK_FLOAT vx12 = (aa21 * delx1) + (aa22 * delx2);
  const KK_FLOAT vy11 = (aa11 * dely1) + (aa12 * dely2);
  const KK_FLOAT vy12 = (aa21 * dely1) + (aa22 * dely2);
  const KK_FLOAT vz11 = (aa11 * delz1) + (aa12 * delz2);
  const KK_FLOAT vz12 = (aa21 * delz1) + (aa22 * delz2);

  aa11 = aa1 * c / rsq2;
  aa21 = aa2 * c / rsq2;

  const KK_FLOAT vx21 = (aa11 * delx2) + (aa12 * delx1);
  const KK_FLOAT vx22 = (aa21 * delx2) + (aa22 * delx1);
  const KK_FLOAT vy21 = (aa11 * dely2) + (aa12 * dely1);
  const KK_FLOAT vy22 = (aa21 * dely2) + (aa22 * dely1);
  const KK_FLOAT vz21 = (aa11 * delz2) + (aa12 * delz1);
  const KK_FLOAT vz22 = (aa21 * delz2) + (aa22 * delz1);

  const KK_FLOAT b1 = d_ba_k1[type] * dtheta / r1;
  const KK_FLOAT b2 = d_ba_k2[type] * dtheta / r2;

  f1[0] -= vx11 + b1*delx1 + vx12;
  f1[1] -= vy11 + b1*dely1 + vy12;
  f1[2] -= vz11 + b1*delz1 + vz12;

  f3[0] -= vx21 + b2*delx2 + vx22;
  f3[1] -= vy21 + b2*dely2 + vy22;
  f3[2] -= vz21 + b2*delz2 + vz22;

  if (eflag) eangle += d_ba_k1[type]*dr1*dtheta + d_ba_k2[type]*dr2*dtheta;

  // apply force to each of 3 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    a_f(i1,0) += f1[0];
    a_f(i1,1) += f1[1];
    a_f(i1,2) += f1[2];
  }

  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) -= f1[0] + f3[0];
    a_f(i2,1) -= f1[1] + f3[1];
    a_f(i2,2) -= f1[2] + f3[2];
  }

  if (NEWTON_BOND || i3 < nlocal) {
    a_f(i3,0) += f3[0];
    a_f(i3,1) += f3[1];
    a_f(i3,2) += f3[2];
  }

  if (EVFLAG) ev_tally(ev,i1,i2,i3,eangle,f1,f3,
                       delx1,dely1,delz1,delx2,dely2,delz2);
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void AngleClass2Kokkos<Space>::operator()(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<ExecutionSpace Space>
void AngleClass2Kokkos<Space>::allocate()
{
  AngleClass2::allocate();

}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void AngleClass2Kokkos<Space>::coeff(int narg, char **arg)
{
  AngleClass2::coeff(narg, arg);

  int n = atom->nangletypes;
  k_k2 = DAT::tdual_float_1d("AngleClass2::k2",n+1);
  k_k3 = DAT::tdual_float_1d("AngleClass2::k3",n+1);
  k_k4 = DAT::tdual_float_1d("AngleClass2::k4",n+1);
  k_bb_k = DAT::tdual_float_1d("AngleClass2::bb_k",n+1);
  k_bb_r1 = DAT::tdual_float_1d("AngleClass2::bb_r1",n+1);
  k_bb_r2 = DAT::tdual_float_1d("AngleClass2::bb_r2",n+1);
  k_ba_k1 = DAT::tdual_float_1d("AngleClass2::ba_k1",n+1);
  k_ba_k2 = DAT::tdual_float_1d("AngleClass2::ba_k2",n+1);
  k_ba_r1 = DAT::tdual_float_1d("AngleClass2::ba_r1",n+1);
  k_ba_r2 = DAT::tdual_float_1d("AngleClass2::ba_r2",n+1);
  k_setflag = DAT::tdual_float_1d("AngleClass2::setflag",n+1);
  k_setflag_a = DAT::tdual_float_1d("AngleClass2::setflag_a",n+1);
  k_setflag_bb = DAT::tdual_float_1d("AngleClass2::setflag_bb",n+1);
  k_setflag_ba = DAT::tdual_float_1d("AngleClass2::setflag_ba",n+1);
  k_theta0 = DAT::tdual_float_1d("AngleClass2::theta0",n+1);

  d_k2 = DualViewHelper<Space>::view(k_k2);
  d_k3 = DualViewHelper<Space>::view(k_k3);
  d_k4 = DualViewHelper<Space>::view(k_k4);
  d_bb_k = DualViewHelper<Space>::view(k_bb_k);
  d_bb_r1 = DualViewHelper<Space>::view(k_bb_r1);
  d_bb_r2 = DualViewHelper<Space>::view(k_bb_r2);
  d_ba_k1 = DualViewHelper<Space>::view(k_ba_k1);
  d_ba_k2 = DualViewHelper<Space>::view(k_ba_k2);
  d_ba_r1 = DualViewHelper<Space>::view(k_ba_r1);
  d_ba_r2 = DualViewHelper<Space>::view(k_ba_r2);
  d_ba_r2 = DualViewHelper<Space>::view(k_ba_r2);
  d_setflag = DualViewHelper<Space>::view(k_setflag);
  d_setflag_a = DualViewHelper<Space>::view(k_setflag_a);
  d_setflag_bb = DualViewHelper<Space>::view(k_setflag_bb);
  d_setflag_ba = DualViewHelper<Space>::view(k_setflag_ba);
  d_theta0 = DualViewHelper<Space>::view(k_theta0);

  //int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_k2.h_view[i] = k2[i];
    k_k3.h_view[i] = k3[i];
    k_k4.h_view[i] = k4[i];
    k_bb_k.h_view[i] = bb_k[i];
    k_bb_r1.h_view[i] = bb_r1[i];
    k_bb_r2.h_view[i] = bb_r2[i];
    k_ba_k1.h_view[i] = ba_k1[i];
    k_ba_k2.h_view[i] = ba_k2[i];
    k_ba_r1.h_view[i] = ba_r1[i];
    k_ba_r2.h_view[i] = ba_r2[i];
    k_setflag.h_view[i] = setflag[i];
    k_setflag_a.h_view[i] = setflag_a[i];
    k_setflag_bb.h_view[i] = setflag_bb[i];
    k_setflag_ba.h_view[i] = setflag_ba[i];
    k_theta0.h_view[i] = theta0[i];
  }

  k_k2.modify_host();
  k_k3.modify_host();
  k_k4.modify_host();
  k_bb_k.modify_host();
  k_bb_r1.modify_host();
  k_bb_r2.modify_host();
  k_ba_k1.modify_host();
  k_ba_k2.modify_host();
  k_ba_r1.modify_host();
  k_ba_r2.modify_host();
  k_setflag.modify_host();
  k_setflag_a.modify_host();
  k_setflag_bb.modify_host();
  k_setflag_ba.modify_host();
  k_theta0.modify_host();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
void AngleClass2Kokkos<Space>::read_restart(FILE *fp)
{
  AngleClass2::read_restart(fp);

  int n = atom->nangletypes;
  k_k2 = DAT::tdual_float_1d("AngleClass2::k2",n+1);
  k_k3 = DAT::tdual_float_1d("AngleClass2::k3",n+1);
  k_k4 = DAT::tdual_float_1d("AngleClass2::k4",n+1);
  k_bb_k = DAT::tdual_float_1d("AngleClass2::bb_k",n+1);
  k_bb_r1 = DAT::tdual_float_1d("AngleClass2::bb_r1",n+1);
  k_bb_r2 = DAT::tdual_float_1d("AngleClass2::bb_r2",n+1);
  k_ba_k1 = DAT::tdual_float_1d("AngleClass2::ba_k1",n+1);
  k_ba_k2 = DAT::tdual_float_1d("AngleClass2::ba_k2",n+1);
  k_ba_r1 = DAT::tdual_float_1d("AngleClass2::ba_r1",n+1);
  k_ba_r2 = DAT::tdual_float_1d("AngleClass2::ba_r2",n+1);
  k_setflag = DAT::tdual_float_1d("AngleClass2::setflag",n+1);
  k_setflag_a = DAT::tdual_float_1d("AngleClass2::setflag_a",n+1);
  k_setflag_bb = DAT::tdual_float_1d("AngleClass2::setflag_bb",n+1);
  k_setflag_ba = DAT::tdual_float_1d("AngleClass2::setflag_ba",n+1);
  k_theta0 = DAT::tdual_float_1d("AngleClass2::theta0",n+1);

  d_k2 = DualViewHelper<Space>::view(k_k2);
  d_k3 = DualViewHelper<Space>::view(k_k3);
  d_k4 = DualViewHelper<Space>::view(k_k4);
  d_bb_k = DualViewHelper<Space>::view(k_bb_k);
  d_bb_r1 = DualViewHelper<Space>::view(k_bb_r1);
  d_bb_r2 = DualViewHelper<Space>::view(k_bb_r2);
  d_ba_k1 = DualViewHelper<Space>::view(k_ba_k1);
  d_ba_k2 = DualViewHelper<Space>::view(k_ba_k2);
  d_ba_r1 = DualViewHelper<Space>::view(k_ba_r1);
  d_ba_r2 = DualViewHelper<Space>::view(k_ba_r2);
  d_ba_r2 = DualViewHelper<Space>::view(k_ba_r2);
  d_setflag = DualViewHelper<Space>::view(k_setflag);
  d_setflag_a = DualViewHelper<Space>::view(k_setflag_a);
  d_setflag_bb = DualViewHelper<Space>::view(k_setflag_bb);
  d_setflag_ba = DualViewHelper<Space>::view(k_setflag_ba);
  d_theta0 = DualViewHelper<Space>::view(k_theta0);

  //int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_k2.h_view[i] = k2[i];
    k_k3.h_view[i] = k3[i];
    k_k4.h_view[i] = k4[i];
    k_bb_k.h_view[i] = bb_k[i];
    k_bb_r1.h_view[i] = bb_r1[i];
    k_bb_r2.h_view[i] = bb_r2[i];
    k_ba_k1.h_view[i] = ba_k1[i];
    k_ba_k2.h_view[i] = ba_k2[i];
    k_ba_r1.h_view[i] = ba_r1[i];
    k_ba_r2.h_view[i] = ba_r2[i];
    k_setflag.h_view[i] = setflag[i];
    k_setflag_a.h_view[i] = setflag_a[i];
    k_setflag_bb.h_view[i] = setflag_bb[i];
    k_setflag_ba.h_view[i] = setflag_ba[i];
    k_theta0.h_view[i] = theta0[i];
  }

  k_k2.modify_host();
  k_k3.modify_host();
  k_k4.modify_host();
  k_bb_k.modify_host();
  k_bb_r1.modify_host();
  k_bb_r2.modify_host();
  k_ba_k1.modify_host();
  k_ba_k2.modify_host();
  k_ba_r1.modify_host();
  k_ba_r2.modify_host();
  k_setflag.modify_host();
  k_setflag_a.modify_host();
  k_setflag_bb.modify_host();
  k_setflag_ba.modify_host();
  k_theta0.modify_host();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

template<ExecutionSpace Space>
//template<int NEWTON_BOND>
KOKKOS_INLINE_FUNCTION
void AngleClass2Kokkos<Space>::ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
                     const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
                     const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const
{
  KK_FLOAT eanglethird;
  KK_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<typename AT::t_float_1d::data_type, typename AT::t_float_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = DualViewHelper<Space>::view(k_eatom);
  Kokkos::View<typename AT::t_float_1d_6::data_type, typename AT::t_float_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = DualViewHelper<Space>::view(k_vatom);

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += eangle;
      else {
        eanglethird = THIRD*eangle;

        if (i < nlocal) ev.evdwl += eanglethird;
        if (j < nlocal) ev.evdwl += eanglethird;
        if (k < nlocal) ev.evdwl += eanglethird;
      }
    }
    if (eflag_atom) {
      eanglethird = THIRD*eangle;

      if (newton_bond || i < nlocal) v_eatom[i] += eanglethird;
      if (newton_bond || j < nlocal) v_eatom[j] += eanglethird;
      if (newton_bond || k < nlocal) v_eatom[k] += eanglethird;
    }
  }

  if (vflag_either) {
    v[0] = delx1*f1[0] + delx2*f3[0];
    v[1] = dely1*f1[1] + dely2*f3[1];
    v[2] = delz1*f1[2] + delz2*f3[2];
    v[3] = delx1*f1[1] + delx2*f3[1];
    v[4] = delx1*f1[2] + delx2*f3[2];
    v[5] = dely1*f1[2] + dely2*f3[2];

    if (vflag_global) {
      if (newton_bond) {
        ev.v[0] += v[0];
        ev.v[1] += v[1];
        ev.v[2] += v[2];
        ev.v[3] += v[3];
        ev.v[4] += v[4];
        ev.v[5] += v[5];
      } else {
        if (i < nlocal) {
          ev.v[0] += THIRD*v[0];
          ev.v[1] += THIRD*v[1];
          ev.v[2] += THIRD*v[2];
          ev.v[3] += THIRD*v[3];
          ev.v[4] += THIRD*v[4];
          ev.v[5] += THIRD*v[5];
        }
        if (j < nlocal) {
          ev.v[0] += THIRD*v[0];
          ev.v[1] += THIRD*v[1];
          ev.v[2] += THIRD*v[2];
          ev.v[3] += THIRD*v[3];
          ev.v[4] += THIRD*v[4];
          ev.v[5] += THIRD*v[5];
        }
        if (k < nlocal) {
          ev.v[0] += THIRD*v[0];

          ev.v[1] += THIRD*v[1];
          ev.v[2] += THIRD*v[2];
          ev.v[3] += THIRD*v[3];
          ev.v[4] += THIRD*v[4];
          ev.v[5] += THIRD*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        v_vatom(i,0) += THIRD*v[0];
        v_vatom(i,1) += THIRD*v[1];
        v_vatom(i,2) += THIRD*v[2];
        v_vatom(i,3) += THIRD*v[3];
        v_vatom(i,4) += THIRD*v[4];
        v_vatom(i,5) += THIRD*v[5];
      }
      if (newton_bond || j < nlocal) {
        v_vatom(j,0) += THIRD*v[0];
        v_vatom(j,1) += THIRD*v[1];
        v_vatom(j,2) += THIRD*v[2];
        v_vatom(j,3) += THIRD*v[3];
        v_vatom(j,4) += THIRD*v[4];
        v_vatom(j,5) += THIRD*v[5];
      }
      if (newton_bond || k < nlocal) {
        v_vatom(k,0) += THIRD*v[0];
        v_vatom(k,1) += THIRD*v[1];
        v_vatom(k,2) += THIRD*v[2];
        v_vatom(k,3) += THIRD*v[3];
        v_vatom(k,4) += THIRD*v[4];
        v_vatom(k,5) += THIRD*v[5];

      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class AngleClass2Kokkos<Device>;
template class AngleClass2Kokkos<Host>;
}

