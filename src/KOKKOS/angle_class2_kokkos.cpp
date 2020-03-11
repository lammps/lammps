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

template<class DeviceType>
AngleClass2Kokkos<DeviceType>::AngleClass2Kokkos(LAMMPS *lmp) : AngleClass2(lmp)
{
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleClass2Kokkos<DeviceType>::~AngleClass2Kokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleClass2Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"angle:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"angle:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  //atomKK->sync(execution_space,datamask_read);
  //if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  //else atomKK->modified(execution_space,F_MASK);

  k_theta0.template sync<DeviceType>();
  k_k2.template sync<DeviceType>();
  k_k3.template sync<DeviceType>();
  k_k4.template sync<DeviceType>();
  k_bb_k.template sync<DeviceType>();
  k_bb_r1.template sync<DeviceType>();
  k_bb_r2.template sync<DeviceType>();
  k_ba_k1.template sync<DeviceType>();
  k_ba_k2.template sync<DeviceType>();
  k_ba_r1.template sync<DeviceType>();
  k_ba_r2.template sync<DeviceType>();
  k_setflag.template sync<DeviceType>();
  k_setflag_a.template sync<DeviceType>();
  k_setflag_bb.template sync<DeviceType>();
  k_setflag_ba.template sync<DeviceType>();

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_anglelist.template sync<DeviceType>();
  anglelist = neighborKK->k_anglelist.template view<DeviceType>();
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
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void AngleClass2Kokkos<DeviceType>::operator()(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = anglelist(n,0);
  const int i2 = anglelist(n,1);
  const int i3 = anglelist(n,2);
  const int type = anglelist(n,3);

  // 1st bond

  const F_FLOAT delx1 = x(i1,0) - x(i2,0);
  const F_FLOAT dely1 = x(i1,1) - x(i2,1);
  const F_FLOAT delz1 = x(i1,2) - x(i2,2);

  const F_FLOAT rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
  const F_FLOAT r1 = sqrt(rsq1);

  // 2nd bond

  const F_FLOAT delx2 = x(i3,0) - x(i2,0);
  const F_FLOAT dely2 = x(i3,1) - x(i2,1);
  const F_FLOAT delz2 = x(i3,2) - x(i2,2);

  const F_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
  const F_FLOAT r2 = sqrt(rsq2);

  // angle (cos and sin)

  F_FLOAT c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  F_FLOAT s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;
  s = 1.0/s;

  // force & energy for angle term

  const F_FLOAT dtheta = acos(c) - d_theta0[type];
  const F_FLOAT dtheta2 = dtheta*dtheta;
  const F_FLOAT dtheta3 = dtheta2*dtheta;
  const F_FLOAT dtheta4 = dtheta3*dtheta;

  const F_FLOAT de_angle = 2.0*d_k2[type]*dtheta + 3.0*d_k3[type]*dtheta2 + 4.0*d_k4[type]*dtheta3;

  const F_FLOAT a = -de_angle*s;
  const F_FLOAT a11 = a*c / rsq1;
  const F_FLOAT a12 = -a / (r1*r2);
  const F_FLOAT a22 = a*c / rsq2;

  F_FLOAT f1[3],f3[3];
  f1[0] = a11*delx1 + a12*delx2;
  f1[1] = a11*dely1 + a12*dely2;
  f1[2] = a11*delz1 + a12*delz2;
  f3[0] = a22*delx2 + a12*delx1;
  f3[1] = a22*dely2 + a12*dely1;
  f3[2] = a22*delz2 + a12*delz1;

  F_FLOAT eangle = 0.0;
  if (eflag) eangle = d_k2[type]*dtheta2 + d_k3[type]*dtheta3 + d_k4[type]*dtheta4;

  // force & energy for bond-bond term

  const F_FLOAT dr1 = r1 - d_bb_r1[type];
  const F_FLOAT dr2 = r2 - d_bb_r2[type];
  const F_FLOAT tk1 = d_bb_k[type] * dr1;
  const F_FLOAT tk2 = d_bb_k[type] * dr2;

  f1[0] -= delx1*tk2/r1;
  f1[1] -= dely1*tk2/r1;
  f1[2] -= delz1*tk2/r1;

  f3[0] -= delx2*tk1/r2;
  f3[1] -= dely2*tk1/r2;
  f3[2] -= delz2*tk1/r2;

  if (eflag) eangle += d_bb_k[type]*dr1*dr2;

  // force & energy for bond-angle term

  const F_FLOAT aa1 = s * dr1 * d_ba_k1[type];
  const F_FLOAT aa2 = s * dr2 * d_ba_k2[type];

  F_FLOAT aa11 = aa1 * c / rsq1;
  F_FLOAT aa12 = -aa1 / (r1 * r2);
  F_FLOAT aa21 = aa2 * c / rsq1;
  F_FLOAT aa22 = -aa2 / (r1 * r2);

  const F_FLOAT vx11 = (aa11 * delx1) + (aa12 * delx2);
  const F_FLOAT vx12 = (aa21 * delx1) + (aa22 * delx2);
  const F_FLOAT vy11 = (aa11 * dely1) + (aa12 * dely2);
  const F_FLOAT vy12 = (aa21 * dely1) + (aa22 * dely2);
  const F_FLOAT vz11 = (aa11 * delz1) + (aa12 * delz2);
  const F_FLOAT vz12 = (aa21 * delz1) + (aa22 * delz2);

  aa11 = aa1 * c / rsq2;
  aa21 = aa2 * c / rsq2;

  const F_FLOAT vx21 = (aa11 * delx2) + (aa12 * delx1);
  const F_FLOAT vx22 = (aa21 * delx2) + (aa22 * delx1);
  const F_FLOAT vy21 = (aa11 * dely2) + (aa12 * dely1);
  const F_FLOAT vy22 = (aa21 * dely2) + (aa22 * dely1);
  const F_FLOAT vz21 = (aa11 * delz2) + (aa12 * delz1);
  const F_FLOAT vz22 = (aa21 * delz2) + (aa22 * delz1);

  const F_FLOAT b1 = d_ba_k1[type] * dtheta / r1;
  const F_FLOAT b2 = d_ba_k2[type] * dtheta / r2;

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

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void AngleClass2Kokkos<DeviceType>::operator()(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleClass2Kokkos<DeviceType>::allocate()
{
  AngleClass2::allocate();

}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleClass2Kokkos<DeviceType>::coeff(int narg, char **arg)
{
  AngleClass2::coeff(narg, arg);

  int n = atom->nangletypes;
  k_k2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::k2",n+1);
  k_k3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::k3",n+1);
  k_k4 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::k4",n+1);
  k_bb_k = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::bb_k",n+1);
  k_bb_r1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::bb_r1",n+1);
  k_bb_r2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::bb_r2",n+1);
  k_ba_k1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::ba_k1",n+1);
  k_ba_k2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::ba_k2",n+1);
  k_ba_r1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::ba_r1",n+1);
  k_ba_r2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::ba_r2",n+1);
  k_setflag = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag",n+1);
  k_setflag_a = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_a",n+1);
  k_setflag_bb = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_bb",n+1);
  k_setflag_ba = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_ba",n+1);
  k_theta0 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::theta0",n+1);

  d_k2 = k_k2.template view<DeviceType>();
  d_k3 = k_k3.template view<DeviceType>();
  d_k4 = k_k4.template view<DeviceType>();
  d_bb_k = k_bb_k.template view<DeviceType>();
  d_bb_r1 = k_bb_r1.template view<DeviceType>();
  d_bb_r2 = k_bb_r2.template view<DeviceType>();
  d_ba_k1 = k_ba_k1.template view<DeviceType>();
  d_ba_k2 = k_ba_k2.template view<DeviceType>();
  d_ba_r1 = k_ba_r1.template view<DeviceType>();
  d_ba_r2 = k_ba_r2.template view<DeviceType>();
  d_ba_r2 = k_ba_r2.template view<DeviceType>();
  d_setflag = k_setflag.template view<DeviceType>();
  d_setflag_a = k_setflag_a.template view<DeviceType>();
  d_setflag_bb = k_setflag_bb.template view<DeviceType>();
  d_setflag_ba = k_setflag_ba.template view<DeviceType>();
  d_theta0 = k_theta0.template view<DeviceType>();

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

  k_k2.template modify<LMPHostType>();
  k_k3.template modify<LMPHostType>();
  k_k4.template modify<LMPHostType>();
  k_bb_k.template modify<LMPHostType>();
  k_bb_r1.template modify<LMPHostType>();
  k_bb_r2.template modify<LMPHostType>();
  k_ba_k1.template modify<LMPHostType>();
  k_ba_k2.template modify<LMPHostType>();
  k_ba_r1.template modify<LMPHostType>();
  k_ba_r2.template modify<LMPHostType>();
  k_setflag.template modify<LMPHostType>();
  k_setflag_a.template modify<LMPHostType>();
  k_setflag_bb.template modify<LMPHostType>();
  k_setflag_ba.template modify<LMPHostType>();
  k_theta0.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleClass2Kokkos<DeviceType>::read_restart(FILE *fp)
{
  AngleClass2::read_restart(fp);

  int n = atom->nangletypes;
  k_k2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::k2",n+1);
  k_k3 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::k3",n+1);
  k_k4 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::k4",n+1);
  k_bb_k = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::bb_k",n+1);
  k_bb_r1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::bb_r1",n+1);
  k_bb_r2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::bb_r2",n+1);
  k_ba_k1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::ba_k1",n+1);
  k_ba_k2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::ba_k2",n+1);
  k_ba_r1 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::ba_r1",n+1);
  k_ba_r2 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::ba_r2",n+1);
  k_setflag = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag",n+1);
  k_setflag_a = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_a",n+1);
  k_setflag_bb = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_bb",n+1);
  k_setflag_ba = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::setflag_ba",n+1);
  k_theta0 = typename ArrayTypes<DeviceType>::tdual_ffloat_1d("AngleClass2::theta0",n+1);

  d_k2 = k_k2.template view<DeviceType>();
  d_k3 = k_k3.template view<DeviceType>();
  d_k4 = k_k4.template view<DeviceType>();
  d_bb_k = k_bb_k.template view<DeviceType>();
  d_bb_r1 = k_bb_r1.template view<DeviceType>();
  d_bb_r2 = k_bb_r2.template view<DeviceType>();
  d_ba_k1 = k_ba_k1.template view<DeviceType>();
  d_ba_k2 = k_ba_k2.template view<DeviceType>();
  d_ba_r1 = k_ba_r1.template view<DeviceType>();
  d_ba_r2 = k_ba_r2.template view<DeviceType>();
  d_ba_r2 = k_ba_r2.template view<DeviceType>();
  d_setflag = k_setflag.template view<DeviceType>();
  d_setflag_a = k_setflag_a.template view<DeviceType>();
  d_setflag_bb = k_setflag_bb.template view<DeviceType>();
  d_setflag_ba = k_setflag_ba.template view<DeviceType>();
  d_theta0 = k_theta0.template view<DeviceType>();

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

  k_k2.template modify<LMPHostType>();
  k_k3.template modify<LMPHostType>();
  k_k4.template modify<LMPHostType>();
  k_bb_k.template modify<LMPHostType>();
  k_bb_r1.template modify<LMPHostType>();
  k_bb_r2.template modify<LMPHostType>();
  k_ba_k1.template modify<LMPHostType>();
  k_ba_k2.template modify<LMPHostType>();
  k_ba_r1.template modify<LMPHostType>();
  k_ba_r2.template modify<LMPHostType>();
  k_setflag.template modify<LMPHostType>();
  k_setflag_a.template modify<LMPHostType>();
  k_setflag_bb.template modify<LMPHostType>();
  k_setflag_ba.template modify<LMPHostType>();
  k_theta0.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
KOKKOS_INLINE_FUNCTION
void AngleClass2Kokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     F_FLOAT &eangle, F_FLOAT *f1, F_FLOAT *f3,
                     const F_FLOAT &delx1, const F_FLOAT &dely1, const F_FLOAT &delz1,
                     const F_FLOAT &delx2, const F_FLOAT &dely2, const F_FLOAT &delz2) const
{
  E_FLOAT eanglethird;
  F_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = k_eatom.template view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = k_vatom.template view<DeviceType>();

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
template class AngleClass2Kokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class AngleClass2Kokkos<LMPHostType>;
#endif
}

