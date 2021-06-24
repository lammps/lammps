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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "angle_charmm_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleCharmmKokkos<DeviceType>::AngleCharmmKokkos(LAMMPS *lmp) : AngleCharmm(lmp)
{
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleCharmmKokkos<DeviceType>::~AngleCharmmKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleCharmmKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    //if(k_eatom.extent(0)<maxeatom) { // won't work without adding zero functor
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"improper:eatom");
      d_eatom = k_eatom.template view<KKDeviceType>();
    //}
  }
  if (vflag_atom) {
    //if(k_vatom.extent(0)<maxvatom) { // won't work without adding zero functor
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"improper:vatom");
      d_vatom = k_vatom.template view<KKDeviceType>();
    //}
  }

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  neighborKK->k_anglelist.template sync<DeviceType>();
  anglelist = neighborKK->k_anglelist.view<KKDeviceType>();
  int nanglelist = neighborKK->nanglelist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleCharmmCompute<1,1> >(0,nanglelist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleCharmmCompute<0,1> >(0,nanglelist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleCharmmCompute<1,0> >(0,nanglelist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleCharmmCompute<0,0> >(0,nanglelist),*this);
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

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void AngleCharmmKokkos<DeviceType>::operator()(TagAngleCharmmCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

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

  // Urey-Bradley bond

  const F_FLOAT delxUB = x(i3,0) - x(i1,0);
  const F_FLOAT delyUB = x(i3,1) - x(i1,1);
  const F_FLOAT delzUB = x(i3,2) - x(i1,2);

  const F_FLOAT rsqUB = delxUB*delxUB + delyUB*delyUB + delzUB*delzUB;
  const F_FLOAT rUB = sqrt(rsqUB);

  // Urey-Bradley force & energy

  const F_FLOAT dr = rUB - d_r_ub[type];
  const F_FLOAT rk = d_k_ub[type] * dr;

  F_FLOAT forceUB = 0.0;
  if (rUB > 0.0) forceUB = -2.0*rk/rUB;

  E_FLOAT eangle = 0.0;
  if (eflag) eangle = rk*dr;

  // angle (cos and sin)

  F_FLOAT c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  F_FLOAT s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;
  s = 1.0/s;

  // harmonic force & energy

  const F_FLOAT dtheta = acos(c) - d_theta0[type];
  const F_FLOAT tk = d_k[type] * dtheta;

  if (eflag) eangle += tk*dtheta;

  const F_FLOAT a = -2.0 * tk * s;
  const F_FLOAT a11 = a*c / rsq1;
  const F_FLOAT a12 = -a / (r1*r2);
  const F_FLOAT a22 = a*c / rsq2;

  F_FLOAT f1[3],f3[3];
  f1[0] = a11*delx1 + a12*delx2 - delxUB*forceUB;
  f1[1] = a11*dely1 + a12*dely2 - delyUB*forceUB;
  f1[2] = a11*delz1 + a12*delz2 - delzUB*forceUB;

  f3[0] = a22*delx2 + a12*delx1 + delxUB*forceUB;
  f3[1] = a22*dely2 + a12*dely1 + delyUB*forceUB;
  f3[2] = a22*delz2 + a12*delz1 + delzUB*forceUB;

  // apply force to each of 3 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    f(i1,0) += f1[0];
    f(i1,1) += f1[1];
    f(i1,2) += f1[2];
  }

  if (NEWTON_BOND || i2 < nlocal) {
    f(i2,0) -= f1[0] + f3[0];
    f(i2,1) -= f1[1] + f3[1];
    f(i2,2) -= f1[2] + f3[2];
  }

  if (NEWTON_BOND || i3 < nlocal) {
    f(i3,0) += f3[0];
    f(i3,1) += f3[1];
    f(i3,2) += f3[2];
  }

  if (EVFLAG) ev_tally(ev,i1,i2,i3,eangle,f1,f3,
                       delx1,dely1,delz1,delx2,dely2,delz2);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void AngleCharmmKokkos<DeviceType>::operator()(TagAngleCharmmCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleCharmmCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleCharmmKokkos<DeviceType>::allocate()
{
  AngleCharmm::allocate();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleCharmmKokkos<DeviceType>::coeff(int narg, char **arg)
{
  AngleCharmm::coeff(narg, arg);

  int n = atom->nangletypes;
  typename AT::tdual_ffloat_1d k_k("AngleCharmm::k",n+1);
  typename AT::tdual_ffloat_1d k_theta0("AngleCharmm::theta0",n+1);
  typename AT::tdual_ffloat_1d k_k_ub("AngleCharmm::k_ub",n+1);
  typename AT::tdual_ffloat_1d k_r_ub("AngleCharmm::r_ub",n+1);

  d_k = k_k.template view<DeviceType>();
  d_theta0 = k_theta0.template view<DeviceType>();
  d_k_ub = k_k_ub.template view<DeviceType>();
  d_r_ub = k_r_ub.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_theta0.h_view[i] = theta0[i];
    k_k_ub.h_view[i] = k_ub[i];
    k_r_ub.h_view[i] = r_ub[i];
  }

  k_k.template modify<LMPHostType>();
  k_theta0.template modify<LMPHostType>();
  k_k_ub.template modify<LMPHostType>();
  k_r_ub.template modify<LMPHostType>();

  k_k.template sync<DeviceType>();
  k_theta0.template sync<DeviceType>();
  k_k_ub.template sync<DeviceType>();
  k_r_ub.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleCharmmKokkos<DeviceType>::read_restart(FILE *fp)
{
  AngleCharmm::read_restart(fp);

  int n = atom->nangletypes;
  typename AT::tdual_ffloat_1d k_k("AngleCharmm::k",n+1);
  typename AT::tdual_ffloat_1d k_theta0("AngleCharmm::theta0",n+1);
  typename AT::tdual_ffloat_1d k_k_ub("AngleCharmm::k_ub",n+1);
  typename AT::tdual_ffloat_1d k_r_ub("AngleCharmm::r_ub",n+1);

  d_k = k_k.template view<DeviceType>();
  d_theta0 = k_theta0.template view<DeviceType>();
  d_k_ub = k_k_ub.template view<DeviceType>();
  d_r_ub = k_r_ub.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_theta0.h_view[i] = theta0[i];
    k_k_ub.h_view[i] = k_ub[i];
    k_r_ub.h_view[i] = r_ub[i];
  }

  k_k.template modify<LMPHostType>();
  k_theta0.template modify<LMPHostType>();
  k_k_ub.template modify<LMPHostType>();
  k_r_ub.template modify<LMPHostType>();

  k_k.template sync<DeviceType>();
  k_theta0.template sync<DeviceType>();
  k_k_ub.template sync<DeviceType>();
  k_r_ub.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
KOKKOS_INLINE_FUNCTION
void AngleCharmmKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     F_FLOAT &eangle, F_FLOAT *f1, F_FLOAT *f3,
                     const F_FLOAT &delx1, const F_FLOAT &dely1, const F_FLOAT &delz1,
                     const F_FLOAT &delx2, const F_FLOAT &dely2, const F_FLOAT &delz2) const
{
  E_FLOAT eanglethird;
  F_FLOAT v[6];

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

      if (newton_bond || i < nlocal) d_eatom[i] += eanglethird;
      if (newton_bond || j < nlocal) d_eatom[j] += eanglethird;
      if (newton_bond || k < nlocal) d_eatom[k] += eanglethird;
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
        d_vatom(i,0) += THIRD*v[0];
        d_vatom(i,1) += THIRD*v[1];
        d_vatom(i,2) += THIRD*v[2];
        d_vatom(i,3) += THIRD*v[3];
        d_vatom(i,4) += THIRD*v[4];
        d_vatom(i,5) += THIRD*v[5];
      }
      if (newton_bond || j < nlocal) {
        d_vatom(j,0) += THIRD*v[0];
        d_vatom(j,1) += THIRD*v[1];
        d_vatom(j,2) += THIRD*v[2];
        d_vatom(j,3) += THIRD*v[3];
        d_vatom(j,4) += THIRD*v[4];
        d_vatom(j,5) += THIRD*v[5];
      }
      if (newton_bond || k < nlocal) {
        d_vatom(k,0) += THIRD*v[0];
        d_vatom(k,1) += THIRD*v[1];
        d_vatom(k,2) += THIRD*v[2];
        d_vatom(k,3) += THIRD*v[3];
        d_vatom(k,4) += THIRD*v[4];
        d_vatom(k,5) += THIRD*v[5];

      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class AngleCharmmKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class AngleCharmmKokkos<LMPHostType>;
#endif
}

