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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "improper_harmonic_kokkos.h"
#include <cmath>
#include "atom_kokkos.h"
#include "comm.h"
#include "neighbor_kokkos.h"
#include "force.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

static constexpr double TOLERANCE = 0.05;
static constexpr double SMALL =     0.001;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ImproperHarmonicKokkos<DeviceType>::ImproperHarmonicKokkos(LAMMPS *lmp) : ImproperHarmonic(lmp)
{
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_warning_flag = Kokkos::DualView<int,DeviceType>("Dihedral:warning_flag");
  d_warning_flag = k_warning_flag.template view<DeviceType>();
  h_warning_flag = k_warning_flag.h_view;

  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ImproperHarmonicKokkos<DeviceType>::~ImproperHarmonicKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperHarmonicKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  //atomKK->sync(execution_space,datamask_read);
  k_k.template sync<DeviceType>();
  k_chi.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  neighborKK->k_improperlist.template sync<DeviceType>();
  improperlist = neighborKK->k_improperlist.view<DeviceType>();
  int nimproperlist = neighborKK->nimproperlist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  h_warning_flag() = 0;
  k_warning_flag.template modify<LMPHostType>();
  k_warning_flag.template sync<DeviceType>();

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperHarmonicCompute<1,1> >(0,nimproperlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperHarmonicCompute<0,1> >(0,nimproperlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperHarmonicCompute<1,0> >(0,nimproperlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperHarmonicCompute<0,0> >(0,nimproperlist),*this);
    }
  }

  // error check

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.template sync<LMPHostType>();
  if (h_warning_flag())
    error->warning(FLERR,"Dihedral problem");

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
void ImproperHarmonicKokkos<DeviceType>::operator()(TagImproperHarmonicCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  const int i1 = improperlist(n,0);
  const int i2 = improperlist(n,1);
  const int i3 = improperlist(n,2);
  const int i4 = improperlist(n,3);
  const int type = improperlist(n,4);

  // geometry of 4-body

  const F_FLOAT vb1x = x(i1,0) - x(i2,0);
  const F_FLOAT vb1y = x(i1,1) - x(i2,1);
  const F_FLOAT vb1z = x(i1,2) - x(i2,2);

  const F_FLOAT vb2x = x(i3,0) - x(i2,0);
  const F_FLOAT vb2y = x(i3,1) - x(i2,1);
  const F_FLOAT vb2z = x(i3,2) - x(i2,2);

  const F_FLOAT vb3x = x(i4,0) - x(i3,0);
  const F_FLOAT vb3y = x(i4,1) - x(i3,1);
  const F_FLOAT vb3z = x(i4,2) - x(i3,2);

  const F_FLOAT ss1 = 1.0 / (vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  const F_FLOAT ss2 = 1.0 / (vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
  const F_FLOAT ss3 = 1.0 / (vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);

  const F_FLOAT r1 = sqrt(ss1);
  const F_FLOAT r2 = sqrt(ss2);
  const F_FLOAT r3 = sqrt(ss3);

  // sin and cos of improper

  const F_FLOAT c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
  const F_FLOAT c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
  const F_FLOAT c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

  F_FLOAT s1 = 1.0 - c1*c1;
  if (s1 < SMALL) s1 = SMALL;
  s1 = 1.0 / s1;

  F_FLOAT s2 = 1.0 - c2*c2;
  if (s2 < SMALL) s2 = SMALL;
  s2 = 1.0 / s2;

  F_FLOAT s12 = sqrt(s1*s2);
  F_FLOAT c = (c1*c2 + c0) * s12;

  // error check

  if ((c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) && !d_warning_flag())
    d_warning_flag() = 1;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  F_FLOAT s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;

  // force & energy

  const F_FLOAT domega = acos(c) - d_chi[type];
  F_FLOAT a = d_k[type] * domega;

  F_FLOAT eimproper = 0.0;
  if (eflag) eimproper = a*domega;

  a = -a * 2.0/s;
  c = c * a;
  s12 = s12 * a;
  const F_FLOAT a11 = c*ss1*s1;
  const F_FLOAT a22 = -ss2 * (2.0*c0*s12 - c*(s1+s2));
  const F_FLOAT a33 = c*ss3*s2;
  const F_FLOAT a12 = -r1*r2*(c1*c*s1 + c2*s12);
  const F_FLOAT a13 = -r1*r3*s12;
  const F_FLOAT a23 = r2*r3*(c2*c*s2 + c1*s12);

  const F_FLOAT sx2  = a22*vb2x + a23*vb3x + a12*vb1x;
  const F_FLOAT sy2  = a22*vb2y + a23*vb3y + a12*vb1y;
  const F_FLOAT sz2  = a22*vb2z + a23*vb3z + a12*vb1z;

  F_FLOAT f1[3],f2[3],f3[3],f4[3];
  f1[0] = a12*vb2x + a13*vb3x + a11*vb1x;
  f1[1] = a12*vb2y + a13*vb3y + a11*vb1y;
  f1[2] = a12*vb2z + a13*vb3z + a11*vb1z;

  f2[0] = -sx2 - f1[0];
  f2[1] = -sy2 - f1[1];
  f2[2] = -sz2 - f1[2];

  f4[0] = a23*vb2x + a33*vb3x + a13*vb1x;
  f4[1] = a23*vb2y + a33*vb3y + a13*vb1y;
  f4[2] = a23*vb2z + a33*vb3z + a13*vb1z;

  f3[0] = sx2 - f4[0];
  f3[1] = sy2 - f4[1];
  f3[2] = sz2 - f4[2];

  // apply force to each of 4 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    f(i1,0) += f1[0];
    f(i1,1) += f1[1];
    f(i1,2) += f1[2];
  }

  if (NEWTON_BOND || i2 < nlocal) {
    f(i2,0) += f2[0];
    f(i2,1) += f2[1];
    f(i2,2) += f2[2];
  }

  if (NEWTON_BOND || i3 < nlocal) {
    f(i3,0) += f3[0];
    f(i3,1) += f3[1];
    f(i3,2) += f3[2];
  }

  if (NEWTON_BOND || i4 < nlocal) {
    f(i4,0) += f4[0];
    f(i4,1) += f4[1];
    f(i4,2) += f4[2];
  }

  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,i4,eimproper,f1,f3,f4,
             vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void ImproperHarmonicKokkos<DeviceType>::operator()(TagImproperHarmonicCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagImproperHarmonicCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperHarmonicKokkos<DeviceType>::allocate()
{
  ImproperHarmonic::allocate();

  int n = atom->nimpropertypes;
  k_k = Kokkos::DualView<F_FLOAT*,DeviceType>("ImproperHarmonic::k",n+1);
  k_chi = Kokkos::DualView<F_FLOAT*,DeviceType>("ImproperHarmonic::chi",n+1);

  d_k = k_k.template view<DeviceType>();
  d_chi = k_chi.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void ImproperHarmonicKokkos<DeviceType>::coeff(int narg, char **arg)
{
  ImproperHarmonic::coeff(narg, arg);

  int n = atom->nimpropertypes;
  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_chi.h_view[i] = chi[i];
  }

  k_k.template modify<LMPHostType>();
  k_chi.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void ImproperHarmonicKokkos<DeviceType>::read_restart(FILE *fp)
{
  ImproperHarmonic::read_restart(fp);

  int n = atom->nimpropertypes;
  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_chi.h_view[i] = chi[i];
  }

  k_k.template modify<LMPHostType>();
  k_chi.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 + r4F4 = (r1-r2) F1 + (r3-r2) F3 + (r4-r2) F4
          = (r1-r2) F1 + (r3-r2) F3 + (r4-r3 + r3-r2) F4
          = vb1*f1 + vb2*f3 + (vb3+vb2)*f4
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
KOKKOS_INLINE_FUNCTION
void ImproperHarmonicKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                        F_FLOAT &eimproper, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                        const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                        const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                        const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const
{
  E_FLOAT eimproperquarter;
  F_FLOAT v[6];


  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += eimproper;
      else {
        eimproperquarter = 0.25*eimproper;
        if (i1 < nlocal) ev.evdwl += eimproperquarter;
        if (i2 < nlocal) ev.evdwl += eimproperquarter;
        if (i3 < nlocal) ev.evdwl += eimproperquarter;
        if (i4 < nlocal) ev.evdwl += eimproperquarter;
      }
    }
    if (eflag_atom) {
      eimproperquarter = 0.25*eimproper;
      if (newton_bond || i1 < nlocal) d_eatom[i1] += eimproperquarter;
      if (newton_bond || i2 < nlocal) d_eatom[i2] += eimproperquarter;
      if (newton_bond || i3 < nlocal) d_eatom[i3] += eimproperquarter;
      if (newton_bond || i4 < nlocal) d_eatom[i4] += eimproperquarter;
    }
  }

  if (vflag_either) {
    v[0] = vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0];
    v[1] = vb1y*f1[1] + vb2y*f3[1] + (vb3y+vb2y)*f4[1];
    v[2] = vb1z*f1[2] + vb2z*f3[2] + (vb3z+vb2z)*f4[2];
    v[3] = vb1x*f1[1] + vb2x*f3[1] + (vb3x+vb2x)*f4[1];
    v[4] = vb1x*f1[2] + vb2x*f3[2] + (vb3x+vb2x)*f4[2];
    v[5] = vb1y*f1[2] + vb2y*f3[2] + (vb3y+vb2y)*f4[2];

    if (vflag_global) {
      if (newton_bond) {
        ev.v[0] += v[0];
        ev.v[1] += v[1];
        ev.v[2] += v[2];
        ev.v[3] += v[3];
        ev.v[4] += v[4];
        ev.v[5] += v[5];
      } else {
        if (i1 < nlocal) {
          ev.v[0] += 0.25*v[0];
          ev.v[1] += 0.25*v[1];
          ev.v[2] += 0.25*v[2];
          ev.v[3] += 0.25*v[3];
          ev.v[4] += 0.25*v[4];
          ev.v[5] += 0.25*v[5];
        }
        if (i2 < nlocal) {
          ev.v[0] += 0.25*v[0];
          ev.v[1] += 0.25*v[1];
          ev.v[2] += 0.25*v[2];
          ev.v[3] += 0.25*v[3];
          ev.v[4] += 0.25*v[4];
          ev.v[5] += 0.25*v[5];
        }
        if (i3 < nlocal) {
          ev.v[0] += 0.25*v[0];
          ev.v[1] += 0.25*v[1];
          ev.v[2] += 0.25*v[2];
          ev.v[3] += 0.25*v[3];
          ev.v[4] += 0.25*v[4];
          ev.v[5] += 0.25*v[5];
        }
        if (i4 < nlocal) {
          ev.v[0] += 0.25*v[0];
          ev.v[1] += 0.25*v[1];
          ev.v[2] += 0.25*v[2];
          ev.v[3] += 0.25*v[3];
          ev.v[4] += 0.25*v[4];
          ev.v[5] += 0.25*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i1 < nlocal) {
        d_vatom(i1,0) += 0.25*v[0];
        d_vatom(i1,1) += 0.25*v[1];
        d_vatom(i1,2) += 0.25*v[2];
        d_vatom(i1,3) += 0.25*v[3];
        d_vatom(i1,4) += 0.25*v[4];
        d_vatom(i1,5) += 0.25*v[5];
      }
      if (newton_bond || i2 < nlocal) {
        d_vatom(i2,0) += 0.25*v[0];
        d_vatom(i2,1) += 0.25*v[1];
        d_vatom(i2,2) += 0.25*v[2];
        d_vatom(i2,3) += 0.25*v[3];
        d_vatom(i2,4) += 0.25*v[4];
        d_vatom(i2,5) += 0.25*v[5];
      }
      if (newton_bond || i3 < nlocal) {
        d_vatom(i3,0) += 0.25*v[0];
        d_vatom(i3,1) += 0.25*v[1];
        d_vatom(i3,2) += 0.25*v[2];
        d_vatom(i3,3) += 0.25*v[3];
        d_vatom(i3,4) += 0.25*v[4];
        d_vatom(i3,5) += 0.25*v[5];
      }
      if (newton_bond || i4 < nlocal) {
        d_vatom(i4,0) += 0.25*v[0];
        d_vatom(i4,1) += 0.25*v[1];
        d_vatom(i4,2) += 0.25*v[2];
        d_vatom(i4,3) += 0.25*v[3];
        d_vatom(i4,4) += 0.25*v[4];
        d_vatom(i4,5) += 0.25*v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ImproperHarmonicKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ImproperHarmonicKokkos<LMPHostType>;
#endif
}

