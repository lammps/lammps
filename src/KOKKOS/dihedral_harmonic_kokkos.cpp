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
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "dihedral_harmonic_kokkos.h"
#include <cmath>
#include <cstdlib>
#include "atom_kokkos.h"
#include "comm.h"
#include "neighbor_kokkos.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

#define TOLERANCE 0.05
#define SMALL     0.001
#define SMALLER   0.00001

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralHarmonicKokkos<DeviceType>::DihedralHarmonicKokkos(LAMMPS *lmp) : DihedralHarmonic(lmp)
{
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_warning_flag = DAT::tdual_int_scalar("Dihedral:warning_flag");
  d_warning_flag = k_warning_flag.view<DeviceType>();
  h_warning_flag = k_warning_flag.h_view;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
DihedralHarmonicKokkos<DeviceType>::~DihedralHarmonicKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralHarmonicKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"dihedral:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,6,"dihedral:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  k_k.template sync<DeviceType>();
  k_cos_shift.template sync<DeviceType>();
  k_sin_shift.template sync<DeviceType>();
  k_sign.template sync<DeviceType>();
  k_multiplicity.template sync<DeviceType>();

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  neighborKK->k_dihedrallist.template sync<DeviceType>();
  dihedrallist = neighborKK->k_dihedrallist.view<DeviceType>();
  int ndihedrallist = neighborKK->ndihedrallist;
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
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralHarmonicCompute<1,1> >(0,ndihedrallist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagDihedralHarmonicCompute<0,1> >(0,ndihedrallist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralHarmonicCompute<1,0> >(0,ndihedrallist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagDihedralHarmonicCompute<0,0> >(0,ndihedrallist),*this);
    }
  }

  // error check

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.template sync<LMPHostType>();
  if (h_warning_flag())
    error->warning(FLERR,"Dihedral problem",0);

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
void DihedralHarmonicKokkos<DeviceType>::operator()(TagDihedralHarmonicCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  // The f array is atomic
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = dihedrallist(n,0);
  const int i2 = dihedrallist(n,1);
  const int i3 = dihedrallist(n,2);
  const int i4 = dihedrallist(n,3);
  const int type = dihedrallist(n,4);

  // 1st bond

  const F_FLOAT vb1x = x(i1,0) - x(i2,0);
  const F_FLOAT vb1y = x(i1,1) - x(i2,1);
  const F_FLOAT vb1z = x(i1,2) - x(i2,2);

  // 2nd bond

  const F_FLOAT vb2x = x(i3,0) - x(i2,0);
  const F_FLOAT vb2y = x(i3,1) - x(i2,1);
  const F_FLOAT vb2z = x(i3,2) - x(i2,2);

  const F_FLOAT vb2xm = -vb2x;
  const F_FLOAT vb2ym = -vb2y;
  const F_FLOAT vb2zm = -vb2z;

  // 3rd bond

  const F_FLOAT vb3x = x(i4,0) - x(i3,0);
  const F_FLOAT vb3y = x(i4,1) - x(i3,1);
  const F_FLOAT vb3z = x(i4,2) - x(i3,2);

  // c,s calculation

  const F_FLOAT ax = vb1y*vb2zm - vb1z*vb2ym;
  const F_FLOAT ay = vb1z*vb2xm - vb1x*vb2zm;
  const F_FLOAT az = vb1x*vb2ym - vb1y*vb2xm;
  const F_FLOAT bx = vb3y*vb2zm - vb3z*vb2ym;
  const F_FLOAT by = vb3z*vb2xm - vb3x*vb2zm;
  const F_FLOAT bz = vb3x*vb2ym - vb3y*vb2xm;

  const F_FLOAT rasq = ax*ax + ay*ay + az*az;
  const F_FLOAT rbsq = bx*bx + by*by + bz*bz;
  const F_FLOAT rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
  const F_FLOAT rg = sqrt(rgsq);

  F_FLOAT rginv,ra2inv,rb2inv;
  rginv = ra2inv = rb2inv = 0.0;
  if (rg > 0) rginv = 1.0/rg;
  if (rasq > 0) ra2inv = 1.0/rasq;
  if (rbsq > 0) rb2inv = 1.0/rbsq;
  const F_FLOAT rabinv = sqrt(ra2inv*rb2inv);

  F_FLOAT c = (ax*bx + ay*by + az*bz)*rabinv;
  const F_FLOAT s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

  // error check

  if ((c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) && !d_warning_flag())
    Kokkos::atomic_fetch_add(&d_warning_flag(),1);

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  const int m = d_multiplicity[type];
  F_FLOAT p = 1.0;
  F_FLOAT ddf1,df1;
  ddf1 = df1 = 0.0;

  for (int i = 0; i < m; i++) {
    ddf1 = p*c - df1*s;
    df1 = p*s + df1*c;
    p = ddf1;
  }

  p = p*d_cos_shift[type] + df1*d_sin_shift[type];
  df1 = df1*d_cos_shift[type] - ddf1*d_sin_shift[type];
  df1 *= -m;
  p += 1.0;

  if (m == 0) {
    p = 1.0 + d_cos_shift[type];
    df1 = 0.0;
  }

  E_FLOAT edihedral = 0.0;
  if (eflag) edihedral = d_k[type] * p;

  const F_FLOAT fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
  const F_FLOAT hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
  const F_FLOAT fga = fg*ra2inv*rginv;
  const F_FLOAT hgb = hg*rb2inv*rginv;
  const F_FLOAT gaa = -ra2inv*rg;
  const F_FLOAT gbb = rb2inv*rg;

  const F_FLOAT dtfx = gaa*ax;
  const F_FLOAT dtfy = gaa*ay;
  const F_FLOAT dtfz = gaa*az;
  const F_FLOAT dtgx = fga*ax - hgb*bx;
  const F_FLOAT dtgy = fga*ay - hgb*by;
  const F_FLOAT dtgz = fga*az - hgb*bz;
  const F_FLOAT dthx = gbb*bx;
  const F_FLOAT dthy = gbb*by;
  const F_FLOAT dthz = gbb*bz;

  const F_FLOAT df = -d_k[type] * df1;

  const F_FLOAT sx2  = df*dtgx;;
  const F_FLOAT sy2  = df*dtgy;;
  const F_FLOAT sz2  = df*dtgz;;

  F_FLOAT f1[3],f2[3],f3[3],f4[3];
  f1[0] = df*dtfx;
  f1[1] = df*dtfy;
  f1[2] = df*dtfz;

  f2[0] = sx2 - f1[0];
  f2[1] = sy2 - f1[1];
  f2[2] = sz2 - f1[2];

  f4[0] = df*dthx;
  f4[1] = df*dthy;
  f4[2] = df*dthz;

  f3[0] = -sx2 - f4[0];
  f3[1] = -sy2 - f4[1];
  f3[2] = -sz2 - f4[2];

  // apply force to each of 4 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    a_f(i1,0) += f1[0];
    a_f(i1,1) += f1[1];
    a_f(i1,2) += f1[2];
  }

  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) += f2[0];
    a_f(i2,1) += f2[1];
    a_f(i2,2) += f2[2];
  }

  if (NEWTON_BOND || i3 < nlocal) {
    a_f(i3,0) += f3[0];
    a_f(i3,1) += f3[1];
    a_f(i3,2) += f3[2];
  }

  if (NEWTON_BOND || i4 < nlocal) {
    a_f(i4,0) += f4[0];
    a_f(i4,1) += f4[1];
    a_f(i4,2) += f4[2];
  }

  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,i4,edihedral,f1,f3,f4,
             vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void DihedralHarmonicKokkos<DeviceType>::operator()(TagDihedralHarmonicCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagDihedralHarmonicCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void DihedralHarmonicKokkos<DeviceType>::allocate()
{
  DihedralHarmonic::allocate();

  int n = atom->ndihedraltypes;
  k_k = DAT::tdual_ffloat_1d("DihedralHarmonic::k",n+1);
  k_cos_shift = DAT::tdual_ffloat_1d("DihedralHarmonic::cos_shift",n+1);
  k_sin_shift = DAT::tdual_ffloat_1d("DihedralHarmonic::sin_shift",n+1);
  k_sign = DAT::tdual_int_1d("DihedralHarmonic::sign",n+1);
  k_multiplicity = DAT::tdual_int_1d("DihedralHarmonic::multiplicity",n+1);

  d_k = k_k.template view<DeviceType>();
  d_cos_shift = k_cos_shift.template view<DeviceType>();
  d_sin_shift = k_sin_shift.template view<DeviceType>();
  d_sign = k_sign.template view<DeviceType>();
  d_multiplicity = k_multiplicity.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralHarmonicKokkos<DeviceType>::coeff(int narg, char **arg)
{
  DihedralHarmonic::coeff(narg, arg);

  int n = atom->ndihedraltypes;
  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_cos_shift.h_view[i] = cos_shift[i];
    k_sin_shift.h_view[i] = sin_shift[i];
    k_sign.h_view[i] = sign[i];
    k_multiplicity.h_view[i] = multiplicity[i];
  }

  k_k.template modify<LMPHostType>();
  k_cos_shift.template modify<LMPHostType>();
  k_sin_shift.template modify<LMPHostType>();
  k_sign.template modify<LMPHostType>();
  k_multiplicity.template modify<LMPHostType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void DihedralHarmonicKokkos<DeviceType>::read_restart(FILE *fp)
{
  DihedralHarmonic::read_restart(fp);

  int n = atom->ndihedraltypes;
  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_cos_shift.h_view[i] = cos_shift[i];
    k_sin_shift.h_view[i] = sin_shift[i];
    k_sign.h_view[i] = sign[i];
    k_multiplicity.h_view[i] = multiplicity[i];
  }

  k_k.template modify<LMPHostType>();
  k_cos_shift.template modify<LMPHostType>();
  k_sin_shift.template modify<LMPHostType>();
  k_sign.template modify<LMPHostType>();
  k_multiplicity.template modify<LMPHostType>();
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
void DihedralHarmonicKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                        F_FLOAT &edihedral, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                        const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                        const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                        const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const
{
  E_FLOAT edihedralquarter;
  F_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = k_vatom.view<DeviceType>();

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += edihedral;
      else {
        edihedralquarter = 0.25*edihedral;
        if (i1 < nlocal) ev.evdwl += edihedralquarter;
        if (i2 < nlocal) ev.evdwl += edihedralquarter;
        if (i3 < nlocal) ev.evdwl += edihedralquarter;
        if (i4 < nlocal) ev.evdwl += edihedralquarter;
      }
    }
    if (eflag_atom) {
      edihedralquarter = 0.25*edihedral;
      if (newton_bond || i1 < nlocal) v_eatom[i1] += edihedralquarter;
      if (newton_bond || i2 < nlocal) v_eatom[i2] += edihedralquarter;
      if (newton_bond || i3 < nlocal) v_eatom[i3] += edihedralquarter;
      if (newton_bond || i4 < nlocal) v_eatom[i4] += edihedralquarter;
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
        v_vatom(i1,0) += 0.25*v[0];
        v_vatom(i1,1) += 0.25*v[1];
        v_vatom(i1,2) += 0.25*v[2];
        v_vatom(i1,3) += 0.25*v[3];
        v_vatom(i1,4) += 0.25*v[4];
        v_vatom(i1,5) += 0.25*v[5];
      }
      if (newton_bond || i2 < nlocal) {
        v_vatom(i2,0) += 0.25*v[0];
        v_vatom(i2,1) += 0.25*v[1];
        v_vatom(i2,2) += 0.25*v[2];
        v_vatom(i2,3) += 0.25*v[3];
        v_vatom(i2,4) += 0.25*v[4];
        v_vatom(i2,5) += 0.25*v[5];
      }
      if (newton_bond || i3 < nlocal) {
        v_vatom(i3,0) += 0.25*v[0];
        v_vatom(i3,1) += 0.25*v[1];
        v_vatom(i3,2) += 0.25*v[2];
        v_vatom(i3,3) += 0.25*v[3];
        v_vatom(i3,4) += 0.25*v[4];
        v_vatom(i3,5) += 0.25*v[5];
      }
      if (newton_bond || i4 < nlocal) {
        v_vatom(i4,0) += 0.25*v[0];
        v_vatom(i4,1) += 0.25*v[1];
        v_vatom(i4,2) += 0.25*v[2];
        v_vatom(i4,3) += 0.25*v[3];
        v_vatom(i4,4) += 0.25*v[4];
        v_vatom(i4,5) += 0.25*v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class DihedralHarmonicKokkos<LMPDeviceType>;
#ifdef KOKKOS_ENABLE_CUDA
template class DihedralHarmonicKokkos<LMPHostType>;
#endif
}

