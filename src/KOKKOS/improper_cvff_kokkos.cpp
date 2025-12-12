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

#include "improper_cvff_kokkos.h"
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
ImproperCvffKokkos<DeviceType>::ImproperCvffKokkos(LAMMPS *lmp) : ImproperCvff(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  k_warning_flag = DAT::tdual_int_scalar("Improper:warning_flag");
  d_warning_flag = k_warning_flag.template view<DeviceType>();
  h_warning_flag = k_warning_flag.view_host();

  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ImproperCvffKokkos<DeviceType>::~ImproperCvffKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperCvffKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    if ((int)k_eatom.extent(0) < maxeatom) {
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"improper:eatom");
      d_eatom = k_eatom.template view<DeviceType>();
    } else Kokkos::deep_copy(d_eatom,0.0);
  }
  if (vflag_atom) {
    if ((int)k_vatom.extent(0) < maxvatom) {
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"improper:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    } else Kokkos::deep_copy(d_vatom,0.0);
  }

  //atomKK->sync(execution_space,datamask_read);
  k_k.template sync<DeviceType>();
  k_sign.template sync<DeviceType>();
  k_multiplicity.template sync<DeviceType>();
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
  k_warning_flag.modify_host();
  k_warning_flag.template sync<DeviceType>();

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperCvffCompute<1,1> >(0,nimproperlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperCvffCompute<0,1> >(0,nimproperlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperCvffCompute<1,0> >(0,nimproperlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperCvffCompute<0,0> >(0,nimproperlist),*this);
    }
  }

  // error check

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.sync_host();
  if (h_warning_flag())
    error->warning(FLERR,"ImproperCvff problem");

  if (eflag_global) energy += static_cast<double>(ev.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  copymode = 0;
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void ImproperCvffKokkos<DeviceType>::operator()(TagImproperCvffCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  const int i1 = improperlist(n,0);
  const int i2 = improperlist(n,1);
  const int i3 = improperlist(n,2);
  const int i4 = improperlist(n,3);
  const int type = improperlist(n,4);

  // 1st bond

  const KK_FLOAT vb1x = x(i1,0) - x(i2,0);
  const KK_FLOAT vb1y = x(i1,1) - x(i2,1);
  const KK_FLOAT vb1z = x(i1,2) - x(i2,2);

  // 2nd bond

  const KK_FLOAT vb2x = x(i3,0) - x(i2,0);
  const KK_FLOAT vb2y = x(i3,1) - x(i2,1);
  const KK_FLOAT vb2z = x(i3,2) - x(i2,2);

  const KK_FLOAT vb2xm = -vb2x;
  const KK_FLOAT vb2ym = -vb2y;
  const KK_FLOAT vb2zm = -vb2z;

  // 3rd bond

  const KK_FLOAT vb3x = x(i4,0) - x(i3,0);
  const KK_FLOAT vb3y = x(i4,1) - x(i3,1);
  const KK_FLOAT vb3z = x(i4,2) - x(i3,2);

  // c0 calculation

  const KK_FLOAT sb1 = static_cast<KK_FLOAT>(1.0) / (vb1x * vb1x + vb1y * vb1y + vb1z * vb1z);
  const KK_FLOAT sb2 = static_cast<KK_FLOAT>(1.0) / (vb2x * vb2x + vb2y * vb2y + vb2z * vb2z);
  const KK_FLOAT sb3 = static_cast<KK_FLOAT>(1.0) / (vb3x * vb3x + vb3y * vb3y + vb3z * vb3z);

  const KK_FLOAT rb1 = sqrt(sb1);
  const KK_FLOAT rb3 = sqrt(sb3);

  const KK_FLOAT c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * rb1 * rb3;

  // 1st and 2nd angle

  const KK_FLOAT b1mag2 = vb1x * vb1x + vb1y * vb1y + vb1z * vb1z;
  const KK_FLOAT b1mag = sqrt(b1mag2);
  const KK_FLOAT b2mag2 = vb2x * vb2x + vb2y * vb2y + vb2z * vb2z;
  const KK_FLOAT b2mag = sqrt(b2mag2);
  const KK_FLOAT b3mag2 = vb3x * vb3x + vb3y * vb3y + vb3z * vb3z;
  const KK_FLOAT b3mag = sqrt(b3mag2);

  KK_FLOAT ctmp = vb1x * vb2x + vb1y * vb2y + vb1z * vb2z;
  const KK_FLOAT r12c1 = static_cast<KK_FLOAT>(1.0) / (b1mag * b2mag);
  const KK_FLOAT c1mag = ctmp * r12c1;

  ctmp = vb2xm * vb3x + vb2ym * vb3y + vb2zm * vb3z;
  const KK_FLOAT r12c2 = static_cast<KK_FLOAT>(1.0) / (b2mag * b3mag);
  const KK_FLOAT c2mag = ctmp * r12c2;

  // cos and sin of 2 angles and final c

  KK_FLOAT sc1 = sqrt(static_cast<KK_FLOAT>(1.0) - c1mag * c1mag);
  if (sc1 <  static_cast<KK_FLOAT>(SMALL)) sc1 = static_cast<KK_FLOAT>(SMALL);
  sc1 = static_cast<KK_FLOAT>(1.0) / sc1;

  KK_FLOAT sc2 = sqrt(static_cast<KK_FLOAT>(1.0) - c2mag * c2mag);
  if (sc2 <  static_cast<KK_FLOAT>(SMALL)) sc2 = static_cast<KK_FLOAT>(SMALL);
  sc2 = static_cast<KK_FLOAT>(1.0) / sc2;

  const KK_FLOAT s1 = sc1 * sc1;
  const KK_FLOAT s2 = sc2 * sc2;
  KK_FLOAT s12 = sc1 * sc2;
  KK_FLOAT c = (c0 + c1mag * c2mag) * s12;

  // error check

  if ((c > static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(TOLERANCE) || c < (-static_cast<KK_FLOAT>(1.0) - static_cast<KK_FLOAT>(TOLERANCE))) && !d_warning_flag())
    d_warning_flag() = 1;

  if (c > static_cast<KK_FLOAT>(1.0)) c = static_cast<KK_FLOAT>(1.0);
  if (c < -static_cast<KK_FLOAT>(1.0)) c = -static_cast<KK_FLOAT>(1.0);

  // force & energy
  // p = 1 + cos(n*phi) for d = 1
  // p = 1 - cos(n*phi) for d = -1
  // pd = dp/dc / 2

  const int m = d_multiplicity[type];

  KK_FLOAT p,pd,rc2;

  if (m == 2) {
    p = static_cast<KK_FLOAT>(2.0) * c * c;
    pd = static_cast<KK_FLOAT>(2.0) * c;
  } else if (m == 3) {
    rc2 = c * c;
    p = (static_cast<KK_FLOAT>(4.0) * rc2 - static_cast<KK_FLOAT>(3.0)) * c + static_cast<KK_FLOAT>(1.0);
    pd = static_cast<KK_FLOAT>(6.0) * rc2 - static_cast<KK_FLOAT>(1.5);
  } else if (m == 4) {
    rc2 = c * c;
    p = static_cast<KK_FLOAT>(8.0) * (rc2 - static_cast<KK_FLOAT>(1.0)) * rc2 + static_cast<KK_FLOAT>(2.0);
    pd = (static_cast<KK_FLOAT>(16.0) * rc2 - static_cast<KK_FLOAT>(8.0)) * c;
  } else if (m == 6) {
    rc2 = c * c;
    p = ((static_cast<KK_FLOAT>(32.0) * rc2 - static_cast<KK_FLOAT>(48.0)) * rc2 + static_cast<KK_FLOAT>(18.0)) * rc2;
    pd = (static_cast<KK_FLOAT>(96.0) * (rc2 - static_cast<KK_FLOAT>(1.0)) * rc2 + static_cast<KK_FLOAT>(18.0)) * c;
  } else if (m == 1) {
    p = c + static_cast<KK_FLOAT>(1.0);
    pd = static_cast<KK_FLOAT>(0.5);
  } else if (m == 5) {
    rc2 = c * c;
    p = ((static_cast<KK_FLOAT>(16.0) * rc2 - static_cast<KK_FLOAT>(20.0)) * rc2 + static_cast<KK_FLOAT>(5.0)) * c + static_cast<KK_FLOAT>(1.0);
    pd = (static_cast<KK_FLOAT>(40.0) * rc2 - static_cast<KK_FLOAT>(30.0)) * rc2 + static_cast<KK_FLOAT>(2.5);
  } else if (m == 0) {
    p = static_cast<KK_FLOAT>(2.0);
    pd = static_cast<KK_FLOAT>(0.0);
  }

  if (sign[type] == -1) {
    p = static_cast<KK_FLOAT>(2.0) - p;
    pd = -pd;
  }

  KK_FLOAT eimproper = 0;
  if (eflag) eimproper = d_k[type] * p;

  const KK_FLOAT a = static_cast<KK_FLOAT>(2.0) * d_k[type] * pd;
  c = c * a;
  s12 = s12 * a;
  const KK_FLOAT a11 = c * sb1 * s1;
  const KK_FLOAT a22 = -sb2 * (static_cast<KK_FLOAT>(2.0) * c0 * s12 - c * (s1 + s2));
  const KK_FLOAT a33 = c * sb3 * s2;
  const KK_FLOAT a12 = -r12c1 * (c1mag * c * s1 + c2mag * s12);
  const KK_FLOAT a13 = -rb1 * rb3 * s12;
  const KK_FLOAT a23 = r12c2 * (c2mag * c * s2 + c1mag * s12);

  const KK_FLOAT sx2 = a12 * vb1x + a22 * vb2x + a23 * vb3x;
  const KK_FLOAT sy2 = a12 * vb1y + a22 * vb2y + a23 * vb3y;
  const KK_FLOAT sz2 = a12 * vb1z + a22 * vb2z + a23 * vb3z;

  KK_FLOAT f1[3],f2[3],f3[3],f4[3];
  f1[0] = a11 * vb1x + a12 * vb2x + a13 * vb3x;
  f1[1] = a11 * vb1y + a12 * vb2y + a13 * vb3y;
  f1[2] = a11 * vb1z + a12 * vb2z + a13 * vb3z;

  f2[0] = -sx2 - f1[0];
  f2[1] = -sy2 - f1[1];
  f2[2] = -sz2 - f1[2];

  f4[0] = a13 * vb1x + a23 * vb2x + a33 * vb3x;
  f4[1] = a13 * vb1y + a23 * vb2y + a33 * vb3y;
  f4[2] = a13 * vb1z + a23 * vb2z + a33 * vb3z;

  f3[0] = sx2 - f4[0];
  f3[1] = sy2 - f4[1];
  f3[2] = sz2 - f4[2];

  // apply force to each of 4 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    f(i1,0) += static_cast<KK_ACC_FLOAT>(f1[0]);
    f(i1,1) += static_cast<KK_ACC_FLOAT>(f1[1]);
    f(i1,2) += static_cast<KK_ACC_FLOAT>(f1[2]);
  }

  if (NEWTON_BOND || i2 < nlocal) {
    f(i2,0) += static_cast<KK_ACC_FLOAT>(f2[0]);
    f(i2,1) += static_cast<KK_ACC_FLOAT>(f2[1]);
    f(i2,2) += static_cast<KK_ACC_FLOAT>(f2[2]);
  }

  if (NEWTON_BOND || i3 < nlocal) {
    f(i3,0) += static_cast<KK_ACC_FLOAT>(f3[0]);
    f(i3,1) += static_cast<KK_ACC_FLOAT>(f3[1]);
    f(i3,2) += static_cast<KK_ACC_FLOAT>(f3[2]);
  }

  if (NEWTON_BOND || i4 < nlocal) {
    f(i4,0) += static_cast<KK_ACC_FLOAT>(f4[0]);
    f(i4,1) += static_cast<KK_ACC_FLOAT>(f4[1]);
    f(i4,2) += static_cast<KK_ACC_FLOAT>(f4[2]);
  }

  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,i4,eimproper,f1,f3,f4,
             vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void ImproperCvffKokkos<DeviceType>::operator()(TagImproperCvffCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagImproperCvffCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperCvffKokkos<DeviceType>::allocate()
{
  ImproperCvff::allocate();

  int n = atom->nimpropertypes;
  k_k = DAT::tdual_kkfloat_1d("ImproperCvff::k",n+1);
  k_sign = DAT::tdual_int_1d("ImproperCvff::sign",n+1);
  k_multiplicity = DAT::tdual_int_1d("ImproperCvff::multiplicity",n+1);

  d_k = k_k.template view<DeviceType>();
  d_sign = k_sign.template view<DeviceType>();
  d_multiplicity = k_multiplicity.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void ImproperCvffKokkos<DeviceType>::coeff(int narg, char **arg)
{
  ImproperCvff::coeff(narg, arg);

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nimpropertypes,ilo,ihi,error);

  for (int i = ilo; i <= ihi; i++) {
    k_k.view_host()[i] = static_cast<KK_FLOAT>(k[i]);
    k_sign.view_host()[i] = sign[i];
    k_multiplicity.view_host()[i] = multiplicity[i];
  }

  k_k.modify_host();
  k_sign.modify_host();
  k_multiplicity.modify_host();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void ImproperCvffKokkos<DeviceType>::read_restart(FILE *fp)
{
  ImproperCvff::read_restart(fp);

  int n = atom->nimpropertypes;
  for (int i = 1; i <= n; i++) {
    k_k.view_host()[i] = static_cast<KK_FLOAT>(k[i]);
    k_sign.view_host()[i] = sign[i];
    k_multiplicity.view_host()[i] = multiplicity[i];
  }

  k_k.modify_host();
  k_sign.modify_host();
  k_multiplicity.modify_host();
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
void ImproperCvffKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                        KK_FLOAT &eimproper, KK_FLOAT *f1, KK_FLOAT *f3, KK_FLOAT *f4,
                        const KK_FLOAT &vb1x, const KK_FLOAT &vb1y, const KK_FLOAT &vb1z,
                        const KK_FLOAT &vb2x, const KK_FLOAT &vb2y, const KK_FLOAT &vb2z,
                        const KK_FLOAT &vb3x, const KK_FLOAT &vb3y, const KK_FLOAT &vb3z) const
{
  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(eimproper);
      else {
        KK_ACC_FLOAT eimproperquarter = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*eimproper);
        if (i1 < nlocal) ev.evdwl += eimproperquarter;
        if (i2 < nlocal) ev.evdwl += eimproperquarter;
        if (i3 < nlocal) ev.evdwl += eimproperquarter;
        if (i4 < nlocal) ev.evdwl += eimproperquarter;
      }
    }
    if (eflag_atom) {
      KK_ACC_FLOAT eimproperquarter = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*eimproper);
      if (newton_bond || i1 < nlocal) d_eatom[i1] += eimproperquarter;
      if (newton_bond || i2 < nlocal) d_eatom[i2] += eimproperquarter;
      if (newton_bond || i3 < nlocal) d_eatom[i3] += eimproperquarter;
      if (newton_bond || i4 < nlocal) d_eatom[i4] += eimproperquarter;
    }
  }

  if (vflag_either) {
    KK_ACC_FLOAT v_quarter_acc[6];
    v_quarter_acc[0] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0]));
    v_quarter_acc[1] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1y*f1[1] + vb2y*f3[1] + (vb3y+vb2y)*f4[1]));
    v_quarter_acc[2] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1z*f1[2] + vb2z*f3[2] + (vb3z+vb2z)*f4[2]));
    v_quarter_acc[3] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1x*f1[1] + vb2x*f3[1] + (vb3x+vb2x)*f4[1]));
    v_quarter_acc[4] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1x*f1[2] + vb2x*f3[2] + (vb3x+vb2x)*f4[2]));
    v_quarter_acc[5] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1y*f1[2] + vb2y*f3[2] + (vb3y+vb2y)*f4[2]));

    if (vflag_global) {
      if (newton_bond) {
        for (int n = 0; n < 6; n++)
          ev.v[n] += static_cast<KK_ACC_FLOAT>(4.0)*v_quarter_acc[n];
      } else {
        if (i1 < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_quarter_acc[n];
        }
        if (i2 < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_quarter_acc[n];
        }
        if (i3 < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_quarter_acc[n];
        }
        if (i4 < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_quarter_acc[n];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i1 < nlocal) {
        for (int n = 0; n < 6; n++)
          d_vatom(i1,n) += v_quarter_acc[n];
      }
      if (newton_bond || i2 < nlocal) {
        for (int n = 0; n < 6; n++)
          d_vatom(i2,n) += v_quarter_acc[n];
      }
      if (newton_bond || i3 < nlocal) {
        for (int n = 0; n < 6; n++)
          d_vatom(i3,n) += v_quarter_acc[n];
      }
      if (newton_bond || i4 < nlocal) {
        for (int n = 0; n < 6; n++)
          d_vatom(i4,n) += v_quarter_acc[n];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ImproperCvffKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ImproperCvffKokkos<LMPHostType>;
#endif
}

