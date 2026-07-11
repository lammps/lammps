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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "improper_fourier_kokkos.h"
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
ImproperFourierKokkos<DeviceType>::ImproperFourierKokkos(LAMMPS *lmp) : ImproperFourier(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ImproperFourierKokkos<DeviceType>::~ImproperFourierKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperFourierKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

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

  k_k.template sync<DeviceType>();
  k_C0.template sync<DeviceType>();
  k_C1.template sync<DeviceType>();
  k_C2.template sync<DeviceType>();
  k_all.template sync<DeviceType>();
  atomKK->sync(execution_space,datamask_read);
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  neighborKK->k_improperlist.template sync<DeviceType>();
  improperlist = neighborKK->k_improperlist.view<DeviceType>();
  int nimproperlist = neighborKK->nimproperlist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  // zero warning flag
  k_warning_flag = DAT::tdual_int_scalar("ImproperFourier::warning_flag");
  d_warning_flag = k_warning_flag.template view<DeviceType>();
  h_warning_flag = k_warning_flag.view_host();
  h_warning_flag() = 0;
  k_warning_flag.modify_host();
  k_warning_flag.template sync<DeviceType>();

  copymode = 1;

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperFourierCompute<1,1> >(0,nimproperlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperFourierCompute<0,1> >(0,nimproperlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperFourierCompute<1,0> >(0,nimproperlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperFourierCompute<0,0> >(0,nimproperlist),*this);
    }
  }

  k_warning_flag.template modify<DeviceType>();
  k_warning_flag.sync_host();
  if (h_warning_flag())
    error->warning(FLERR,"Improper problem: see log");

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

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperFourierKokkos<DeviceType>::operator()(TagImproperFourierCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const
{
  const int i1 = improperlist(n,0);
  const int i2 = improperlist(n,1);
  const int i3 = improperlist(n,2);
  const int i4 = improperlist(n,3);
  const int type = improperlist(n,4);

  const KK_FLOAT vb1x = x(i2,0) - x(i1,0);
  const KK_FLOAT vb1y = x(i2,1) - x(i1,1);
  const KK_FLOAT vb1z = x(i2,2) - x(i1,2);

  const KK_FLOAT vb2x = x(i3,0) - x(i1,0);
  const KK_FLOAT vb2y = x(i3,1) - x(i1,1);
  const KK_FLOAT vb2z = x(i3,2) - x(i1,2);

  const KK_FLOAT vb3x = x(i4,0) - x(i1,0);
  const KK_FLOAT vb3y = x(i4,1) - x(i1,1);
  const KK_FLOAT vb3z = x(i4,2) - x(i1,2);

  addone<NEWTON_BOND,EVFLAG>(ev, i1,i2,i3,i4, type,
                              vb1x,vb1y,vb1z, vb2x,vb2y,vb2z, vb3x,vb3y,vb3z);
  if (d_all[type]) {
    addone<NEWTON_BOND,EVFLAG>(ev, i1,i4,i2,i3, type,
                                vb3x,vb3y,vb3z, vb1x,vb1y,vb1z, vb2x,vb2y,vb2z);
    addone<NEWTON_BOND,EVFLAG>(ev, i1,i3,i4,i2, type,
                                vb2x,vb2y,vb2z, vb3x,vb3y,vb3z, vb1x,vb1y,vb1z);
  }
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperFourierKokkos<DeviceType>::operator()(TagImproperFourierCompute<NEWTON_BOND,EVFLAG>, const int &n) const
{
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagImproperFourierCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperFourierKokkos<DeviceType>::addone(EV_FLOAT &ev,
    const int i1, const int i2, const int i3, const int i4,
    const int type,
    const KK_FLOAT vb1x, const KK_FLOAT vb1y, const KK_FLOAT vb1z,
    const KK_FLOAT vb2x, const KK_FLOAT vb2y, const KK_FLOAT vb2z,
    const KK_FLOAT vb3x, const KK_FLOAT vb3y, const KK_FLOAT vb3z) const
{
  KK_FLOAT eimproper,f1[3],f2[3],f3[3],f4[3];
  KK_FLOAT c,c2,a,s,projhfg,dhax,dhay,dhaz,dahx,dahy,dahz,cotphi;
  KK_FLOAT ax,ay,az,ra2,rh2,ra,rh,rar,rhr,arx,ary,arz,hrx,hry,hrz;

  eimproper = static_cast<KK_FLOAT>(0.0);

  ax = vb1y*vb2z - vb1z*vb2y;
  ay = vb1z*vb2x - vb1x*vb2z;
  az = vb1x*vb2y - vb1y*vb2x;
  ra2 = ax*ax + ay*ay + az*az;
  rh2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
  ra = sqrt(ra2);
  rh = sqrt(rh2);
  if (ra < static_cast<KK_FLOAT>(SMALL)) ra = static_cast<KK_FLOAT>(SMALL);
  if (rh < static_cast<KK_FLOAT>(SMALL)) rh = static_cast<KK_FLOAT>(SMALL);

  rar = static_cast<KK_FLOAT>(1.0)/ra;
  rhr = static_cast<KK_FLOAT>(1.0)/rh;
  arx = ax*rar;
  ary = ay*rar;
  arz = az*rar;
  hrx = vb3x*rhr;
  hry = vb3y*rhr;
  hrz = vb3z*rhr;

  c = arx*hrx + ary*hry + arz*hrz;

  if (c > static_cast<KK_FLOAT>(1.0) + static_cast<KK_FLOAT>(TOLERANCE) ||
      c < static_cast<KK_FLOAT>(-1.0) - static_cast<KK_FLOAT>(TOLERANCE))
    d_warning_flag() = 1;

  if (c > static_cast<KK_FLOAT>(1.0)) c = static_cast<KK_FLOAT>(1.0);
  if (c < static_cast<KK_FLOAT>(-1.0)) c = static_cast<KK_FLOAT>(-1.0);

  s = sqrt(static_cast<KK_FLOAT>(1.0) - c*c);
  if (s < static_cast<KK_FLOAT>(SMALL)) s = static_cast<KK_FLOAT>(SMALL);
  cotphi = c/s;

  projhfg = (vb3x*vb1x + vb3y*vb1y + vb3z*vb1z) /
    sqrt(vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  projhfg += (vb3x*vb2x + vb3y*vb2y + vb3z*vb2z) /
    sqrt(vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
  if (projhfg > static_cast<KK_FLOAT>(0.0)) {
    s *= static_cast<KK_FLOAT>(-1.0);
    cotphi *= static_cast<KK_FLOAT>(-1.0);
  }

  c2 = static_cast<KK_FLOAT>(2.0)*s*s - static_cast<KK_FLOAT>(1.0);
  if (EVFLAG && eflag)
    eimproper = d_k[type]*(d_C0[type] + d_C1[type]*s + d_C2[type]*c2);

  a = d_k[type]*(d_C1[type] + static_cast<KK_FLOAT>(4.0)*d_C2[type]*s)*cotphi;
  dhax = hrx - c*arx;
  dhay = hry - c*ary;
  dhaz = hrz - c*arz;

  dahx = arx - c*hrx;
  dahy = ary - c*hry;
  dahz = arz - c*hrz;

  f2[0] = (dhay*vb1z - dhaz*vb1y)*rar*a;
  f2[1] = (dhaz*vb1x - dhax*vb1z)*rar*a;
  f2[2] = (dhax*vb1y - dhay*vb1x)*rar*a;

  f3[0] = (-dhay*vb2z + dhaz*vb2y)*rar*a;
  f3[1] = (-dhaz*vb2x + dhax*vb2z)*rar*a;
  f3[2] = (-dhax*vb2y + dhay*vb2x)*rar*a;

  f4[0] = dahx*rhr*a;
  f4[1] = dahy*rhr*a;
  f4[2] = dahz*rhr*a;

  f1[0] = -(f2[0] + f3[0] + f4[0]);
  f1[1] = -(f2[1] + f3[1] + f4[1]);
  f1[2] = -(f2[2] + f3[2] + f4[2]);

  // note: i2 gets f3, i3 gets f2 (swapped relative to local variable names)
  if (NEWTON_BOND || i1 < nlocal) {
    f(i1,0) += static_cast<KK_ACC_FLOAT>(f1[0]);
    f(i1,1) += static_cast<KK_ACC_FLOAT>(f1[1]);
    f(i1,2) += static_cast<KK_ACC_FLOAT>(f1[2]);
  }
  if (NEWTON_BOND || i2 < nlocal) {
    f(i2,0) += static_cast<KK_ACC_FLOAT>(f3[0]);
    f(i2,1) += static_cast<KK_ACC_FLOAT>(f3[1]);
    f(i2,2) += static_cast<KK_ACC_FLOAT>(f3[2]);
  }
  if (NEWTON_BOND || i3 < nlocal) {
    f(i3,0) += static_cast<KK_ACC_FLOAT>(f2[0]);
    f(i3,1) += static_cast<KK_ACC_FLOAT>(f2[1]);
    f(i3,2) += static_cast<KK_ACC_FLOAT>(f2[2]);
  }
  if (NEWTON_BOND || i4 < nlocal) {
    f(i4,0) += static_cast<KK_ACC_FLOAT>(f4[0]);
    f(i4,1) += static_cast<KK_ACC_FLOAT>(f4[1]);
    f(i4,2) += static_cast<KK_ACC_FLOAT>(f4[2]);
  }

  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,i4,eimproper,f1,f2,f4,
             -vb1x,-vb1y,-vb1z,vb2x-vb1x,vb2y-vb1y,vb2z-vb1z,vb3x-vb2x,vb3y-vb2y,vb3z-vb2z);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperFourierKokkos<DeviceType>::allocate()
{
  ImproperFourier::allocate();
  int n = atom->nimpropertypes;
  k_k   = DAT::tdual_kkfloat_1d("ImproperFourier::k",n+1);
  k_C0  = DAT::tdual_kkfloat_1d("ImproperFourier::C0",n+1);
  k_C1  = DAT::tdual_kkfloat_1d("ImproperFourier::C1",n+1);
  k_C2  = DAT::tdual_kkfloat_1d("ImproperFourier::C2",n+1);
  k_all = DAT::tdual_int_1d("ImproperFourier::all",n+1);
  d_k   = k_k.template view<DeviceType>();
  d_C0  = k_C0.template view<DeviceType>();
  d_C1  = k_C1.template view<DeviceType>();
  d_C2  = k_C2.template view<DeviceType>();
  d_all = k_all.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperFourierKokkos<DeviceType>::coeff(int narg, char **arg)
{
  ImproperFourier::coeff(narg, arg);
  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nimpropertypes,ilo,ihi,error);
  for (int i = ilo; i <= ihi; i++) {
    k_k.view_host()[i]   = static_cast<KK_FLOAT>(k[i]);
    k_C0.view_host()[i]  = static_cast<KK_FLOAT>(C0[i]);
    k_C1.view_host()[i]  = static_cast<KK_FLOAT>(C1[i]);
    k_C2.view_host()[i]  = static_cast<KK_FLOAT>(C2[i]);
    k_all.view_host()[i] = all[i];
  }
  k_k.modify_host();
  k_C0.modify_host();
  k_C1.modify_host();
  k_C2.modify_host();
  k_all.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperFourierKokkos<DeviceType>::read_restart(FILE *fp)
{
  ImproperFourier::read_restart(fp);
  int n = atom->nimpropertypes;
  for (int i = 1; i <= n; i++) {
    k_k.view_host()[i]   = static_cast<KK_FLOAT>(k[i]);
    k_C0.view_host()[i]  = static_cast<KK_FLOAT>(C0[i]);
    k_C1.view_host()[i]  = static_cast<KK_FLOAT>(C1[i]);
    k_C2.view_host()[i]  = static_cast<KK_FLOAT>(C2[i]);
    k_all.view_host()[i] = all[i];
  }
  k_k.modify_host();
  k_C0.modify_host();
  k_C1.modify_host();
  k_C2.modify_host();
  k_all.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperFourierKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2,
    const int i3, const int i4, KK_FLOAT &eimproper, KK_FLOAT *f1, KK_FLOAT *f3, KK_FLOAT *f4,
    const KK_FLOAT &vb1x, const KK_FLOAT &vb1y, const KK_FLOAT &vb1z,
    const KK_FLOAT &vb2x, const KK_FLOAT &vb2y, const KK_FLOAT &vb2z,
    const KK_FLOAT &vb3x, const KK_FLOAT &vb3y, const KK_FLOAT &vb3z) const
{
  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(eimproper);
      else {
        KK_ACC_FLOAT eq = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*eimproper);
        if (i1 < nlocal) ev.evdwl += eq;
        if (i2 < nlocal) ev.evdwl += eq;
        if (i3 < nlocal) ev.evdwl += eq;
        if (i4 < nlocal) ev.evdwl += eq;
      }
    }
    if (eflag_atom) {
      KK_ACC_FLOAT eq = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*eimproper);
      if (newton_bond || i1 < nlocal) d_eatom[i1] += eq;
      if (newton_bond || i2 < nlocal) d_eatom[i2] += eq;
      if (newton_bond || i3 < nlocal) d_eatom[i3] += eq;
      if (newton_bond || i4 < nlocal) d_eatom[i4] += eq;
    }
  }
  if (vflag_either) {
    KK_ACC_FLOAT v[6];
    v[0] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1x*f1[0] + vb2x*f3[0] + (vb3x+vb2x)*f4[0]));
    v[1] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1y*f1[1] + vb2y*f3[1] + (vb3y+vb2y)*f4[1]));
    v[2] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1z*f1[2] + vb2z*f3[2] + (vb3z+vb2z)*f4[2]));
    v[3] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1x*f1[1] + vb2x*f3[1] + (vb3x+vb2x)*f4[1]));
    v[4] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1x*f1[2] + vb2x*f3[2] + (vb3x+vb2x)*f4[2]));
    v[5] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.25)*(vb1y*f1[2] + vb2y*f3[2] + (vb3y+vb2y)*f4[2]));
    if (vflag_global) {
      if (newton_bond) {
        for (int m = 0; m < 6; m++) ev.v[m] += static_cast<KK_ACC_FLOAT>(4.0)*v[m];
      } else {
        if (i1 < nlocal) for (int m = 0; m < 6; m++) ev.v[m] += v[m];
        if (i2 < nlocal) for (int m = 0; m < 6; m++) ev.v[m] += v[m];
        if (i3 < nlocal) for (int m = 0; m < 6; m++) ev.v[m] += v[m];
        if (i4 < nlocal) for (int m = 0; m < 6; m++) ev.v[m] += v[m];
      }
    }
    if (vflag_atom) {
      if (newton_bond || i1 < nlocal) for (int m = 0; m < 6; m++) d_vatom(i1,m) += v[m];
      if (newton_bond || i2 < nlocal) for (int m = 0; m < 6; m++) d_vatom(i2,m) += v[m];
      if (newton_bond || i3 < nlocal) for (int m = 0; m < 6; m++) d_vatom(i3,m) += v[m];
      if (newton_bond || i4 < nlocal) for (int m = 0; m < 6; m++) d_vatom(i4,m) += v[m];
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ImproperFourierKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ImproperFourierKokkos<LMPHostType>;
#endif
}
