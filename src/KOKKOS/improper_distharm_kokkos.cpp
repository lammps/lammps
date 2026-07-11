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

#include "improper_distharm_kokkos.h"
#include <cmath>
#include "atom_kokkos.h"
#include "comm.h"
#include "neighbor_kokkos.h"
#include "force.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ImproperDistHarmKokkos<DeviceType>::ImproperDistHarmKokkos(LAMMPS *lmp) : ImproperDistHarm(lmp)
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
ImproperDistHarmKokkos<DeviceType>::~ImproperDistHarmKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperDistHarmKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_chi.template sync<DeviceType>();
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

  copymode = 1;

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperDistHarmCompute<1,1> >(0,nimproperlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperDistHarmCompute<0,1> >(0,nimproperlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperDistHarmCompute<1,0> >(0,nimproperlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperDistHarmCompute<0,0> >(0,nimproperlist),*this);
    }
  }

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
void ImproperDistHarmKokkos<DeviceType>::operator()(TagImproperDistHarmCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const
{
  const int i1 = improperlist(n,0);
  const int i2 = improperlist(n,1);
  const int i3 = improperlist(n,2);
  const int i4 = improperlist(n,3);
  const int type = improperlist(n,4);

  // i4 is the central atom; bonds from i1-i2, i1-i3, i1-i4, i2-i3, i2-i4, i3-i4

  const KK_FLOAT xab = x(i2,0) - x(i1,0);
  const KK_FLOAT yab = x(i2,1) - x(i1,1);
  const KK_FLOAT zab = x(i2,2) - x(i1,2);

  const KK_FLOAT xac = x(i3,0) - x(i1,0);
  const KK_FLOAT yac = x(i3,1) - x(i1,1);
  const KK_FLOAT zac = x(i3,2) - x(i1,2);

  const KK_FLOAT xad = x(i4,0) - x(i1,0);
  const KK_FLOAT yad = x(i4,1) - x(i1,1);
  const KK_FLOAT zad = x(i4,2) - x(i1,2);

  const KK_FLOAT xbc = x(i3,0) - x(i2,0);
  const KK_FLOAT ybc = x(i3,1) - x(i2,1);
  const KK_FLOAT zbc = x(i3,2) - x(i2,2);

  const KK_FLOAT xcd = x(i4,0) - x(i3,0);
  const KK_FLOAT ycd = x(i4,1) - x(i3,1);
  const KK_FLOAT zcd = x(i4,2) - x(i3,2);

  KK_FLOAT xna =   ybc*zcd - zbc*ycd;
  KK_FLOAT yna = -(xbc*zcd - zbc*xcd);
  KK_FLOAT zna =   xbc*ycd - ybc*xcd;
  const KK_FLOAT rna = static_cast<KK_FLOAT>(1.0) / sqrt(xna*xna + yna*yna + zna*zna);
  xna *= rna;
  yna *= rna;
  zna *= rna;

  const KK_FLOAT da = -(xna*xad + yna*yad + zna*zad);

  const KK_FLOAT domega = da - d_chi[type];

  KK_FLOAT eimproper = static_cast<KK_FLOAT>(0.0);
  if (EVFLAG && eflag)
    eimproper = d_k[type]*domega*domega;

  const KK_FLOAT a = static_cast<KK_FLOAT>(2.0)*d_k[type]*domega;

  KK_FLOAT f1[3],f2[3],f3[3],f4[3];
  f1[0] = -a*xna;
  f1[1] = -a*yna;
  f1[2] = -a*zna;
  f4[0] = a*xna;
  f4[1] = a*yna;
  f4[2] = a*zna;

  f2[0] =  a*(yad*zcd - zad*ycd)*rna + a*da*rna*(yna*zcd - zna*ycd);
  f2[1] =  a*(zad*xcd - xad*zcd)*rna + a*da*rna*(zna*xcd - xna*zcd);
  f2[2] =  a*(xad*ycd - yad*xcd)*rna + a*da*rna*(xna*ycd - yna*xcd);

  f3[0] = -a*(yad*zcd - zad*ycd)*rna - a*da*rna*(yna*zcd - zna*ycd);
  f3[1] = -a*(zad*xcd - xad*zcd)*rna - a*da*rna*(zna*xcd - xna*zcd);
  f3[2] = -a*(xad*ycd - yad*xcd)*rna - a*da*rna*(xna*ycd - yna*xcd);

  f3[0] += -a*(yad*zbc - zad*ybc)*rna - a*da*rna*(yna*zbc - zna*ybc);
  f3[1] += -a*(zad*xbc - xad*zbc)*rna - a*da*rna*(zna*xbc - xna*zbc);
  f3[2] += -a*(xad*ybc - yad*xbc)*rna - a*da*rna*(xna*ybc - yna*xbc);
  f4[0] +=  a*(yad*zbc - zad*ybc)*rna + a*da*rna*(yna*zbc - zna*ybc);
  f4[1] +=  a*(zad*xbc - xad*zbc)*rna + a*da*rna*(zna*xbc - xna*zbc);
  f4[2] +=  a*(xad*ybc - yad*xbc)*rna + a*da*rna*(xna*ybc - yna*xbc);

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
    ev_tally(ev,i1,i2,i3,i4,eimproper,f2,f3,f4,
             xab,yab,zab,xac,yac,zac,xad-xac,yad-yac,zad-zac);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperDistHarmKokkos<DeviceType>::operator()(TagImproperDistHarmCompute<NEWTON_BOND,EVFLAG>, const int &n) const
{
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagImproperDistHarmCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperDistHarmKokkos<DeviceType>::allocate()
{
  ImproperDistHarm::allocate();
  int n = atom->nimpropertypes;
  k_k = DAT::tdual_kkfloat_1d("ImproperDistHarm::k",n+1);
  k_chi = DAT::tdual_kkfloat_1d("ImproperDistHarm::chi",n+1);
  d_k = k_k.template view<DeviceType>();
  d_chi = k_chi.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperDistHarmKokkos<DeviceType>::coeff(int narg, char **arg)
{
  ImproperDistHarm::coeff(narg, arg);
  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nimpropertypes,ilo,ihi,error);
  for (int i = ilo; i <= ihi; i++) {
    k_k.view_host()[i] = static_cast<KK_FLOAT>(k[i]);
    k_chi.view_host()[i] = static_cast<KK_FLOAT>(chi[i]);
  }
  k_k.modify_host();
  k_chi.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperDistHarmKokkos<DeviceType>::read_restart(FILE *fp)
{
  ImproperDistHarm::read_restart(fp);
  int n = atom->nimpropertypes;
  for (int i = 1; i <= n; i++) {
    k_k.view_host()[i] = static_cast<KK_FLOAT>(k[i]);
    k_chi.view_host()[i] = static_cast<KK_FLOAT>(chi[i]);
  }
  k_k.modify_host();
  k_chi.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperDistHarmKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2,
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
template class ImproperDistHarmKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ImproperDistHarmKokkos<LMPHostType>;
#endif
}
