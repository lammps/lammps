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

#include "improper_ring_kokkos.h"
#include <cmath>
#include "atom_kokkos.h"
#include "comm.h"
#include "neighbor_kokkos.h"
#include "force.h"
#include "memory_kokkos.h"
#include "error.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ImproperRingKokkos<DeviceType>::ImproperRingKokkos(LAMMPS *lmp) : ImproperRing(lmp)
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
ImproperRingKokkos<DeviceType>::~ImproperRingKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperRingKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperRingCompute<1,1> >(0,nimproperlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperRingCompute<0,1> >(0,nimproperlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperRingCompute<1,0> >(0,nimproperlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperRingCompute<0,0> >(0,nimproperlist),*this);
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
void ImproperRingKokkos<DeviceType>::operator()(TagImproperRingCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const
{
  const int i1 = improperlist(n,0);
  const int i2 = improperlist(n,1);
  const int i3 = improperlist(n,2);
  const int i4 = improperlist(n,3);
  const int type = improperlist(n,4);

  // vb1,vb2,vb3 for ev_tally compatibility
  const KK_FLOAT vb1x = x(i1,0) - x(i2,0);
  const KK_FLOAT vb1y = x(i1,1) - x(i2,1);
  const KK_FLOAT vb1z = x(i1,2) - x(i2,2);

  const KK_FLOAT vb2x = x(i3,0) - x(i2,0);
  const KK_FLOAT vb2y = x(i3,1) - x(i2,1);
  const KK_FLOAT vb2z = x(i3,2) - x(i2,2);

  const KK_FLOAT vb3x = x(i4,0) - x(i3,0);
  const KK_FLOAT vb3y = x(i4,1) - x(i3,1);
  const KK_FLOAT vb3z = x(i4,2) - x(i3,2);

  // three angle combinations: at1-at2-at3
  // icomb 0: i1-i2-i4
  // icomb 1: i1-i2-i3
  // icomb 2: i4-i2-i3
  const int at1[3] = {i1, i1, i4};
  const int at2[3] = {i2, i2, i2};
  const int at3[3] = {i4, i3, i3};

  KK_FLOAT bvec1x[3],bvec1y[3],bvec1z[3];
  KK_FLOAT bvec2x[3],bvec2y[3],bvec2z[3];
  KK_FLOAT bvec1n[3],bvec2n[3],bend_angle[3];

  KK_FLOAT angle_summer = static_cast<KK_FLOAT>(0.0);
  const KK_FLOAT chi_type = d_chi[type];

  for (int icomb = 0; icomb < 3; icomb++) {
    bvec1x[icomb] = x(at2[icomb],0) - x(at1[icomb],0);
    bvec1y[icomb] = x(at2[icomb],1) - x(at1[icomb],1);
    bvec1z[icomb] = x(at2[icomb],2) - x(at1[icomb],2);
    bvec1n[icomb] = sqrt(bvec1x[icomb]*bvec1x[icomb] + bvec1y[icomb]*bvec1y[icomb] + bvec1z[icomb]*bvec1z[icomb]);

    bvec2x[icomb] = x(at3[icomb],0) - x(at2[icomb],0);
    bvec2y[icomb] = x(at3[icomb],1) - x(at2[icomb],1);
    bvec2z[icomb] = x(at3[icomb],2) - x(at2[icomb],2);
    bvec2n[icomb] = sqrt(bvec2x[icomb]*bvec2x[icomb] + bvec2y[icomb]*bvec2y[icomb] + bvec2z[icomb]*bvec2z[icomb]);

    bend_angle[icomb] = (bvec2x[icomb]*bvec1x[icomb] + bvec2y[icomb]*bvec1y[icomb] + bvec2z[icomb]*bvec1z[icomb]);
    bend_angle[icomb] /= (bvec1n[icomb] * bvec2n[icomb]);
    if (bend_angle[icomb] >  static_cast<KK_FLOAT>(1.0)) bend_angle[icomb] -= static_cast<KK_FLOAT>(SMALL);
    if (bend_angle[icomb] < static_cast<KK_FLOAT>(-1.0)) bend_angle[icomb] += static_cast<KK_FLOAT>(SMALL);

    angle_summer += (bend_angle[icomb] - chi_type);
  }

  KK_FLOAT eimproper = static_cast<KK_FLOAT>(0.0);
  if (EVFLAG && eflag) {
    KK_FLOAT as2 = angle_summer*angle_summer;
    eimproper = (static_cast<KK_FLOAT>(1.0)/static_cast<KK_FLOAT>(6.0))*d_k[type]*as2*as2*as2;
  }

  const KK_FLOAT as4 = angle_summer*angle_summer*angle_summer*angle_summer;
  const KK_FLOAT angfac = d_k[type]*as4*angle_summer;

  KK_FLOAT f1[3],f3[3],f4[3];
  f1[0] = static_cast<KK_FLOAT>(0.0); f1[1] = static_cast<KK_FLOAT>(0.0); f1[2] = static_cast<KK_FLOAT>(0.0);
  f3[0] = static_cast<KK_FLOAT>(0.0); f3[1] = static_cast<KK_FLOAT>(0.0); f3[2] = static_cast<KK_FLOAT>(0.0);
  f4[0] = static_cast<KK_FLOAT>(0.0); f4[1] = static_cast<KK_FLOAT>(0.0); f4[2] = static_cast<KK_FLOAT>(0.0);

  for (int icomb = 0; icomb < 3; icomb++) {
    const KK_FLOAT cjiji = bvec1n[icomb]*bvec1n[icomb];
    const KK_FLOAT ckjkj = bvec2n[icomb]*bvec2n[icomb];
    const KK_FLOAT ckjji = bvec2x[icomb]*bvec1x[icomb] + bvec2y[icomb]*bvec1y[icomb] + bvec2z[icomb]*bvec1z[icomb];

    const KK_FLOAT cfact1 = angfac / sqrt(ckjkj * cjiji);
    const KK_FLOAT cfact2 = ckjji / ckjkj;
    const KK_FLOAT cfact3 = ckjji / cjiji;

    const KK_FLOAT fkx = cfact2*bvec2x[icomb] - bvec1x[icomb];
    const KK_FLOAT fky = cfact2*bvec2y[icomb] - bvec1y[icomb];
    const KK_FLOAT fkz = cfact2*bvec2z[icomb] - bvec1z[icomb];

    const KK_FLOAT fix = bvec2x[icomb] - cfact3*bvec1x[icomb];
    const KK_FLOAT fiy = bvec2y[icomb] - cfact3*bvec1y[icomb];
    const KK_FLOAT fiz = bvec2z[icomb] - cfact3*bvec1z[icomb];

    const KK_FLOAT fjx = -fix - fkx;
    const KK_FLOAT fjy = -fiy - fky;
    const KK_FLOAT fjz = -fiz - fkz;

    const KK_FLOAT fix_s = cfact1*fix;
    const KK_FLOAT fiy_s = cfact1*fiy;
    const KK_FLOAT fiz_s = cfact1*fiz;
    const KK_FLOAT fkx_s = cfact1*fkx;
    const KK_FLOAT fky_s = cfact1*fky;
    const KK_FLOAT fkz_s = cfact1*fkz;

    if (at1[icomb] == i1)       { f1[0] += cfact1*fix; f1[1] += cfact1*fiy; f1[2] += cfact1*fiz; }
    else if (at2[icomb] == i1)  { f1[0] += cfact1*fjx; f1[1] += cfact1*fjy; f1[2] += cfact1*fjz; }
    else if (at3[icomb] == i1)  { f1[0] += cfact1*fkx; f1[1] += cfact1*fky; f1[2] += cfact1*fkz; }

    if (at1[icomb] == i3)       { f3[0] += cfact1*fix; f3[1] += cfact1*fiy; f3[2] += cfact1*fiz; }
    else if (at2[icomb] == i3)  { f3[0] += cfact1*fjx; f3[1] += cfact1*fjy; f3[2] += cfact1*fjz; }
    else if (at3[icomb] == i3)  { f3[0] += cfact1*fkx; f3[1] += cfact1*fky; f3[2] += cfact1*fkz; }

    if (at1[icomb] == i4)       { f4[0] += cfact1*fix; f4[1] += cfact1*fiy; f4[2] += cfact1*fiz; }
    else if (at2[icomb] == i4)  { f4[0] += cfact1*fjx; f4[1] += cfact1*fjy; f4[2] += cfact1*fjz; }
    else if (at3[icomb] == i4)  { f4[0] += cfact1*fkx; f4[1] += cfact1*fky; f4[2] += cfact1*fkz; }

    // also apply to i2 (atom index at2[icomb] == i2 always)
    if (NEWTON_BOND || i2 < nlocal) {
      f(i2,0) += static_cast<KK_ACC_FLOAT>(cfact1*fjx);
      f(i2,1) += static_cast<KK_ACC_FLOAT>(cfact1*fjy);
      f(i2,2) += static_cast<KK_ACC_FLOAT>(cfact1*fjz);
    }
    if (NEWTON_BOND || at1[icomb] < nlocal) {
      f(at1[icomb],0) += static_cast<KK_ACC_FLOAT>(fix_s);
      f(at1[icomb],1) += static_cast<KK_ACC_FLOAT>(fiy_s);
      f(at1[icomb],2) += static_cast<KK_ACC_FLOAT>(fiz_s);
    }
    if (NEWTON_BOND || at3[icomb] < nlocal) {
      f(at3[icomb],0) += static_cast<KK_ACC_FLOAT>(fkx_s);
      f(at3[icomb],1) += static_cast<KK_ACC_FLOAT>(fky_s);
      f(at3[icomb],2) += static_cast<KK_ACC_FLOAT>(fkz_s);
    }
  }

  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,i4,eimproper,f1,f3,f4,
             vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperRingKokkos<DeviceType>::operator()(TagImproperRingCompute<NEWTON_BOND,EVFLAG>, const int &n) const
{
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagImproperRingCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperRingKokkos<DeviceType>::allocate()
{
  ImproperRing::allocate();
  int n = atom->nimpropertypes;
  k_k = DAT::tdual_kkfloat_1d("ImproperRing::k",n+1);
  k_chi = DAT::tdual_kkfloat_1d("ImproperRing::chi",n+1);
  d_k = k_k.template view<DeviceType>();
  d_chi = k_chi.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperRingKokkos<DeviceType>::coeff(int narg, char **arg)
{
  ImproperRing::coeff(narg, arg);
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
void ImproperRingKokkos<DeviceType>::read_restart(FILE *fp)
{
  ImproperRing::read_restart(fp);
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
void ImproperRingKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2,
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
template class ImproperRingKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ImproperRingKokkos<LMPHostType>;
#endif
}
