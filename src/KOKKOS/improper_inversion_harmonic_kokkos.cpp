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

#include "improper_inversion_harmonic_kokkos.h"
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
ImproperInversionHarmonicKokkos<DeviceType>::ImproperInversionHarmonicKokkos(LAMMPS *lmp) : ImproperInversionHarmonic(lmp)
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
ImproperInversionHarmonicKokkos<DeviceType>::~ImproperInversionHarmonicKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperInversionHarmonicKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  k_kw.template sync<DeviceType>();
  k_w0.template sync<DeviceType>();
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
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperInversionHarmonicCompute<1,1> >(0,nimproperlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagImproperInversionHarmonicCompute<0,1> >(0,nimproperlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperInversionHarmonicCompute<1,0> >(0,nimproperlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagImproperInversionHarmonicCompute<0,0> >(0,nimproperlist),*this);
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
void ImproperInversionHarmonicKokkos<DeviceType>::operator()(TagImproperInversionHarmonicCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const
{
  const int i1 = improperlist(n,0);
  const int i2 = improperlist(n,1);
  const int i3 = improperlist(n,2);
  const int i4 = improperlist(n,3);
  const int type = improperlist(n,4);

  // 1st bond: IJ
  const KK_FLOAT vb1x = x(i2,0) - x(i1,0);
  const KK_FLOAT vb1y = x(i2,1) - x(i1,1);
  const KK_FLOAT vb1z = x(i2,2) - x(i1,2);
  const KK_FLOAT rrvb1 = static_cast<KK_FLOAT>(1.0)/sqrt(vb1x*vb1x + vb1y*vb1y + vb1z*vb1z);
  const KK_FLOAT rr2vb1 = rrvb1*rrvb1;

  // 2nd bond: IK
  const KK_FLOAT vb2x = x(i3,0) - x(i1,0);
  const KK_FLOAT vb2y = x(i3,1) - x(i1,1);
  const KK_FLOAT vb2z = x(i3,2) - x(i1,2);
  const KK_FLOAT rrvb2 = static_cast<KK_FLOAT>(1.0)/sqrt(vb2x*vb2x + vb2y*vb2y + vb2z*vb2z);
  const KK_FLOAT rr2vb2 = rrvb2*rrvb2;

  // 3rd bond: IL
  const KK_FLOAT vb3x = x(i4,0) - x(i1,0);
  const KK_FLOAT vb3y = x(i4,1) - x(i1,1);
  const KK_FLOAT vb3z = x(i4,2) - x(i1,2);
  const KK_FLOAT rrvb3 = static_cast<KK_FLOAT>(1.0)/sqrt(vb3x*vb3x + vb3y*vb3y + vb3z*vb3z);
  const KK_FLOAT rr2vb3 = rrvb3*rrvb3;

  invang<NEWTON_BOND,EVFLAG>(ev, i1,i2,i3,i4, type,
      vb3x,vb3y,vb3z, rrvb3,rr2vb3,
      vb2x,vb2y,vb2z, rrvb2,rr2vb2,
      vb1x,vb1y,vb1z, rrvb1,rr2vb1);
  invang<NEWTON_BOND,EVFLAG>(ev, i1,i3,i4,i2, type,
      vb1x,vb1y,vb1z, rrvb1,rr2vb1,
      vb3x,vb3y,vb3z, rrvb3,rr2vb3,
      vb2x,vb2y,vb2z, rrvb2,rr2vb2);
  invang<NEWTON_BOND,EVFLAG>(ev, i1,i4,i2,i3, type,
      vb2x,vb2y,vb2z, rrvb2,rr2vb2,
      vb1x,vb1y,vb1z, rrvb1,rr2vb1,
      vb3x,vb3y,vb3z, rrvb3,rr2vb3);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperInversionHarmonicKokkos<DeviceType>::operator()(TagImproperInversionHarmonicCompute<NEWTON_BOND,EVFLAG>, const int &n) const
{
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagImproperInversionHarmonicCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperInversionHarmonicKokkos<DeviceType>::invang(EV_FLOAT &ev,
    const int i1, const int i2, const int i3, const int i4,
    const int type,
    const KK_FLOAT vb1x, const KK_FLOAT vb1y, const KK_FLOAT vb1z,
    const KK_FLOAT rrvb1, const KK_FLOAT rr2vb1,
    const KK_FLOAT vb2x, const KK_FLOAT vb2y, const KK_FLOAT vb2z,
    const KK_FLOAT rrvb2, const KK_FLOAT rr2vb2,
    const KK_FLOAT vb3x, const KK_FLOAT vb3y, const KK_FLOAT vb3z,
    const KK_FLOAT rrvb3, const KK_FLOAT rr2vb3) const
{
  KK_FLOAT eimproper,f1[3],f2[3],f3[3],f4[3];
  KK_FLOAT omega,cosomega,domega,gomega,rjk,rjl;
  KK_FLOAT upx,upy,upz,upn,rup,umx,umy,umz,umn,rum,wwr;
  KK_FLOAT rucb,rudb,rvcb,rvdb,rupupn,rumumn;

  eimproper = static_cast<KK_FLOAT>(0.0);

  rjk = vb3x*vb2x + vb3y*vb2y + vb3z*vb2z;
  rjl = vb1x*vb3x + vb1y*vb3y + vb1z*vb3z;

  upx = vb2x*rrvb2 + vb1x*rrvb1;
  upy = vb2y*rrvb2 + vb1y*rrvb1;
  upz = vb2z*rrvb2 + vb1z*rrvb1;
  upn = static_cast<KK_FLOAT>(1.0)/sqrt(upx*upx + upy*upy + upz*upz);
  upx *= upn; upy *= upn; upz *= upn;
  rup = vb3x*upx + vb3y*upy + vb3z*upz;

  umx = vb2x*rrvb2 - vb1x*rrvb1;
  umy = vb2y*rrvb2 - vb1y*rrvb1;
  umz = vb2z*rrvb2 - vb1z*rrvb1;
  umn = static_cast<KK_FLOAT>(1.0)/sqrt(umx*umx + umy*umy + umz*umz);
  umx *= umn; umy *= umn; umz *= umn;
  rum = vb3x*umx + vb3y*umy + vb3z*umz;

  wwr = sqrt(rup*rup + rum*rum);
  cosomega = wwr*rrvb3;
  if (cosomega > static_cast<KK_FLOAT>(1.0)) cosomega = static_cast<KK_FLOAT>(1.0);
  omega = acos(cosomega);

  domega = omega - d_w0[type];
  if (EVFLAG && eflag) eimproper = d_kw[type]*(domega*domega);

  gomega = static_cast<KK_FLOAT>(0.0);
  if (omega*omega > static_cast<KK_FLOAT>(1.0e-24))
    gomega = static_cast<KK_FLOAT>(2.0)*d_kw[type]*domega/sin(omega);

  rucb = rjk - rup*(vb2x*upx + vb2y*upy + vb2z*upz);
  rudb = rjl - rup*(vb1x*upx + vb1y*upy + vb1z*upz);
  rvcb = rjk - rum*(vb2x*umx + vb2y*umy + vb2z*umz);
  rvdb = rjl - rum*(vb1x*umx + vb1y*umy + vb1z*umz);

  rupupn = rup*upn;
  rumumn = rum*umn;

  f2[0] = gomega*(-cosomega*vb3x*rr2vb3 + rrvb3*(rup*upx + rum*umx)/wwr);
  f2[1] = gomega*(-cosomega*vb3y*rr2vb3 + rrvb3*(rup*upy + rum*umy)/wwr);
  f2[2] = gomega*(-cosomega*vb3z*rr2vb3 + rrvb3*(rup*upz + rum*umz)/wwr);

  f3[0] = gomega*rrvb3*(rupupn*rrvb2*(vb3x - rup*upx - rucb*vb2x*rr2vb2) +
          rumumn*rrvb2*(vb3x - rum*umx - rvcb*vb2x*rr2vb2))/wwr;
  f3[1] = gomega*rrvb3*(rupupn*rrvb2*(vb3y - rup*upy - rucb*vb2y*rr2vb2) +
          rumumn*rrvb2*(vb3y - rum*umy - rvcb*vb2y*rr2vb2))/wwr;
  f3[2] = gomega*rrvb3*(rupupn*rrvb2*(vb3z - rup*upz - rucb*vb2z*rr2vb2) +
          rumumn*rrvb2*(vb3z - rum*umz - rvcb*vb2z*rr2vb2))/wwr;

  f4[0] = gomega*rrvb3*(rupupn*rrvb1*(vb3x - rup*upx - rudb*vb1x*rr2vb1) -
          rumumn*rrvb1*(vb3x - rum*umx - rvdb*vb1x*rr2vb1))/wwr;
  f4[1] = gomega*rrvb3*(rupupn*rrvb1*(vb3y - rup*upy - rudb*vb1y*rr2vb1) -
          rumumn*rrvb1*(vb3y - rum*umy - rvdb*vb1y*rr2vb1))/wwr;
  f4[2] = gomega*rrvb3*(rupupn*rrvb1*(vb3z - rup*upz - rudb*vb1z*rr2vb1) -
          rumumn*rrvb1*(vb3z - rum*umz - rvdb*vb1z*rr2vb1))/wwr;

  f1[0] = -(f2[0] + f3[0] + f4[0]);
  f1[1] = -(f2[1] + f3[1] + f4[1]);
  f1[2] = -(f2[2] + f3[2] + f4[2]);

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

  if (EVFLAG) {
    const KK_FLOAT rb3x = vb1x - vb2x;
    const KK_FLOAT rb3y = vb1y - vb2y;
    const KK_FLOAT rb3z = vb1z - vb2z;
    ev_tally(ev,i1,i2,i3,i4,eimproper,f2,f3,f4,
             vb3x,vb3y,vb3z, vb2x,vb2y,vb2z, rb3x,rb3y,rb3z);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperInversionHarmonicKokkos<DeviceType>::allocate()
{
  ImproperInversionHarmonic::allocate();
  int n = atom->nimpropertypes;
  k_kw = DAT::tdual_kkfloat_1d("ImproperInversionHarmonic::kw",n+1);
  k_w0 = DAT::tdual_kkfloat_1d("ImproperInversionHarmonic::w0",n+1);
  d_kw = k_kw.template view<DeviceType>();
  d_w0 = k_w0.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperInversionHarmonicKokkos<DeviceType>::coeff(int narg, char **arg)
{
  ImproperInversionHarmonic::coeff(narg, arg);
  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nimpropertypes,ilo,ihi,error);
  for (int i = ilo; i <= ihi; i++) {
    k_kw.view_host()[i] = static_cast<KK_FLOAT>(kw[i]);
    k_w0.view_host()[i] = static_cast<KK_FLOAT>(w0[i]);
  }
  k_kw.modify_host();
  k_w0.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ImproperInversionHarmonicKokkos<DeviceType>::read_restart(FILE *fp)
{
  ImproperInversionHarmonic::read_restart(fp);
  int n = atom->nimpropertypes;
  for (int i = 1; i <= n; i++) {
    k_kw.view_host()[i] = static_cast<KK_FLOAT>(kw[i]);
    k_w0.view_host()[i] = static_cast<KK_FLOAT>(w0[i]);
  }
  k_kw.modify_host();
  k_w0.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ImproperInversionHarmonicKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i1, const int i2,
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
template class ImproperInversionHarmonicKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ImproperInversionHarmonicKokkos<LMPHostType>;
#endif
}
