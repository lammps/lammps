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

#include "angle_mm3_kokkos.h"

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

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleMM3Kokkos<DeviceType>::AngleMM3Kokkos(LAMMPS *lmp) : AngleMM3(lmp)
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
AngleMM3Kokkos<DeviceType>::~AngleMM3Kokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleMM3Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"angle:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"angle:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  k_theta0.template sync<DeviceType>();
  k_k2.template sync<DeviceType>();

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_anglelist.template sync<DeviceType>();
  anglelist = neighborKK->k_anglelist.template view<DeviceType>();
  int nanglelist = neighborKK->nanglelist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
        TagAngleMM3Compute<1,1> >(0,nanglelist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
        TagAngleMM3Compute<0,1> >(0,nanglelist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
        TagAngleMM3Compute<1,0> >(0,nanglelist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
        TagAngleMM3Compute<0,0> >(0,nanglelist),*this);
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

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleMM3Kokkos<DeviceType>::operator()(TagAngleMM3Compute<NEWTON_BOND,EVFLAG>,
                                            const int &n, EV_FLOAT& ev) const
{
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,
    typename KKDevice<DeviceType>::value,
    Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

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

  // cos and sin of angle

  KK_FLOAT c = (delx1*delx2 + dely1*dely2 + delz1*delz2) / (r1*r2);
  if (c > static_cast<KK_FLOAT>(1.0))  c = static_cast<KK_FLOAT>(1.0);
  if (c < static_cast<KK_FLOAT>(-1.0)) c = static_cast<KK_FLOAT>(-1.0);

  KK_FLOAT s = sqrt(static_cast<KK_FLOAT>(1.0) - c*c);
  if (s < static_cast<KK_FLOAT>(SMALL)) s = static_cast<KK_FLOAT>(SMALL);
  s = static_cast<KK_FLOAT>(1.0)/s;

  // force & energy for MM3 angle term (dtheta in radians)

  const KK_FLOAT dtheta  = acos(c) - d_theta0[type];
  const KK_FLOAT dtheta2 = dtheta*dtheta;
  const KK_FLOAT dtheta3 = dtheta2*dtheta;
  const KK_FLOAT dtheta4 = dtheta3*dtheta;

  const KK_FLOAT de_angle = static_cast<KK_FLOAT>(2.0)*d_k2[type]*dtheta*(
    static_cast<KK_FLOAT>(1.0)
    - static_cast<KK_FLOAT>(1.203211)*dtheta
    + static_cast<KK_FLOAT>(0.367674)*dtheta2
    - static_cast<KK_FLOAT>(0.3239159)*dtheta3
    + static_cast<KK_FLOAT>(0.711270)*dtheta4);

  const KK_FLOAT a   = -de_angle*s;
  const KK_FLOAT a11 = a*c/rsq1;
  const KK_FLOAT a12 = -a/(r1*r2);
  const KK_FLOAT a22 = a*c/rsq2;

  KK_FLOAT f1[3], f3[3];
  f1[0] = a11*delx1 + a12*delx2;
  f1[1] = a11*dely1 + a12*dely2;
  f1[2] = a11*delz1 + a12*delz2;
  f3[0] = a22*delx2 + a12*delx1;
  f3[1] = a22*dely2 + a12*dely1;
  f3[2] = a22*delz2 + a12*delz1;

  KK_FLOAT eangle = static_cast<KK_FLOAT>(0.0);
  if (EVFLAG && eflag)
    eangle = d_k2[type]*dtheta2*(
      static_cast<KK_FLOAT>(1.0)
      - static_cast<KK_FLOAT>(0.802141)*dtheta
      + static_cast<KK_FLOAT>(0.183837)*dtheta2
      - static_cast<KK_FLOAT>(0.131664)*dtheta3
      + static_cast<KK_FLOAT>(0.237090)*dtheta4);

  // apply force to each of 3 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    a_f(i1,0) += static_cast<KK_ACC_FLOAT>(f1[0]);
    a_f(i1,1) += static_cast<KK_ACC_FLOAT>(f1[1]);
    a_f(i1,2) += static_cast<KK_ACC_FLOAT>(f1[2]);
  }
  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) -= static_cast<KK_ACC_FLOAT>(f1[0] + f3[0]);
    a_f(i2,1) -= static_cast<KK_ACC_FLOAT>(f1[1] + f3[1]);
    a_f(i2,2) -= static_cast<KK_ACC_FLOAT>(f1[2] + f3[2]);
  }
  if (NEWTON_BOND || i3 < nlocal) {
    a_f(i3,0) += static_cast<KK_ACC_FLOAT>(f3[0]);
    a_f(i3,1) += static_cast<KK_ACC_FLOAT>(f3[1]);
    a_f(i3,2) += static_cast<KK_ACC_FLOAT>(f3[2]);
  }

  if (EVFLAG)
    ev_tally(ev,i1,i2,i3,eangle,f1,f3,
             delx1,dely1,delz1,delx2,dely2,delz2);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleMM3Kokkos<DeviceType>::operator()(TagAngleMM3Compute<NEWTON_BOND,EVFLAG>,
                                            const int &n) const
{
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleMM3Compute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleMM3Kokkos<DeviceType>::allocate()
{
  AngleMM3::allocate();
  int n = atom->nangletypes;
  k_theta0 = DAT::tdual_kkfloat_1d("AngleMM3::theta0",n+1);
  k_k2     = DAT::tdual_kkfloat_1d("AngleMM3::k2",n+1);
  d_theta0 = k_theta0.template view<DeviceType>();
  d_k2     = k_k2.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleMM3Kokkos<DeviceType>::coeff(int narg, char **arg)
{
  AngleMM3::coeff(narg, arg);
  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);
  for (int i = ilo; i <= ihi; i++) {
    k_theta0.view_host()[i] = static_cast<KK_FLOAT>(theta0[i]);
    k_k2.view_host()[i]     = static_cast<KK_FLOAT>(k2[i]);
  }
  k_theta0.modify_host();
  k_k2.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleMM3Kokkos<DeviceType>::read_restart(FILE *fp)
{
  AngleMM3::read_restart(fp);
  int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_theta0.view_host()[i] = static_cast<KK_FLOAT>(theta0[i]);
    k_k2.view_host()[i]     = static_cast<KK_FLOAT>(k2[i]);
  }
  k_theta0.modify_host();
  k_k2.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleMM3Kokkos<DeviceType>::ev_tally(EV_FLOAT &ev,
    const int i, const int j, const int k,
    KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
    const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
    const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const
{
  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,
    typename KKDevice<DeviceType>::value,
    Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,
    typename KKDevice<DeviceType>::value,
    Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = d_vatom;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(eangle);
      else {
        const KK_ACC_FLOAT et = static_cast<KK_ACC_FLOAT>(THIRD*eangle);
        if (i < nlocal) ev.evdwl += et;
        if (j < nlocal) ev.evdwl += et;
        if (k < nlocal) ev.evdwl += et;
      }
    }
    if (eflag_atom) {
      const KK_ACC_FLOAT et = static_cast<KK_ACC_FLOAT>(THIRD*eangle);
      if (newton_bond || i < nlocal) v_eatom[i] += et;
      if (newton_bond || j < nlocal) v_eatom[j] += et;
      if (newton_bond || k < nlocal) v_eatom[k] += et;
    }
  }

  if (vflag_either) {
    KK_FLOAT v[6];
    v[0] = delx1*f1[0] + delx2*f3[0];
    v[1] = dely1*f1[1] + dely2*f3[1];
    v[2] = delz1*f1[2] + delz2*f3[2];
    v[3] = delx1*f1[1] + delx2*f3[1];
    v[4] = delx1*f1[2] + delx2*f3[2];
    v[5] = dely1*f1[2] + dely2*f3[2];

    if (vflag_global) {
      if (newton_bond) {
        for (int m = 0; m < 6; m++) ev.v[m] += static_cast<KK_ACC_FLOAT>(v[m]);
      } else {
        const KK_FLOAT vt = THIRD;
        if (i < nlocal) for (int m = 0; m < 6; m++) ev.v[m] += static_cast<KK_ACC_FLOAT>(vt*v[m]);
        if (j < nlocal) for (int m = 0; m < 6; m++) ev.v[m] += static_cast<KK_ACC_FLOAT>(vt*v[m]);
        if (k < nlocal) for (int m = 0; m < 6; m++) ev.v[m] += static_cast<KK_ACC_FLOAT>(vt*v[m]);
      }
    }
    if (vflag_atom) {
      const KK_FLOAT vt = THIRD;
      if (newton_bond || i < nlocal)
        for (int m = 0; m < 6; m++) v_vatom(i,m) += static_cast<KK_ACC_FLOAT>(vt*v[m]);
      if (newton_bond || j < nlocal)
        for (int m = 0; m < 6; m++) v_vatom(j,m) += static_cast<KK_ACC_FLOAT>(vt*v[m]);
      if (newton_bond || k < nlocal)
        for (int m = 0; m < 6; m++) v_vatom(k,m) += static_cast<KK_ACC_FLOAT>(vt*v[m]);
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class AngleMM3Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class AngleMM3Kokkos<LMPHostType>;
#endif
}
