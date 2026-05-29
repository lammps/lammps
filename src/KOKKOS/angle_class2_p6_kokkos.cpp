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

#include "angle_class2_p6_kokkos.h"

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
AngleClass2P6Kokkos<DeviceType>::AngleClass2P6Kokkos(LAMMPS *lmp) : AngleClass2P6(lmp)
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
AngleClass2P6Kokkos<DeviceType>::~AngleClass2P6Kokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleClass2P6Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  k_k3.template sync<DeviceType>();
  k_k4.template sync<DeviceType>();
  k_k5.template sync<DeviceType>();
  k_k6.template sync<DeviceType>();
  k_bb_k.template sync<DeviceType>();
  k_bb_r1.template sync<DeviceType>();
  k_bb_r2.template sync<DeviceType>();
  k_ba_k1.template sync<DeviceType>();
  k_ba_k2.template sync<DeviceType>();
  k_ba_r1.template sync<DeviceType>();
  k_ba_r2.template sync<DeviceType>();

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
        TagAngleClass2P6Compute<1,1> >(0,nanglelist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType,
        TagAngleClass2P6Compute<0,1> >(0,nanglelist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
        TagAngleClass2P6Compute<1,0> >(0,nanglelist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,
        TagAngleClass2P6Compute<0,0> >(0,nanglelist),*this);
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
void AngleClass2P6Kokkos<DeviceType>::operator()(TagAngleClass2P6Compute<NEWTON_BOND,EVFLAG>,
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

  // angle term

  const KK_FLOAT dtheta  = acos(c) - d_theta0[type];
  const KK_FLOAT dtheta2 = dtheta*dtheta;
  const KK_FLOAT dtheta3 = dtheta2*dtheta;
  const KK_FLOAT dtheta4 = dtheta3*dtheta;
  const KK_FLOAT dtheta5 = dtheta4*dtheta;
  const KK_FLOAT dtheta6 = dtheta5*dtheta;

  const KK_FLOAT de_angle =
    static_cast<KK_FLOAT>(2.0)*d_k2[type]*dtheta
    + static_cast<KK_FLOAT>(3.0)*d_k3[type]*dtheta2
    + static_cast<KK_FLOAT>(4.0)*d_k4[type]*dtheta3
    + static_cast<KK_FLOAT>(5.0)*d_k5[type]*dtheta4
    + static_cast<KK_FLOAT>(6.0)*d_k6[type]*dtheta5;

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
    eangle = d_k2[type]*dtheta2 + d_k3[type]*dtheta3
           + d_k4[type]*dtheta4 + d_k5[type]*dtheta5 + d_k6[type]*dtheta6;

  // bond-bond term

  const KK_FLOAT dr1_bb = r1 - d_bb_r1[type];
  const KK_FLOAT dr2_bb = r2 - d_bb_r2[type];
  const KK_FLOAT tk1 = d_bb_k[type]*dr1_bb;
  const KK_FLOAT tk2 = d_bb_k[type]*dr2_bb;

  f1[0] -= delx1*tk2/r1;
  f1[1] -= dely1*tk2/r1;
  f1[2] -= delz1*tk2/r1;
  f3[0] -= delx2*tk1/r2;
  f3[1] -= dely2*tk1/r2;
  f3[2] -= delz2*tk1/r2;

  if (EVFLAG && eflag) eangle += d_bb_k[type]*dr1_bb*dr2_bb;

  // bond-angle term

  const KK_FLOAT dr1_ba = r1 - d_ba_r1[type];
  const KK_FLOAT dr2_ba = r2 - d_ba_r2[type];
  const KK_FLOAT aa1 = s*dr1_ba*d_ba_k1[type];
  const KK_FLOAT aa2 = s*dr2_ba*d_ba_k2[type];

  const KK_FLOAT aa11 = aa1*c/rsq1;
  const KK_FLOAT aa12 = -aa1/(r1*r2);
  KK_FLOAT aa21 = aa2*c/rsq1;
  const KK_FLOAT aa22 = -aa2/(r1*r2);

  const KK_FLOAT vx11 = aa11*delx1 + aa12*delx2;
  const KK_FLOAT vx12 = aa21*delx1 + aa22*delx2;
  const KK_FLOAT vy11 = aa11*dely1 + aa12*dely2;
  const KK_FLOAT vy12 = aa21*dely1 + aa22*dely2;
  const KK_FLOAT vz11 = aa11*delz1 + aa12*delz2;
  const KK_FLOAT vz12 = aa21*delz1 + aa22*delz2;

  const KK_FLOAT aa11b = aa1*c/rsq2;
  aa21 = aa2*c/rsq2;

  const KK_FLOAT vx21 = aa11b*delx2 + aa12*delx1;
  const KK_FLOAT vx22 = aa21*delx2 + aa22*delx1;
  const KK_FLOAT vy21 = aa11b*dely2 + aa12*dely1;
  const KK_FLOAT vy22 = aa21*dely2 + aa22*dely1;
  const KK_FLOAT vz21 = aa11b*delz2 + aa12*delz1;
  const KK_FLOAT vz22 = aa21*delz2 + aa22*delz1;

  const KK_FLOAT b1 = d_ba_k1[type]*dtheta/r1;
  const KK_FLOAT b2 = d_ba_k2[type]*dtheta/r2;

  f1[0] -= vx11 + b1*delx1 + vx12;
  f1[1] -= vy11 + b1*dely1 + vy12;
  f1[2] -= vz11 + b1*delz1 + vz12;
  f3[0] -= vx21 + b2*delx2 + vx22;
  f3[1] -= vy21 + b2*dely2 + vy22;
  f3[2] -= vz21 + b2*delz2 + vz22;

  if (EVFLAG && eflag)
    eangle += d_ba_k1[type]*dr1_ba*dtheta + d_ba_k2[type]*dr2_ba*dtheta;

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
void AngleClass2P6Kokkos<DeviceType>::operator()(TagAngleClass2P6Compute<NEWTON_BOND,EVFLAG>,
                                                 const int &n) const
{
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleClass2P6Compute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleClass2P6Kokkos<DeviceType>::allocate()
{
  AngleClass2P6::allocate();
  int n = atom->nangletypes;

  k_theta0 = DAT::tdual_kkfloat_1d("AngleClass2P6::theta0",n+1);
  k_k2  = DAT::tdual_kkfloat_1d("AngleClass2P6::k2",n+1);
  k_k3  = DAT::tdual_kkfloat_1d("AngleClass2P6::k3",n+1);
  k_k4  = DAT::tdual_kkfloat_1d("AngleClass2P6::k4",n+1);
  k_k5  = DAT::tdual_kkfloat_1d("AngleClass2P6::k5",n+1);
  k_k6  = DAT::tdual_kkfloat_1d("AngleClass2P6::k6",n+1);
  k_bb_k  = DAT::tdual_kkfloat_1d("AngleClass2P6::bb_k",n+1);
  k_bb_r1 = DAT::tdual_kkfloat_1d("AngleClass2P6::bb_r1",n+1);
  k_bb_r2 = DAT::tdual_kkfloat_1d("AngleClass2P6::bb_r2",n+1);
  k_ba_k1 = DAT::tdual_kkfloat_1d("AngleClass2P6::ba_k1",n+1);
  k_ba_k2 = DAT::tdual_kkfloat_1d("AngleClass2P6::ba_k2",n+1);
  k_ba_r1 = DAT::tdual_kkfloat_1d("AngleClass2P6::ba_r1",n+1);
  k_ba_r2 = DAT::tdual_kkfloat_1d("AngleClass2P6::ba_r2",n+1);

  d_theta0 = k_theta0.template view<DeviceType>();
  d_k2  = k_k2.template view<DeviceType>();
  d_k3  = k_k3.template view<DeviceType>();
  d_k4  = k_k4.template view<DeviceType>();
  d_k5  = k_k5.template view<DeviceType>();
  d_k6  = k_k6.template view<DeviceType>();
  d_bb_k  = k_bb_k.template view<DeviceType>();
  d_bb_r1 = k_bb_r1.template view<DeviceType>();
  d_bb_r2 = k_bb_r2.template view<DeviceType>();
  d_ba_k1 = k_ba_k1.template view<DeviceType>();
  d_ba_k2 = k_ba_k2.template view<DeviceType>();
  d_ba_r1 = k_ba_r1.template view<DeviceType>();
  d_ba_r2 = k_ba_r2.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleClass2P6Kokkos<DeviceType>::coeff(int narg, char **arg)
{
  AngleClass2P6::coeff(narg, arg);
  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);
  for (int i = ilo; i <= ihi; i++) {
    k_theta0.view_host()[i] = static_cast<KK_FLOAT>(theta0[i]);
    k_k2.view_host()[i]  = static_cast<KK_FLOAT>(k2[i]);
    k_k3.view_host()[i]  = static_cast<KK_FLOAT>(k3[i]);
    k_k4.view_host()[i]  = static_cast<KK_FLOAT>(k4[i]);
    k_k5.view_host()[i]  = static_cast<KK_FLOAT>(k5[i]);
    k_k6.view_host()[i]  = static_cast<KK_FLOAT>(k6[i]);
    k_bb_k.view_host()[i]  = static_cast<KK_FLOAT>(bb_k[i]);
    k_bb_r1.view_host()[i] = static_cast<KK_FLOAT>(bb_r1[i]);
    k_bb_r2.view_host()[i] = static_cast<KK_FLOAT>(bb_r2[i]);
    k_ba_k1.view_host()[i] = static_cast<KK_FLOAT>(ba_k1[i]);
    k_ba_k2.view_host()[i] = static_cast<KK_FLOAT>(ba_k2[i]);
    k_ba_r1.view_host()[i] = static_cast<KK_FLOAT>(ba_r1[i]);
    k_ba_r2.view_host()[i] = static_cast<KK_FLOAT>(ba_r2[i]);
  }
  k_theta0.modify_host(); k_k2.modify_host(); k_k3.modify_host();
  k_k4.modify_host(); k_k5.modify_host(); k_k6.modify_host();
  k_bb_k.modify_host(); k_bb_r1.modify_host(); k_bb_r2.modify_host();
  k_ba_k1.modify_host(); k_ba_k2.modify_host();
  k_ba_r1.modify_host(); k_ba_r2.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleClass2P6Kokkos<DeviceType>::read_restart(FILE *fp)
{
  AngleClass2P6::read_restart(fp);
  int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_theta0.view_host()[i] = static_cast<KK_FLOAT>(theta0[i]);
    k_k2.view_host()[i]  = static_cast<KK_FLOAT>(k2[i]);
    k_k3.view_host()[i]  = static_cast<KK_FLOAT>(k3[i]);
    k_k4.view_host()[i]  = static_cast<KK_FLOAT>(k4[i]);
    k_k5.view_host()[i]  = static_cast<KK_FLOAT>(k5[i]);
    k_k6.view_host()[i]  = static_cast<KK_FLOAT>(k6[i]);
    k_bb_k.view_host()[i]  = static_cast<KK_FLOAT>(bb_k[i]);
    k_bb_r1.view_host()[i] = static_cast<KK_FLOAT>(bb_r1[i]);
    k_bb_r2.view_host()[i] = static_cast<KK_FLOAT>(bb_r2[i]);
    k_ba_k1.view_host()[i] = static_cast<KK_FLOAT>(ba_k1[i]);
    k_ba_k2.view_host()[i] = static_cast<KK_FLOAT>(ba_k2[i]);
    k_ba_r1.view_host()[i] = static_cast<KK_FLOAT>(ba_r1[i]);
    k_ba_r2.view_host()[i] = static_cast<KK_FLOAT>(ba_r2[i]);
  }
  k_theta0.modify_host(); k_k2.modify_host(); k_k3.modify_host();
  k_k4.modify_host(); k_k5.modify_host(); k_k6.modify_host();
  k_bb_k.modify_host(); k_bb_r1.modify_host(); k_bb_r2.modify_host();
  k_ba_k1.modify_host(); k_ba_k2.modify_host();
  k_ba_r1.modify_host(); k_ba_r2.modify_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleClass2P6Kokkos<DeviceType>::ev_tally(EV_FLOAT &ev,
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
template class AngleClass2P6Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class AngleClass2P6Kokkos<LMPHostType>;
#endif
}
