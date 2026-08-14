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

#include "angle_cross_kokkos.h"

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
AngleCrossKokkos<DeviceType>::AngleCrossKokkos(LAMMPS *lmp) : AngleCross(lmp)
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
AngleCrossKokkos<DeviceType>::~AngleCrossKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleCrossKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  if (eflag_atom) {
    if ((int)k_eatom.extent(0) < maxeatom) {
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"angle:eatom");
      d_eatom = k_eatom.template view<DeviceType>();
    } else Kokkos::deep_copy(d_eatom,0.0);
  }
  if (vflag_atom) {
    if ((int)k_vatom.extent(0) < maxvatom) {
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"angle:vatom");
      d_vatom = k_vatom.template view<DeviceType>();
    } else Kokkos::deep_copy(d_vatom,0.0);
  }

  k_kss.template sync<DeviceType>();
  k_kbs0.template sync<DeviceType>();
  k_kbs1.template sync<DeviceType>();
  k_r00.template sync<DeviceType>();
  k_r01.template sync<DeviceType>();
  k_theta0.template sync<DeviceType>();

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
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleCrossCompute<1,1> >(0,nanglelist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleCrossCompute<0,1> >(0,nanglelist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleCrossCompute<1,0> >(0,nanglelist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleCrossCompute<0,0> >(0,nanglelist),*this);
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
void AngleCrossKokkos<DeviceType>::operator()(TagAngleCrossCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

  const int i1 = anglelist(n,0);
  const int i2 = anglelist(n,1);
  const int i3 = anglelist(n,2);
  const int type = anglelist(n,3);

  const KK_FLOAT delx1 = x(i1,0) - x(i2,0);
  const KK_FLOAT dely1 = x(i1,1) - x(i2,1);
  const KK_FLOAT delz1 = x(i1,2) - x(i2,2);

  const KK_FLOAT rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
  const KK_FLOAT r1 = sqrt(rsq1);

  const KK_FLOAT delx2 = x(i3,0) - x(i2,0);
  const KK_FLOAT dely2 = x(i3,1) - x(i2,1);
  const KK_FLOAT delz2 = x(i3,2) - x(i2,2);

  const KK_FLOAT rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
  const KK_FLOAT r2 = sqrt(rsq2);

  KK_FLOAT c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;

  if (c > static_cast<KK_FLOAT>(1.0)) c = static_cast<KK_FLOAT>(1.0);
  if (c < static_cast<KK_FLOAT>(-1.0)) c = static_cast<KK_FLOAT>(-1.0);

  KK_FLOAT s = sqrt(static_cast<KK_FLOAT>(1.0) - c*c);
  if (s < static_cast<KK_FLOAT>(SMALL)) s = static_cast<KK_FLOAT>(SMALL);
  s = static_cast<KK_FLOAT>(1.0)/s;

  // bond-bond contribution
  const KK_FLOAT dr1 = r1 - d_r00[type];
  const KK_FLOAT dr2 = r2 - d_r01[type];
  const KK_FLOAT tk1 = d_kss[type] * dr1;
  const KK_FLOAT tk2 = d_kss[type] * dr2;

  KK_FLOAT f1[3],f3[3];
  f1[0] = -delx1*tk2/r1;
  f1[1] = -dely1*tk2/r1;
  f1[2] = -delz1*tk2/r1;

  f3[0] = -delx2*tk1/r2;
  f3[1] = -dely2*tk1/r2;
  f3[2] = -delz2*tk1/r2;

  KK_FLOAT eangle = static_cast<KK_FLOAT>(0.0);
  if (eflag) eangle = d_kss[type]*dr1*dr2;

  // bond-angle contribution
  const KK_FLOAT dtheta = acos(c) - d_theta0[type];

  const KK_FLOAT aa1 = s * dr1 * d_kbs0[type];
  const KK_FLOAT aa2 = s * dr2 * d_kbs1[type];

  KK_FLOAT aa11 = aa1 * c / rsq1;
  const KK_FLOAT aa12 = -aa1 / (r1 * r2);
  KK_FLOAT aa21 = aa2 * c / rsq1;
  const KK_FLOAT aa22 = -aa2 / (r1 * r2);

  const KK_FLOAT vx11 = aa11*delx1 + aa12*delx2;
  const KK_FLOAT vx12 = aa21*delx1 + aa22*delx2;
  const KK_FLOAT vy11 = aa11*dely1 + aa12*dely2;
  const KK_FLOAT vy12 = aa21*dely1 + aa22*dely2;
  const KK_FLOAT vz11 = aa11*delz1 + aa12*delz2;
  const KK_FLOAT vz12 = aa21*delz1 + aa22*delz2;

  aa11 = aa1 * c / rsq2;
  aa21 = aa2 * c / rsq2;

  const KK_FLOAT vx21 = aa11*delx2 + aa12*delx1;
  const KK_FLOAT vx22 = aa21*delx2 + aa22*delx1;
  const KK_FLOAT vy21 = aa11*dely2 + aa12*dely1;
  const KK_FLOAT vy22 = aa21*dely2 + aa22*dely1;
  const KK_FLOAT vz21 = aa11*delz2 + aa12*delz1;
  const KK_FLOAT vz22 = aa21*delz2 + aa22*delz1;

  const KK_FLOAT b1 = d_kbs0[type] * dtheta / r1;
  const KK_FLOAT b2 = d_kbs1[type] * dtheta / r2;

  f1[0] -= vx11 + b1*delx1 + vx12;
  f1[1] -= vy11 + b1*dely1 + vy12;
  f1[2] -= vz11 + b1*delz1 + vz12;

  f3[0] -= vx21 + b2*delx2 + vx22;
  f3[1] -= vy21 + b2*dely2 + vy22;
  f3[2] -= vz21 + b2*delz2 + vz22;

  if (eflag) eangle += d_kbs0[type]*dr1*dtheta + d_kbs1[type]*dr2*dtheta;

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

  if (EVFLAG) ev_tally(ev,i1,i2,i3,eangle,f1,f3,
                       delx1,dely1,delz1,delx2,dely2,delz2);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleCrossKokkos<DeviceType>::operator()(TagAngleCrossCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagAngleCrossCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleCrossKokkos<DeviceType>::allocate()
{
  AngleCross::allocate();

  int n = atom->nangletypes;
  k_kss    = DAT::tdual_kkfloat_1d("AngleCross::kss",n+1);
  k_kbs0   = DAT::tdual_kkfloat_1d("AngleCross::kbs0",n+1);
  k_kbs1   = DAT::tdual_kkfloat_1d("AngleCross::kbs1",n+1);
  k_r00    = DAT::tdual_kkfloat_1d("AngleCross::r00",n+1);
  k_r01    = DAT::tdual_kkfloat_1d("AngleCross::r01",n+1);
  k_theta0 = DAT::tdual_kkfloat_1d("AngleCross::theta0",n+1);

  d_kss    = k_kss.template view<DeviceType>();
  d_kbs0   = k_kbs0.template view<DeviceType>();
  d_kbs1   = k_kbs1.template view<DeviceType>();
  d_r00    = k_r00.template view<DeviceType>();
  d_r01    = k_r01.template view<DeviceType>();
  d_theta0 = k_theta0.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleCrossKokkos<DeviceType>::coeff(int narg, char **arg)
{
  AngleCross::coeff(narg, arg);

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  for (int i = ilo; i <= ihi; i++) {
    k_kss.view_host()[i]    = static_cast<KK_FLOAT>(kss[i]);
    k_kbs0.view_host()[i]   = static_cast<KK_FLOAT>(kbs0[i]);
    k_kbs1.view_host()[i]   = static_cast<KK_FLOAT>(kbs1[i]);
    k_r00.view_host()[i]    = static_cast<KK_FLOAT>(r00[i]);
    k_r01.view_host()[i]    = static_cast<KK_FLOAT>(r01[i]);
    k_theta0.view_host()[i] = static_cast<KK_FLOAT>(theta0[i]);
  }

  k_kss.modify_host();
  k_kbs0.modify_host();
  k_kbs1.modify_host();
  k_r00.modify_host();
  k_r01.modify_host();
  k_theta0.modify_host();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleCrossKokkos<DeviceType>::read_restart(FILE *fp)
{
  AngleCross::read_restart(fp);

  int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_kss.view_host()[i]    = static_cast<KK_FLOAT>(kss[i]);
    k_kbs0.view_host()[i]   = static_cast<KK_FLOAT>(kbs0[i]);
    k_kbs1.view_host()[i]   = static_cast<KK_FLOAT>(kbs1[i]);
    k_r00.view_host()[i]    = static_cast<KK_FLOAT>(r00[i]);
    k_r01.view_host()[i]    = static_cast<KK_FLOAT>(r01[i]);
    k_theta0.view_host()[i] = static_cast<KK_FLOAT>(theta0[i]);
  }

  k_kss.modify_host();
  k_kbs0.modify_host();
  k_kbs1.modify_host();
  k_r00.modify_host();
  k_r01.modify_host();
  k_theta0.modify_host();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleCrossKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
                     const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
                     const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const
{
  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = d_vatom;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(eangle);
      else {
        KK_ACC_FLOAT eanglethird = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*eangle);
        if (i < nlocal) ev.evdwl += eanglethird;
        if (j < nlocal) ev.evdwl += eanglethird;
        if (k < nlocal) ev.evdwl += eanglethird;
      }
    }
    if (eflag_atom) {
      KK_ACC_FLOAT eanglethird = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*eangle);
      if (newton_bond || i < nlocal) v_eatom[i] += static_cast<KK_ACC_FLOAT>(eanglethird);
      if (newton_bond || j < nlocal) v_eatom[j] += static_cast<KK_ACC_FLOAT>(eanglethird);
      if (newton_bond || k < nlocal) v_eatom[k] += static_cast<KK_ACC_FLOAT>(eanglethird);
    }
  }

  if (vflag_either) {
    KK_ACC_FLOAT v_third_acc[6];
    v_third_acc[0] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(delx1*f1[0] + delx2*f3[0]));
    v_third_acc[1] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(dely1*f1[1] + dely2*f3[1]));
    v_third_acc[2] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(delz1*f1[2] + delz2*f3[2]));
    v_third_acc[3] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(delx1*f1[1] + delx2*f3[1]));
    v_third_acc[4] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(delx1*f1[2] + delx2*f3[2]));
    v_third_acc[5] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(THIRD)*(dely1*f1[2] + dely2*f3[2]));

    if (vflag_global) {
      if (newton_bond) {
        for (int n = 0; n < 6; n++)
          ev.v[n] += static_cast<KK_ACC_FLOAT>(static_cast<KK_ACC_FLOAT>(3.0)*v_third_acc[n]);
      } else {
        if (i < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_third_acc[n];
        }
        if (j < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_third_acc[n];
        }
        if (k < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_third_acc[n];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        for (int n = 0; n < 6; n++)
          v_vatom(i,n) += v_third_acc[n];
      }
      if (newton_bond || j < nlocal) {
        for (int n = 0; n < 6; n++)
          v_vatom(j,n) += v_third_acc[n];
      }
      if (newton_bond || k < nlocal) {
        for (int n = 0; n < 6; n++)
          v_vatom(k,n) += v_third_acc[n];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class AngleCrossKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class AngleCrossKokkos<LMPHostType>;
#endif
}
