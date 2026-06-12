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

#include "angle_dipole_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double SMALL = 1.0e-100;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleDipoleKokkos<DeviceType>::AngleDipoleKokkos(LAMMPS *lmp) : AngleDipole(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TORQUE_MASK | MU_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  centroidstressflag = CENTROID_NOTAVAIL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
AngleDipoleKokkos<DeviceType>::~AngleDipoleKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleDipoleKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  if (!force->newton_bond)
    error->all(FLERR,"'newton' flag for bonded interactions must be 'on'");

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

  k_k.template sync<DeviceType>();
  k_gamma0.template sync<DeviceType>();

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  d_mu     = atomKK->k_mu.template view<DeviceType>();
  d_torque = atomKK->k_torque.template view<DeviceType>();
  neighborKK->k_anglelist.template sync<DeviceType>();
  anglelist = neighborKK->k_anglelist.template view<DeviceType>();
  int nanglelist = neighborKK->nanglelist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  EV_FLOAT ev;

  if (evflag) {
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagAngleDipoleCompute<1> >(0,nanglelist),*this,ev);
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagAngleDipoleCompute<0> >(0,nanglelist),*this);
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

  // torque has been modified on device; mark it
  atomKK->k_torque.template modify<DeviceType>();

  copymode = 0;
}

template<class DeviceType>
template<int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleDipoleKokkos<DeviceType>::operator()(TagAngleDipoleCompute<EVFLAG>, const int &n, EV_FLOAT& ev) const {

  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;
  Kokkos::View<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout,typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_torque = d_torque;

  const int iDip   = anglelist(n,0);
  const int iRef   = anglelist(n,1);
  const int iDummy = anglelist(n,2);
  const int type   = anglelist(n,3);

  const KK_FLOAT delx = x(iRef,0) - x(iDip,0);
  const KK_FLOAT dely = x(iRef,1) - x(iDip,1);
  const KK_FLOAT delz = x(iRef,2) - x(iDip,2);

  const KK_FLOAT r2 = delx*delx + dely*dely + delz*delz;
  if (r2 < static_cast<KK_FLOAT>(SMALL)*static_cast<KK_FLOAT>(SMALL)) return;
  const KK_FLOAT r = sqrt(r2);

  const KK_FLOAT rmu = r * d_mu(iDip,3);
  const KK_FLOAT cosGamma = (d_mu(iDip,0)*delx + d_mu(iDip,1)*dely + d_mu(iDip,2)*delz) / rmu;
  const KK_FLOAT deltaGamma = cosGamma - cos(d_gamma0[type]);
  const KK_FLOAT kdg = d_k[type] * deltaGamma;

  KK_FLOAT eangle = static_cast<KK_FLOAT>(0.0);
  if (eflag) eangle = kdg * deltaGamma;

  const KK_FLOAT tangle = static_cast<KK_FLOAT>(2.0) * kdg / rmu;

  const KK_FLOAT delTx = tangle * (dely*d_mu(iDip,2) - delz*d_mu(iDip,1));
  const KK_FLOAT delTy = tangle * (delz*d_mu(iDip,0) - delx*d_mu(iDip,2));
  const KK_FLOAT delTz = tangle * (delx*d_mu(iDip,1) - dely*d_mu(iDip,0));

  a_torque(iDip,0) += static_cast<KK_ACC_FLOAT>(delTx);
  a_torque(iDip,1) += static_cast<KK_ACC_FLOAT>(delTy);
  a_torque(iDip,2) += static_cast<KK_ACC_FLOAT>(delTz);

  // force couple that counterbalances dipolar torque
  const KK_FLOAT fx = dely*delTz - delz*delTy;
  const KK_FLOAT fy = delz*delTx - delx*delTz;
  const KK_FLOAT fz = delx*delTy - dely*delTx;

  const KK_FLOAT fmod2 = delTx*delTx + delTy*delTy + delTz*delTz;
  const KK_FLOAT len2  = fx*fx + fy*fy + fz*fz;
  if (len2 < static_cast<KK_FLOAT>(SMALL)*static_cast<KK_FLOAT>(SMALL)) return;

  const KK_FLOAT fmod_len = sqrt(fmod2) / (r * sqrt(len2));

  KK_FLOAT fi[3], fj[3];
  fi[0] = fx * fmod_len;
  fi[1] = fy * fmod_len;
  fi[2] = fz * fmod_len;
  fj[0] = -fi[0];
  fj[1] = -fi[1];
  fj[2] = -fi[2];

  // fj on iDip, fi on iRef (newton_bond always on for this style)
  a_f(iDip,0) += static_cast<KK_ACC_FLOAT>(fj[0]);
  a_f(iDip,1) += static_cast<KK_ACC_FLOAT>(fj[1]);
  a_f(iDip,2) += static_cast<KK_ACC_FLOAT>(fj[2]);
  a_f(iRef,0) += static_cast<KK_ACC_FLOAT>(fi[0]);
  a_f(iRef,1) += static_cast<KK_ACC_FLOAT>(fi[1]);
  a_f(iRef,2) += static_cast<KK_ACC_FLOAT>(fi[2]);

  // virial = r*F contributions are zero (all del vectors are zero in ev_tally)
  if (EVFLAG) {
    ev_tally(ev, iRef, iDip, iDummy, eangle, fj, fi,
             static_cast<KK_FLOAT>(0.0), static_cast<KK_FLOAT>(0.0), static_cast<KK_FLOAT>(0.0),
             static_cast<KK_FLOAT>(0.0), static_cast<KK_FLOAT>(0.0), static_cast<KK_FLOAT>(0.0));
  }
}

template<class DeviceType>
template<int EVFLAG>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleDipoleKokkos<DeviceType>::operator()(TagAngleDipoleCompute<EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<EVFLAG>(TagAngleDipoleCompute<EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void AngleDipoleKokkos<DeviceType>::allocate()
{
  AngleDipole::allocate();

  int n = atom->nangletypes;
  k_k      = DAT::tdual_kkfloat_1d("AngleDipole::k",n+1);
  k_gamma0 = DAT::tdual_kkfloat_1d("AngleDipole::gamma0",n+1);

  d_k      = k_k.template view<DeviceType>();
  d_gamma0 = k_gamma0.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleDipoleKokkos<DeviceType>::coeff(int narg, char **arg)
{
  AngleDipole::coeff(narg, arg);

  int ilo,ihi;
  utils::bounds(FLERR,arg[0],1,atom->nangletypes,ilo,ihi,error);

  for (int i = ilo; i <= ihi; i++) {
    k_k.view_host()[i]      = static_cast<KK_FLOAT>(k[i]);
    k_gamma0.view_host()[i] = static_cast<KK_FLOAT>(gamma0[i]);
  }

  k_k.modify_host();
  k_gamma0.modify_host();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void AngleDipoleKokkos<DeviceType>::read_restart(FILE *fp)
{
  AngleDipole::read_restart(fp);

  int n = atom->nangletypes;
  for (int i = 1; i <= n; i++) {
    k_k.view_host()[i]      = static_cast<KK_FLOAT>(k[i]);
    k_gamma0.view_host()[i] = static_cast<KK_FLOAT>(gamma0[i]);
  }

  k_k.modify_host();
  k_gamma0.modify_host();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
   virial = r1F1 + r2F2 + r3F3 = (r1-r2) F1 + (r3-r2) F3 = del1*f1 + del2*f3
   (for dipole, all del vectors are zero so virial contribution is zero)
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void AngleDipoleKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
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
template class AngleDipoleKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class AngleDipoleKokkos<LMPHostType>;
#endif
}
