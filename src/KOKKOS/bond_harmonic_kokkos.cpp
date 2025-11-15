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

#include "bond_harmonic_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "force.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondHarmonicKokkos<DeviceType>::BondHarmonicKokkos(LAMMPS *lmp) : BondHarmonic(lmp)
{
  kokkosable = 1;

  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
BondHarmonicKokkos<DeviceType>::~BondHarmonicKokkos()
{
  if (!copymode) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondHarmonicKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  neighborKK->k_bondlist.template sync<DeviceType>();
  bondlist = neighborKK->k_bondlist.template view<DeviceType>();
  int nbondlist = neighborKK->nbondlist;
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;

  copymode = 1;

  // loop over neighbors of my atoms

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondHarmonicCompute<1,1> >(0,nbondlist),*this,ev);
    } else {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagBondHarmonicCompute<0,1> >(0,nbondlist),*this,ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondHarmonicCompute<1,0> >(0,nbondlist),*this);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagBondHarmonicCompute<0,0> >(0,nbondlist),*this);
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
KOKKOS_INLINE_FUNCTION
void BondHarmonicKokkos<DeviceType>::operator()(TagBondHarmonicCompute<NEWTON_BOND,EVFLAG>, const int &n, EV_FLOAT& ev) const {

  const int i1 = bondlist(n,0);
  const int i2 = bondlist(n,1);
  const int type = bondlist(n,2);

  const KK_FLOAT delx = x(i1,0) - x(i2,0);
  const KK_FLOAT dely = x(i1,1) - x(i2,1);
  const KK_FLOAT delz = x(i1,2) - x(i2,2);

  const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
  const KK_FLOAT r = sqrt(rsq);
  const KK_FLOAT dr = r - d_r0[type];
  const KK_FLOAT rk = d_k[type] * dr;

  // force & energy

  KK_FLOAT fbond = 0;
  if (r > 0) fbond = -static_cast<KK_FLOAT>(2.0) * rk / r;

  KK_FLOAT ebond = 0;
  if (eflag)
    ebond = static_cast<KK_FLOAT>(rk*dr);

  // apply force to each of 2 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    f(i1,0) += static_cast<KK_ACC_FLOAT>(delx*fbond);
    f(i1,1) += static_cast<KK_ACC_FLOAT>(dely*fbond);
    f(i1,2) += static_cast<KK_ACC_FLOAT>(delz*fbond);
  }

  if (NEWTON_BOND || i2 < nlocal) {
    f(i2,0) -= static_cast<KK_ACC_FLOAT>(delx*fbond);
    f(i2,1) -= static_cast<KK_ACC_FLOAT>(dely*fbond);
    f(i2,2) -= static_cast<KK_ACC_FLOAT>(delz*fbond);
  }

  if (EVFLAG) ev_tally(ev,i1,i2,ebond,fbond,delx,dely,delz);
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void BondHarmonicKokkos<DeviceType>::operator()(TagBondHarmonicCompute<NEWTON_BOND,EVFLAG>, const int &n) const {
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND,EVFLAG>(TagBondHarmonicCompute<NEWTON_BOND,EVFLAG>(), n, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondHarmonicKokkos<DeviceType>::allocate()
{
  BondHarmonic::allocate();
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void BondHarmonicKokkos<DeviceType>::coeff(int narg, char **arg)
{
  BondHarmonic::coeff(narg, arg);

  int n = atom->nbondtypes;
  DAT::tdual_kkfloat_1d k_k("BondHarmonic::k",n+1);
  DAT::tdual_kkfloat_1d k_r0("BondHarmonic::r0",n+1);

  d_k = k_k.template view<DeviceType>();
  d_r0 = k_r0.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k.view_host()[i] = static_cast<KK_FLOAT>(k[i]);
    k_r0.view_host()[i] = static_cast<KK_FLOAT>(r0[i]);
  }

  k_k.modify_host();
  k_r0.modify_host();
  k_k.template sync<DeviceType>();
  k_r0.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

template<class DeviceType>
void BondHarmonicKokkos<DeviceType>::read_restart(FILE *fp)
{
  BondHarmonic::read_restart(fp);

  int n = atom->nbondtypes;
  DAT::tdual_kkfloat_1d k_k("BondHarmonic::k",n+1);
  DAT::tdual_kkfloat_1d k_r0("BondHarmonic::r0",n+1);

  d_k = k_k.template view<DeviceType>();
  d_r0 = k_r0.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k.view_host()[i] = static_cast<KK_FLOAT>(k[i]);
    k_r0.view_host()[i] = static_cast<KK_FLOAT>(r0[i]);
  }

  k_k.modify_host();
  k_r0.modify_host();
  k_k.template sync<DeviceType>();
  k_r0.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators
------------------------------------------------------------------------- */

template<class DeviceType>
//template<int NEWTON_BOND>
KOKKOS_INLINE_FUNCTION
void BondHarmonicKokkos<DeviceType>::ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &ebond, const KK_FLOAT &fbond, const KK_FLOAT &delx,
                const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(ebond);
      else {
        KK_ACC_FLOAT ebondhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
        if (i < nlocal) ev.evdwl += static_cast<KK_ACC_FLOAT>(ebondhalf);
        if (j < nlocal) ev.evdwl += static_cast<KK_ACC_FLOAT>(ebondhalf);
      }
    }
    if (eflag_atom) {
      KK_ACC_FLOAT ebondhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*ebond);
      if (newton_bond || i < nlocal) d_eatom[i] += static_cast<KK_ACC_FLOAT>(ebondhalf);
      if (newton_bond || j < nlocal) d_eatom[j] += static_cast<KK_ACC_FLOAT>(ebondhalf);
    }
  }

  if (vflag_either) {
    KK_ACC_FLOAT v_half_acc[6];
    v_half_acc[0] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delx*fbond);
    v_half_acc[1] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*dely*fbond);
    v_half_acc[2] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delz*delz*fbond);
    v_half_acc[3] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*dely*fbond);
    v_half_acc[4] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*delx*delz*fbond);
    v_half_acc[5] = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*dely*delz*fbond);

    if (vflag_global) {
      if (newton_bond) {
        for (int n = 0; n < 6; n++)
          ev.v[n] += static_cast<KK_ACC_FLOAT>(2.0)*v_half_acc[n];
      } else {
        if (i < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_half_acc[n];
        }
        if (j < nlocal) {
          for (int n = 0; n < 6; n++)
            ev.v[n] += v_half_acc[n];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        for (int n = 0; n < 6; n++)
          d_vatom(i,n) += v_half_acc[n];
      }
      if (newton_bond || j < nlocal) {
        for (int n = 0; n < 6; n++)
          d_vatom(j,n) += v_half_acc[n];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class BondHarmonicKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class BondHarmonicKokkos<LMPHostType>;
#endif
}

