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
    //if(k_eatom.extent(0)<maxeatom) { // won't work without adding zero functor
      memoryKK->destroy_kokkos(k_eatom,eatom);
      memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"improper:eatom");
      d_eatom = k_eatom.template view<KKDeviceType>();
    //}
  }
  if (vflag_atom) {
    //if(k_vatom.extent(0)<maxvatom) { // won't work without adding zero functor
      memoryKK->destroy_kokkos(k_vatom,vatom);
      memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"improper:vatom");
      d_vatom = k_vatom.template view<KKDeviceType>();
    //}
  }

//  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
//  else atomKK->modified(execution_space,F_MASK);

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

  if (eflag_global) energy += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
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

  const F_FLOAT delx = x(i1,0) - x(i2,0);
  const F_FLOAT dely = x(i1,1) - x(i2,1);
  const F_FLOAT delz = x(i1,2) - x(i2,2);

  const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
  const F_FLOAT r = sqrt(rsq);
  const F_FLOAT dr = r - d_r0[type];
  const F_FLOAT rk = d_k[type] * dr;

  // force & energy

  F_FLOAT fbond = 0.0;
  if (r > 0.0) fbond = -2.0*rk/r;

  F_FLOAT ebond = 0.0;
  if (eflag)
    ebond = rk*dr;

  // apply force to each of 2 atoms

  if (NEWTON_BOND || i1 < nlocal) {
    f(i1,0) += delx*fbond;
    f(i1,1) += dely*fbond;
    f(i1,2) += delz*fbond;
  }

  if (NEWTON_BOND || i2 < nlocal) {
    f(i2,0) -= delx*fbond;
    f(i2,1) -= dely*fbond;
    f(i2,2) -= delz*fbond;
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
  typename AT::tdual_ffloat_1d k_k("BondHarmonic::k",n+1);
  typename AT::tdual_ffloat_1d k_r0("BondHarmonic::r0",n+1);

  d_k = k_k.template view<DeviceType>();
  d_r0 = k_r0.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_r0.h_view[i] = r0[i];
  }

  k_k.template modify<LMPHostType>();
  k_r0.template modify<LMPHostType>();
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
  typename AT::tdual_ffloat_1d k_k("BondHarmonic::k",n+1);
  typename AT::tdual_ffloat_1d k_r0("BondHarmonic::r0",n+1);

  d_k = k_k.template view<DeviceType>();
  d_r0 = k_r0.template view<DeviceType>();

  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_r0.h_view[i] = r0[i];
  }

  k_k.template modify<LMPHostType>();
  k_r0.template modify<LMPHostType>();
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
      const F_FLOAT &ebond, const F_FLOAT &fbond, const F_FLOAT &delx,
                const F_FLOAT &dely, const F_FLOAT &delz) const
{
  E_FLOAT ebondhalf;
  F_FLOAT v[6];

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += ebond;
      else {
        ebondhalf = 0.5*ebond;
        if (i < nlocal) ev.evdwl += ebondhalf;
        if (j < nlocal) ev.evdwl += ebondhalf;
      }
    }
    if (eflag_atom) {
      ebondhalf = 0.5*ebond;
      if (newton_bond || i < nlocal) d_eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) d_eatom[j] += ebondhalf;
    }
  }

  if (vflag_either) {
    v[0] = delx*delx*fbond;
    v[1] = dely*dely*fbond;
    v[2] = delz*delz*fbond;
    v[3] = delx*dely*fbond;
    v[4] = delx*delz*fbond;
    v[5] = dely*delz*fbond;

    if (vflag_global) {
      if (newton_bond) {
        ev.v[0] += v[0];
        ev.v[1] += v[1];
        ev.v[2] += v[2];
        ev.v[3] += v[3];
        ev.v[4] += v[4];
        ev.v[5] += v[5];
      } else {
        if (i < nlocal) {
          ev.v[0] += 0.5*v[0];
          ev.v[1] += 0.5*v[1];
          ev.v[2] += 0.5*v[2];
          ev.v[3] += 0.5*v[3];
          ev.v[4] += 0.5*v[4];
          ev.v[5] += 0.5*v[5];
        }
        if (j < nlocal) {
          ev.v[0] += 0.5*v[0];
          ev.v[1] += 0.5*v[1];
          ev.v[2] += 0.5*v[2];
          ev.v[3] += 0.5*v[3];
          ev.v[4] += 0.5*v[4];
          ev.v[5] += 0.5*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        d_vatom(i,0) += 0.5*v[0];
        d_vatom(i,1) += 0.5*v[1];
        d_vatom(i,2) += 0.5*v[2];
        d_vatom(i,3) += 0.5*v[3];
        d_vatom(i,4) += 0.5*v[4];
        d_vatom(i,5) += 0.5*v[5];
      }
      if (newton_bond || j < nlocal) {
        d_vatom(j,0) += 0.5*v[0];
        d_vatom(j,1) += 0.5*v[1];
        d_vatom(j,2) += 0.5*v[2];
        d_vatom(j,3) += 0.5*v[3];
        d_vatom(j,4) += 0.5*v[4];
        d_vatom(j,5) += 0.5*v[5];
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

