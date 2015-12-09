/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include "bond_harmonic_kokkos.h"
#include "atom_kokkos.h"
#include "neighbor_kokkos.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "atom_masks.h"

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
    memory->destroy_kokkos(k_eatom,eatom);
    memory->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void BondHarmonicKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memory->destroy_kokkos(k_eatom,eatom);
    memory->create_kokkos(k_eatom,eatom,maxeatom,"bond:eatom");
    d_eatom = k_eatom.d_view;
  }
  if (vflag_atom) {
    memory->destroy_kokkos(k_vatom,vatom);
    memory->create_kokkos(k_vatom,vatom,maxvatom,6,"bond:vatom");
    d_vatom = k_vatom.d_view;
  }

  atomKK->sync(execution_space,datamask_read);
  k_k.template sync<DeviceType>();
  k_r0.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  neighborKK->k_bondlist.template sync<DeviceType>();
  bondlist = neighborKK->k_bondlist.view<DeviceType>();
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
  DeviceType::fence();

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

  // The f array is atomic
  Kokkos::View<F_FLOAT*[3], typename DAT::t_f_array::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > a_f = f;

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
    a_f(i1,0) += delx*fbond;
    a_f(i1,1) += dely*fbond;
    a_f(i1,2) += delz*fbond;
  }

  if (NEWTON_BOND || i2 < nlocal) {
    a_f(i2,0) -= delx*fbond;
    a_f(i2,1) -= dely*fbond;
    a_f(i2,2) -= delz*fbond;
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

  int n = atom->nbondtypes;
  k_k = DAT::tdual_ffloat_1d("BondHarmonic::k",n+1);
  k_r0 = DAT::tdual_ffloat_1d("BondHarmonic::r0",n+1);

  d_k = k_k.d_view;
  d_r0 = k_r0.d_view;
}

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

template<class DeviceType>
void BondHarmonicKokkos<DeviceType>::coeff(int narg, char **arg)
{
  BondHarmonic::coeff(narg, arg);

  int n = atom->nbondtypes;
  for (int i = 1; i <= n; i++) {
    k_k.h_view[i] = k[i];
    k_r0.h_view[i] = r0[i];
  }

  k_k.template modify<LMPHostType>();
  k_r0.template modify<LMPHostType>();
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

  // The eatom and vatom arrays are atomic
  Kokkos::View<E_FLOAT*, typename DAT::t_efloat_1d::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = k_eatom.view<DeviceType>();
  Kokkos::View<F_FLOAT*[6], typename DAT::t_virial_array::array_layout,DeviceType,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = k_vatom.view<DeviceType>();

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
      if (newton_bond || i < nlocal) v_eatom[i] += ebondhalf;
      if (newton_bond || j < nlocal) v_eatom[j] += ebondhalf;
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
        v_vatom(i,0) += 0.5*v[0];
        v_vatom(i,1) += 0.5*v[1];
        v_vatom(i,2) += 0.5*v[2];
        v_vatom(i,3) += 0.5*v[3];
        v_vatom(i,4) += 0.5*v[4];
        v_vatom(i,5) += 0.5*v[5];
      }
      if (newton_bond || j < nlocal) {
        v_vatom(j,0) += 0.5*v[0];
        v_vatom(j,1) += 0.5*v[1];
        v_vatom(j,2) += 0.5*v[2];
        v_vatom(j,3) += 0.5*v[3];
        v_vatom(j,4) += 0.5*v[4];
        v_vatom(j,5) += 0.5*v[5];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template class BondHarmonicKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class BondHarmonicKokkos<LMPHostType>;
#endif