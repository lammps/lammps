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

#include "stdio.h"
#include "string.h"
#include "fix_nve_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVEKokkos<DeviceType>::FixNVEKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = X_MASK | V_MASK | F_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEKokkos<DeviceType>::init()
{
  FixNVE::init();

  atomKK->k_mass.modify<LMPHostType>();
  atomKK->k_mass.sync<LMPDeviceType>();
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEKokkos<DeviceType>::initial_integrate(int vflag)
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  rmass = atomKK->rmass;
  mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  if (rmass) {
    FixNVEKokkosInitialIntegrateFunctor<DeviceType,1> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  } else {
    FixNVEKokkosInitialIntegrateFunctor<DeviceType,0> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  }
  DeviceType::fence();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVEKokkos<DeviceType>::initial_integrate_item(int i) const
{
  if (mask[i] & groupbit) {
    const double dtfm = dtf / mass[type[i]];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVEKokkos<DeviceType>::initial_integrate_rmass_item(int i) const
{
  if (mask[i] & groupbit) {
    const double dtfm = dtf / rmass[type[i]];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(execution_space,datamask_read);
  atomKK->modified(execution_space,datamask_modify);

  v = atomKK->k_v.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  rmass = atomKK->rmass;
  mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  if (rmass) {
    FixNVEKokkosFinalIntegrateFunctor<DeviceType,1> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  } else {
    FixNVEKokkosFinalIntegrateFunctor<DeviceType,0> functor(this);
    Kokkos::parallel_for(nlocal,functor);
  }
  DeviceType::fence();

  // debug
  //atomKK->sync(Host,datamask_read);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVEKokkos<DeviceType>::final_integrate_item(int i) const
{
  if (mask[i] & groupbit) {
    const double dtfm = dtf / mass[type[i]];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVEKokkos<DeviceType>::final_integrate_rmass_item(int i) const
{
  if (mask[i] & groupbit) {
    const double dtfm = dtf / rmass[type[i]];
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVEKokkos<DeviceType>::cleanup_copy()
{
  id = style = NULL;
  vatom = NULL;
}

template class FixNVEKokkos<LMPDeviceType>;
#if DEVICE==2
template class FixNVEKokkos<LMPHostType>;
#endif
