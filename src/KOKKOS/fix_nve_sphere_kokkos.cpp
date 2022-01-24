// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_nve_sphere_kokkos.h"
#include "atom_masks.h"
#include "atom_kokkos.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixNVESphereKokkos<DeviceType>::FixNVESphereKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixNVESphere(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereKokkos<DeviceType>::cleanup_copy()
{
  id = style = nullptr;
  vatom = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereKokkos<DeviceType>::init()
{
  FixNVESphere::init();

  if (extra == DIPOLE) {
    error->all(FLERR,"Fix nve/sphere/kk doesn't yet support dipole");
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereKokkos<DeviceType>::initial_integrate(int /*vflag*/)
{
  atomKK->sync(execution_space, X_MASK | V_MASK | OMEGA_MASK| F_MASK | TORQUE_MASK | RMASS_MASK | RADIUS_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVESphereKokkosInitialIntegrateFunctor<DeviceType> f(this);
  Kokkos::parallel_for(nlocal,f);

  atomKK->modified(execution_space,  X_MASK | V_MASK | OMEGA_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVESphereKokkos<DeviceType>::initial_integrate_item(const int i) const
{
  const double dtfrotate = dtf / inertia;

  if (mask(i) & groupbit) {
    const double dtfm = dtf / rmass(i);
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);
    x(i,0) += dtv * v(i,0);
    x(i,1) += dtv * v(i,1);
    x(i,2) += dtv * v(i,2);

    const double dtirotate = dtfrotate / (radius(i)*radius(i)*rmass(i));
    omega(i,0) += dtirotate * torque(i,0);
    omega(i,1) += dtirotate * torque(i,1);
    omega(i,2) += dtirotate * torque(i,2);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixNVESphereKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(execution_space, V_MASK | OMEGA_MASK| F_MASK | TORQUE_MASK | RMASS_MASK | RADIUS_MASK | MASK_MASK);

  v = atomKK->k_v.view<DeviceType>();
  omega = atomKK->k_omega.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();

  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  FixNVESphereKokkosFinalIntegrateFunctor<DeviceType> f(this);
  Kokkos::parallel_for(nlocal,f);

  atomKK->modified(execution_space, V_MASK | OMEGA_MASK);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixNVESphereKokkos<DeviceType>::final_integrate_item(const int i) const
{
  const double dtfrotate = dtf / inertia;

  if (mask(i) & groupbit) {
    const double dtfm = dtf / rmass(i);
    v(i,0) += dtfm * f(i,0);
    v(i,1) += dtfm * f(i,1);
    v(i,2) += dtfm * f(i,2);

    const double dtirotate = dtfrotate / (radius(i)*radius(i)*rmass(i));
    omega(i,0) += dtirotate * torque(i,0);
    omega(i,1) += dtirotate * torque(i,1);
    omega(i,2) += dtirotate * torque(i,2);
  }
}

namespace LAMMPS_NS {
template class FixNVESphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNVESphereKokkos<LMPHostType>;
#endif
}
