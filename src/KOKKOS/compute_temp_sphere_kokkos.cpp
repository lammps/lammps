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

#include "compute_temp_sphere_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

static constexpr double INERTIA = 0.4;    // moment of inertia prefactor for sphere

enum { ROTATE, ALL };

/* ---------------------------------------------------------------------- */

template <class DeviceType>
ComputeTempSphereKokkos<DeviceType>::ComputeTempSphereKokkos(LAMMPS *lmp, int narg, char **arg) :
    ComputeTempSphere(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read  = V_MASK | OMEGA_MASK | RMASS_MASK | RADIUS_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
double ComputeTempSphereKokkos<DeviceType>::compute_scalar()
{
  atomKK->sync(execution_space, datamask_read);

  invoked_scalar = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_scalar != update->ntimestep) tbias->compute_scalar();
    tbias->remove_bias_all();
    atomKK->sync(execution_space, V_MASK);
  }

  v        = atomKK->k_v.view<DeviceType>();
  omega_kk = atomKK->k_omega.view<DeviceType>();
  rmass_kk = atomKK->k_rmass.view<DeviceType>();
  radius_kk = atomKK->k_radius.view<DeviceType>();
  mask     = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  CTEMP t_kk;

  copymode = 1;
  if (mode == ALL)
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType, TagComputeTempSphereScalar<1>>(0, nlocal), *this, t_kk);
  else
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType, TagComputeTempSphereScalar<0>>(0, nlocal), *this, t_kk);
  copymode = 0;

  if (tempbias) {
    tbias->restore_bias_all();
    atomKK->sync(execution_space, V_MASK);
  }

  double t = t_kk.t0;
  MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  if (dynamic || tempbias == 2) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR, "Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

template <class DeviceType>
template <int MODE>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempSphereKokkos<DeviceType>::operator()(TagComputeTempSphereScalar<MODE>,
                                                     const int &i, CTEMP &t_kk) const
{
  if (mask[i] & groupbit) {
    const KK_FLOAT inertiaone = INERTIA * rmass_kk[i] * radius_kk[i] * radius_kk[i];
    if (MODE) {    // ALL: translational + rotational
      t_kk.t0 += (v(i, 0) * v(i, 0) + v(i, 1) * v(i, 1) + v(i, 2) * v(i, 2)) * rmass_kk[i];
    }
    t_kk.t0 += (omega_kk(i, 0) * omega_kk(i, 0) + omega_kk(i, 1) * omega_kk(i, 1) +
                omega_kk(i, 2) * omega_kk(i, 2)) *
        inertiaone;
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void ComputeTempSphereKokkos<DeviceType>::compute_vector()
{
  atomKK->sync(execution_space, datamask_read);

  invoked_vector = update->ntimestep;

  if (tempbias) {
    if (tbias->invoked_vector != update->ntimestep) tbias->compute_vector();
    tbias->remove_bias_all();
    atomKK->sync(execution_space, V_MASK);
  }

  v        = atomKK->k_v.view<DeviceType>();
  omega_kk = atomKK->k_omega.view<DeviceType>();
  rmass_kk = atomKK->k_rmass.view<DeviceType>();
  radius_kk = atomKK->k_radius.view<DeviceType>();
  mask     = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  CTEMP t_kk;

  copymode = 1;
  if (mode == ALL)
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType, TagComputeTempSphereVector<1>>(0, nlocal), *this, t_kk);
  else
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<DeviceType, TagComputeTempSphereVector<0>>(0, nlocal), *this, t_kk);
  copymode = 0;

  if (tempbias) {
    tbias->restore_bias_all();
    atomKK->sync(execution_space, V_MASK);
  }

  double t[6];
  t[0] = t_kk.t0;
  t[1] = t_kk.t1;
  t[2] = t_kk.t2;
  t[3] = t_kk.t3;
  t[4] = t_kk.t4;
  t[5] = t_kk.t5;
  MPI_Allreduce(t, vector, 6, MPI_DOUBLE, MPI_SUM, world);
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

template <class DeviceType>
template <int MODE>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeTempSphereKokkos<DeviceType>::operator()(TagComputeTempSphereVector<MODE>,
                                                     const int &i, CTEMP &t_kk) const
{
  if (mask[i] & groupbit) {
    const KK_FLOAT massone    = rmass_kk[i];
    const KK_FLOAT inertiaone = INERTIA * massone * radius_kk[i] * radius_kk[i];
    if (MODE) {    // ALL: translational + rotational
      t_kk.t0 += massone * v(i, 0) * v(i, 0);
      t_kk.t1 += massone * v(i, 1) * v(i, 1);
      t_kk.t2 += massone * v(i, 2) * v(i, 2);
      t_kk.t3 += massone * v(i, 0) * v(i, 1);
      t_kk.t4 += massone * v(i, 0) * v(i, 2);
      t_kk.t5 += massone * v(i, 1) * v(i, 2);
    }
    t_kk.t0 += inertiaone * omega_kk(i, 0) * omega_kk(i, 0);
    t_kk.t1 += inertiaone * omega_kk(i, 1) * omega_kk(i, 1);
    t_kk.t2 += inertiaone * omega_kk(i, 2) * omega_kk(i, 2);
    t_kk.t3 += inertiaone * omega_kk(i, 0) * omega_kk(i, 1);
    t_kk.t4 += inertiaone * omega_kk(i, 0) * omega_kk(i, 2);
    t_kk.t5 += inertiaone * omega_kk(i, 1) * omega_kk(i, 2);
  }
}

namespace LAMMPS_NS {
template class ComputeTempSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeTempSphereKokkos<LMPHostType>;
#endif
}
