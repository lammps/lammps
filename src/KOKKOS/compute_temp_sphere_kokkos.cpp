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

static constexpr double INERTIA = 0.4;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeTempSphereKokkos<DeviceType>::ComputeTempSphereKokkos(LAMMPS *lmp, int narg, char **arg)
  : ComputeTempSphere(lmp, narg, arg)
{
  if (tempbias)
    error->all(FLERR, "compute temp/sphere/kk does not support bias option");

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = V_MASK | OMEGA_MASK | RADIUS_MASK | RMASS_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeTempSphereKokkos<DeviceType>::compute_scalar()
{
  atomKK->sync(execution_space, datamask_read);

  invoked_scalar = update->ntimestep;

  v      = atomKK->k_v.view<DeviceType>();
  omega  = atomKK->k_omega.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  rmass  = atomKK->k_rmass.view<DeviceType>();
  mask   = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  CTEMP t_kk;
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempSphereScalar>(0, nlocal),
                          *this, t_kk);
  copymode = 0;

  double t = t_kk.t0;
  MPI_Allreduce(&t, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  if (dynamic) dof_compute();
  if (dof < 0.0 && natoms_temp > 0.0)
    error->all(FLERR, "Temperature compute degrees of freedom < 0");
  scalar *= tfactor;
  return scalar;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeTempSphereKokkos<DeviceType>::operator()(TagComputeTempSphereScalar,
                                                      const int &i, CTEMP &t_kk) const
{
  if (mask[i] & groupbit) {
    t_kk.t0 += rmass[i] * (v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2));
    t_kk.t0 += INERTIA * rmass[i] * radius[i]*radius[i] *
               (omega(i,0)*omega(i,0) + omega(i,1)*omega(i,1) + omega(i,2)*omega(i,2));
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeTempSphereKokkos<DeviceType>::compute_vector()
{
  atomKK->sync(execution_space, datamask_read);

  invoked_vector = update->ntimestep;

  v      = atomKK->k_v.view<DeviceType>();
  omega  = atomKK->k_omega.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  rmass  = atomKK->k_rmass.view<DeviceType>();
  mask   = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  CTEMP t_kk;
  copymode = 1;
  Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagComputeTempSphereVector>(0, nlocal),
                          *this, t_kk);
  copymode = 0;

  double t[6] = {t_kk.t0, t_kk.t1, t_kk.t2, t_kk.t3, t_kk.t4, t_kk.t5};
  MPI_Allreduce(t, vector, 6, MPI_DOUBLE, MPI_SUM, world);
  for (int i = 0; i < 6; i++) vector[i] *= force->mvv2e;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeTempSphereKokkos<DeviceType>::operator()(TagComputeTempSphereVector,
                                                      const int &i, CTEMP &t_kk) const
{
  if (mask[i] & groupbit) {
    const double m  = rmass[i];
    const double ir = INERTIA * m * radius[i]*radius[i];
    t_kk.t0 += m  * v(i,0)*v(i,0) + ir * omega(i,0)*omega(i,0);
    t_kk.t1 += m  * v(i,1)*v(i,1) + ir * omega(i,1)*omega(i,1);
    t_kk.t2 += m  * v(i,2)*v(i,2) + ir * omega(i,2)*omega(i,2);
    t_kk.t3 += m  * v(i,0)*v(i,1) + ir * omega(i,0)*omega(i,1);
    t_kk.t4 += m  * v(i,0)*v(i,2) + ir * omega(i,0)*omega(i,2);
    t_kk.t5 += m  * v(i,1)*v(i,2) + ir * omega(i,1)*omega(i,2);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeTempSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeTempSphereKokkos<LMPHostType>;
#endif
}
