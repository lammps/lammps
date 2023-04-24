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

#include "compute_erotate_sphere_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeERotateSphereKokkos<DeviceType>::ComputeERotateSphereKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeERotateSphere(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = OMEGA_MASK | RADIUS_MASK | MASK_MASK | RMASS_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double ComputeERotateSphereKokkos<DeviceType>::compute_scalar()
{
  atomKK->sync(execution_space,datamask_read);

  invoked_scalar = update->ntimestep;

  omega = atomKK->k_omega.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  // sum rotational energy for each particle
  // point particles will not contribute, due to radius = 0.0

  double erotate = 0.0;

  {
    // local variables for lambda capture

    auto l_omega = omega;
    auto l_radius = radius;
    auto l_rmass = rmass;
    auto l_mask = mask;
    auto l_groupbit = groupbit;

    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType>(0,nlocal), LAMMPS_LAMBDA(int i, double &erotate) {
      if (l_mask[i] & l_groupbit) {
        auto omega0 = l_omega(i,0);
        auto omega1 = l_omega(i,1);
        auto omega2 = l_omega(i,2);
        auto radius = l_radius(i);
        erotate +=
            (omega0 * omega0 + omega1 * omega1 + omega2 * omega2) *
            radius * radius * l_rmass[i];
      }
    },erotate);
  }

  MPI_Allreduce(&erotate, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  scalar *= pfactor;
  return scalar;
}

namespace LAMMPS_NS {
template class ComputeERotateSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeERotateSphereKokkos<LMPHostType>;
#endif
}
