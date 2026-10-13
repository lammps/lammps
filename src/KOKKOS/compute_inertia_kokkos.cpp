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

#include "compute_inertia_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "group_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeInertiaKokkos<DeviceType>::ComputeInertiaKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeInertia(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | MASK_MASK | IMAGE_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeInertiaKokkos<DeviceType>::compute_vector()
{
  invoked_vector = update->ntimestep;

  auto *groupKK = (GroupKokkos *) group;

  double xcm[3], itensor[3][3];
  if (group->dynamic[igroup]) masstotal = groupKK->mass_kk<DeviceType>(igroup);
  groupKK->xcm_kk<DeviceType>(igroup,masstotal,xcm);
  groupKK->inertia_kk<DeviceType>(igroup,xcm,itensor);

  // the extended-particle contribution reads the radius array and the
  // ellipsoid, line, tri and body bonus arrays through the plain host
  // pointers, none of which have a device counterpart.  for point particles
  // it adds nothing, so skip it entirely then and keep the whole compute on
  // the device for the usual case.  note that finite-size spheres carry only
  // a radius and no bonus data, so testing the shape atom styles alone is
  // not enough

  if ((atom->radius_flag) || (atom->ellipsoid_flag) || (atom->line_flag) ||
      (atom->tri_flag) || (atom->body_flag)) {
    atomKK->sync(Host,MASK_MASK|RMASS_MASK|TYPE_MASK|RADIUS_MASK|ELLIPSOID_MASK|BONUS_MASK);
    group->inertia_extended(igroup,itensor);
  }

  vector[0] = itensor[0][0];
  vector[1] = itensor[1][1];
  vector[2] = itensor[2][2];
  vector[3] = itensor[0][1];
  vector[4] = itensor[1][2];
  vector[5] = itensor[0][2];
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeInertiaKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeInertiaKokkos<LMPHostType>;
#endif
}
