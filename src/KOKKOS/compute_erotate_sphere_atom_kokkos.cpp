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

#include "compute_erotate_sphere_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeErotateSphereAtomKokkos<DeviceType>::ComputeErotateSphereAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeErotateSphereAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = OMEGA_MASK | RADIUS_MASK | RMASS_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeErotateSphereAtomKokkos<DeviceType>::~ComputeErotateSphereAtomKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_erot, erot);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeErotateSphereAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow erot array if necessary

  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_erot, erot);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_erot, erot, nmax, "erotate/sphere/atom:erot");
    vector_atom = erot;
    d_erot = k_erot.template view<DeviceType>();
  }

  // compute rotational kinetic energy for each atom in group, on the device
  // point particles will have erot = 0.0, due to radius = 0.0

  atomKK->sync(execution_space, datamask_read);

  omega = atomKK->k_omega.view<DeviceType>();
  radius = atomKK->k_radius.view<DeviceType>();
  rmass = atomKK->k_rmass.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeErotateSphereAtom>(0,nlocal),*this);
  copymode = 0;

  k_erot.modify<DeviceType>();
  k_erot.sync_host();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeErotateSphereAtomKokkos<DeviceType>::operator()(TagComputeErotateSphereAtom, const int &i) const {
  if (mask[i] & groupbit) {
    d_erot[i] =
      (omega(i,0)*omega(i,0) + omega(i,1)*omega(i,1) + omega(i,2)*omega(i,2)) *
      radius[i] * radius[i] * rmass[i];
    d_erot[i] *= static_cast<KK_FLOAT>(pfactor);
  } else {
    d_erot[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeErotateSphereAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeErotateSphereAtomKokkos<LMPHostType>;
#endif
}
