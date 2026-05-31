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

#include "compute_ke_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "force.h"
#include "memory_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeKEAtomKokkos<DeviceType>::ComputeKEAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeKEAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = V_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeKEAtomKokkos<DeviceType>::~ComputeKEAtomKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(k_ke, ke);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeKEAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow ke array if necessary

  if (atom->nmax > nmax) {
    memoryKK->destroy_kokkos(k_ke, ke);
    nmax = atom->nmax;
    memoryKK->create_kokkos(k_ke, ke, nmax, "ke/atom:ke");
    vector_atom = ke;
    d_ke = k_ke.template view<DeviceType>();
  }

  // compute kinetic energy for each atom in group, on the device

  atomKK->sync(execution_space, datamask_read);
  atomKK->k_mass.sync<DeviceType>();

  mvv2e = force->mvv2e;
  v = atomKK->k_v.view<DeviceType>();
  if (atomKK->rmass)
    rmass = atomKK->k_rmass.view<DeviceType>();
  else
    mass = atomKK->k_mass.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  if (atomKK->rmass)
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeKEAtom<1> >(0,nlocal),*this);
  else
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagComputeKEAtom<0> >(0,nlocal),*this);
  copymode = 0;

  k_ke.modify<DeviceType>();
  k_ke.sync_host();
}

template<class DeviceType>
template<int RMASS>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void ComputeKEAtomKokkos<DeviceType>::operator()(TagComputeKEAtom<RMASS>, const int &i) const {
  if (mask[i] & groupbit) {
    KK_FLOAT massone = 0.0;
    if (RMASS) massone = rmass[i];
    else massone = mass[type[i]];
    d_ke[i] = 0.5 * mvv2e * massone *
      (v(i,0)*v(i,0) + v(i,1)*v(i,1) + v(i,2)*v(i,2));
  } else {
    d_ke[i] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeKEAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeKEAtomKokkos<LMPHostType>;
#endif
}
