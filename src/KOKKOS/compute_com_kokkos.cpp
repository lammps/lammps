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

#include "compute_com_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "group_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeCOMKokkos<DeviceType>::ComputeCOMKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeCOM(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | MASK_MASK | IMAGE_MASK | RMASS_MASK | TYPE_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeCOMKokkos<DeviceType>::compute_vector()
{
  invoked_vector = update->ntimestep;

  auto *groupKK = (GroupKokkos *) group;
  if (group->dynamic[igroup]) masstotal = groupKK->mass_kk<DeviceType>(igroup);

  groupKK->xcm_kk<DeviceType>(igroup,masstotal,vector);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class ComputeCOMKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeCOMKokkos<LMPHostType>;
#endif
}
