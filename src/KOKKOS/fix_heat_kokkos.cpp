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

#include "fix_heat_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixHeatKokkos<DeviceType>::FixHeatKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixHeat(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixHeatKokkos<DeviceType>::end_of_step()
{
  // group->ke(), group->vcm(), and the velocity-scaling loop all run on the
  // host.  Sync V (and X for region matching) to host before calling the base
  // class, then mark V as modified on host.

  atomKK->sync(Host, V_MASK | X_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  FixHeat::end_of_step();

  atomKK->modified(Host, V_MASK);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixHeatKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixHeatKokkos<LMPHostType>;
#endif
}
