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

#include "fix_temp_csvr_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixTempCSVRKokkos<DeviceType>::FixTempCSVRKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixTempCSVR(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixTempCSVRKokkos<DeviceType>::end_of_step()
{
  // ensure velocity data is available on the host for the base class
  // the base class temperature compute, resamplekin(), and v-scaling loop
  // all run on the host using the atom->v pointer

  atomKK->sync(Host, V_MASK | MASK_MASK);

  FixTempCSVR::end_of_step();

  atomKK->modified(Host, V_MASK);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixTempCSVRKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixTempCSVRKokkos<LMPHostType>;
#endif
}
