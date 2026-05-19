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

#include "fix_press_berendsen_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "compute.h"
#include "force.h"
#include "kspace.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixPressBerendsenKokkos<DeviceType>::FixPressBerendsenKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixPressBerendsen(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPressBerendsenKokkos<DeviceType>::end_of_step()
{
  // sync atom data for temperature and pressure computes;
  // skip the sync for kokkosable computes since they handle syncing internally

  if (!temperature->kokkosable)
    atomKK->sync(temperature->execution_space, temperature->datamask_read);
  if (!pressure->kokkosable)
    atomKK->sync(pressure->execution_space, pressure->datamask_read);

  // call base class: computes T and P, couples, dilates, calls remap()
  // remap() calls domain->x2lamda(N) and domain->lamda2x(N), which are
  // overridden by DomainKokkos to operate on the device and handle syncing

  FixPressBerendsen::end_of_step();
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixPressBerendsenKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixPressBerendsenKokkos<LMPHostType>;
#endif
}
