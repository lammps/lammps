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

#include "fix_press_langevin_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "compute.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixPressLangevinKokkos<DeviceType>::FixPressLangevinKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixPressLangevin(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPressLangevinKokkos<DeviceType>::post_force(int vflag)
{
  // sync atom data for the pressure compute before it runs

  atomKK->sync(pressure->execution_space, pressure->datamask_read);

  FixPressLangevin::post_force(vflag);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPressLangevinKokkos<DeviceType>::pre_exchange()
{
  if (!pre_exchange_flag) return;

  // pre_exchange reads and writes x and image on host via direct array access

  atomKK->sync(Host, X_MASK | IMAGE_MASK);

  FixPressLangevin::pre_exchange();

  // after remap, DomainKokkos::lamda2x marks x as modified on device;
  // image was modified on host by domain->image_flip + domain->remap
  atomKK->modified(Host, IMAGE_MASK);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixPressLangevinKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixPressLangevinKokkos<LMPHostType>;
#endif
}
