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

#include "fix_indent_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixIndentKokkos<DeviceType>::FixIndentKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixIndent(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixIndentKokkos<DeviceType>::post_force(int vflag)
{
  // variable evaluation and per-atom distance-to-indenter loop run on host
  // sync X and F before calling base class, mark F modified on host

  atomKK->sync(Host, X_MASK | F_MASK | MASK_MASK);

  FixIndent::post_force(vflag);

  atomKK->modified(Host, F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixIndentKokkos<DeviceType>::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixIndentKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixIndentKokkos<LMPHostType>;
#endif
}
