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

#include "fix_move_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixMoveKokkos<DeviceType>::FixMoveKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixMove(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *)atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMoveKokkos<DeviceType>::initial_integrate(int vflag)
{
  // mstyle == VARIABLE requires variable evaluation (host only).
  // All motion modes update x and v, and may read f.
  // xoriginal is a plain CPU array in the base class; the entire loop runs
  // on the host.  Sync X, V, F and per-atom arrays before calling the base
  // class, then mark X and V as modified on host.

  atomKK->sync(Host, X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK |
               TYPE_MASK | OMEGA_MASK | ANGMOM_MASK);

  FixMove::initial_integrate(vflag);

  atomKK->modified(Host, X_MASK | V_MASK | OMEGA_MASK | ANGMOM_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMoveKokkos<DeviceType>::final_integrate()
{
  // final_integrate only updates v from f for non-prescribed directions

  atomKK->sync(Host, V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);

  FixMove::final_integrate();

  atomKK->modified(Host, V_MASK);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixMoveKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixMoveKokkos<LMPHostType>;
#endif
}
