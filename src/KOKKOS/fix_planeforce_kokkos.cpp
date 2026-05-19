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

#include "fix_planeforce_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixPlaneForceKokkos<DeviceType>::FixPlaneForceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixPlaneForce(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPlaneForceKokkos<DeviceType>::init()
{
  FixPlaneForce::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Cannot (yet) use respa with fix planeforce/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPlaneForceKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, F_MASK | MASK_MASK);

  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixPlaneForce>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space, F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixPlaneForceKokkos<DeviceType>::operator()(TagFixPlaneForce, const int &i) const
{
  if (mask[i] & groupbit) {
    const double dot = f(i,0)*xdir + f(i,1)*ydir + f(i,2)*zdir;
    f(i,0) -= dot*xdir;
    f(i,1) -= dot*ydir;
    f(i,2) -= dot*zdir;
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixPlaneForceKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixPlaneForceKokkos<LMPHostType>;
#endif
}
