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

#include "fix_lineforce_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixLineForceKokkos<DeviceType>::FixLineForceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixLineForce(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixLineForceKokkos<DeviceType>::init()
{
  FixLineForce::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Cannot (yet) use respa with fix lineforce/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixLineForceKokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync(execution_space, F_MASK | MASK_MASK);

  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixLineForce>(0,nlocal),*this);
  copymode = 0;

  atomKK->modified(execution_space, F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixLineForceKokkos<DeviceType>::operator()(TagFixLineForce, const int &i) const
{
  if (mask[i] & groupbit) {
    const KK_FLOAT xdir_kk = static_cast<KK_FLOAT>(xdir);
    const KK_FLOAT ydir_kk = static_cast<KK_FLOAT>(ydir);
    const KK_FLOAT zdir_kk = static_cast<KK_FLOAT>(zdir);
    const KK_FLOAT dot = static_cast<KK_FLOAT>(f(i,0))*xdir_kk + static_cast<KK_FLOAT>(f(i,1))*ydir_kk + static_cast<KK_FLOAT>(f(i,2))*zdir_kk;
    f(i,0) = static_cast<KK_ACC_FLOAT>(dot*xdir_kk);
    f(i,1) = static_cast<KK_ACC_FLOAT>(dot*ydir_kk);
    f(i,2) = static_cast<KK_ACC_FLOAT>(dot*zdir_kk);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixLineForceKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixLineForceKokkos<LMPHostType>;
#endif
}
