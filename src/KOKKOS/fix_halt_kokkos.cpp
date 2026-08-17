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

#include "fix_halt_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "neighbor_kokkos.h"

#include <cmath>

using namespace LAMMPS_NS;

enum { BONDMAX, TLIMIT, DISKFREE, VARIABLE };

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixHaltKokkos<DeviceType>::FixHaltKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixHalt(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixHaltKokkos<DeviceType>::end_of_step()
{
  if (attribute == VARIABLE)
    atomKK->sync(Host, ALL_MASK);

  FixHalt::end_of_step();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double FixHaltKokkos<DeviceType>::bondmax()
{
  atomKK->sync(execution_space, X_MASK);
  neighborKK->k_bondlist.template sync<DeviceType>();

  x = atomKK->k_x.template view<DeviceType>();
  d_bondlist = neighborKK->k_bondlist.template view<DeviceType>();
  int nbondlist = neighborKK->nbondlist;

  KK_FLOAT maxone = 0.0;

  if (nbondlist > 0) {
    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixHaltBondmax>(0, nbondlist),
                            *this, Kokkos::Max<KK_FLOAT>(maxone));
    copymode = 0;
  }

  double maxall;
  MPI_Allreduce(&maxone, &maxall, 1, MPI_DOUBLE, MPI_MAX, world);

  return sqrt(maxall);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixHaltKokkos<DeviceType>::operator()(TagFixHaltBondmax, const int &n, KK_FLOAT &maxone) const
{
  const int i1 = d_bondlist(n,0);
  const int i2 = d_bondlist(n,1);
  const KK_FLOAT delx = x(i1,0) - x(i2,0);
  const KK_FLOAT dely = x(i1,1) - x(i2,1);
  const KK_FLOAT delz = x(i1,2) - x(i2,2);
  const KK_FLOAT rsq = delx*delx + dely*dely + delz*delz;
  maxone = MAX(maxone,rsq);
}

namespace LAMMPS_NS {
template class FixHaltKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixHaltKokkos<LMPHostType>;
#endif
}
