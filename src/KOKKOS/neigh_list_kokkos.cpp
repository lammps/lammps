// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "neigh_list_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
NeighListKokkos<DeviceType>::NeighListKokkos(class LAMMPS *lmp):NeighList(lmp)
{
  _stride = 1;
  maxneighs = 16;
  kokkos = 1;
  maxatoms = 0;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NeighListKokkos<DeviceType>::grow(int nmax)
{
  // skip if this list is already long enough to store nmax atoms
  //  and maxneighs neighbors

  if (nmax <= maxatoms && (int)d_neighbors.extent(1) >= maxneighs) return;
  maxatoms = nmax;

  k_ilist = DAT::tdual_int_1d("neighlist:ilist",maxatoms);
  d_ilist = k_ilist.view<DeviceType>();
  d_numneigh = typename ArrayTypes<DeviceType>::t_int_1d("neighlist:numneigh",maxatoms);
  d_neighbors = typename ArrayTypes<DeviceType>::t_neighbors_2d(Kokkos::NoInit("neighlist:neighbors"),maxatoms,maxneighs);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class NeighListKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class NeighListKokkos<LMPHostType>;
#endif
}

