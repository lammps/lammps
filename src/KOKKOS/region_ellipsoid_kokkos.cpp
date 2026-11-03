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

#include "region_ellipsoid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegEllipsoidKokkos<DeviceType>::RegEllipsoidKokkos(LAMMPS *lmp, int narg, char **arg)
  : RegEllipsoid(lmp, narg, arg)
{
  atomKK = (AtomKokkos*) atom;
}

/* ----------------------------------------------------------------------
   cache domain->dimension: a device kernel cannot dereference the pointer
------------------------------------------------------------------------- */

template<class DeviceType>
void RegEllipsoidKokkos<DeviceType>::init()
{
  RegEllipsoid::init();
  dimension = domain->dimension;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void RegEllipsoidKokkos<DeviceType>::match_all_kokkos(int groupbit_in, DAT::tdual_int_1d k_match_in)
{
  groupbit = groupbit_in;
  boxremap.capture(domain);
  d_match = k_match_in.template view<DeviceType>();
  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  atomKK->sync(execution_space, X_MASK | MASK_MASK);
  d_x = atomKK->k_x.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagRegEllipsoidMatchAll>(0,nlocal),*this);
  copymode = 0;
  k_match_in.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void RegEllipsoidKokkos<DeviceType>::operator()(TagRegEllipsoidMatchAll, const int &i) const {
  if (d_mask[i] & groupbit) {
    KK_FLOAT x_tmp = d_x(i,0);
    KK_FLOAT y_tmp = d_x(i,1);
    KK_FLOAT z_tmp = d_x(i,2);
    d_match[i] = match_kokkos(static_cast<double>(x_tmp),static_cast<double>(y_tmp),static_cast<double>(z_tmp));
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class RegEllipsoidKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class RegEllipsoidKokkos<LMPHostType>;
#endif
}

