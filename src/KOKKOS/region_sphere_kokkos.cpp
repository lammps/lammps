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

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio at gmail)
------------------------------------------------------------------------- */

#include "region_sphere_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegSphereKokkos<DeviceType>::RegSphereKokkos(LAMMPS *lmp, int narg, char **arg)
  : RegSphere(lmp, narg, arg)
{
  atomKK = (AtomKokkos*) atom;
  memoryKK->create_kokkos(d_contact,1,"region_sphere:d_contact");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegSphereKokkos<DeviceType>::~RegSphereKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(d_contact);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void RegSphereKokkos<DeviceType>::match_all_kokkos(int groupbit_in, DAT::tdual_int_1d k_match_in)
{
  groupbit = groupbit_in;
  d_match = k_match_in.template view<DeviceType>();
  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  atomKK->sync(execution_space, X_MASK | MASK_MASK);
  d_x = atomKK->k_x.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagRegSphereMatchAll>(0,nlocal),*this);
  copymode = 0;
  k_match_in.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void RegSphereKokkos<DeviceType>::operator()(TagRegSphereMatchAll, const int &i) const {
  if (d_mask[i] & groupbit) {
    double x_tmp = d_x(i,0);
    double y_tmp = d_x(i,1);
    double z_tmp = d_x(i,2);
    d_match[i] = match_kokkos(x_tmp,y_tmp,z_tmp);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class RegSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class RegSphereKokkos<LMPHostType>;
#endif
}

