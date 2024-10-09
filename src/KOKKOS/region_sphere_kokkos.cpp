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

namespace LAMMPS_NS {
template class RegSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class RegSphereKokkos<LMPHostType>;
#endif
}

