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

#include "region_grid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "memory_kokkos.h"

using namespace LAMMPS_NS;
using namespace MathSpecialKokkos;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegGridKokkos<DeviceType>::RegGridKokkos(LAMMPS *lmp, int narg, char **arg)
  : RegGrid(lmp, narg, arg)
{
  atomKK = (AtomKokkos*) atom;
  memoryKK->create_kokkos(d_contact,6,"region_grid:d_contact");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
RegGridKokkos<DeviceType>::~RegGridKokkos()
{
  if (copymode) return;
  memoryKK->destroy_kokkos(d_contact);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void RegGridKokkos<DeviceType>::init()
{
  RegGrid::init();
  allocate_grid_view();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void RegGridKokkos<DeviceType>::allocate_grid_view()
{
  if (!grid3d || !griddata) return;

  RegGrid::resolve_grid_reference();
  memoryKK->create_kokkos(k_griddata,griddata,nx,ny,nz,"RegGridKokkos:k_griddata");

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void RegGridKokkos<DeviceType>::match_all_kokkos(int groupbit_in, DAT::tdual_int_1d k_match_in)
{
  k_griddata.sync_device();

  groupbit = groupbit_in;
  d_match = k_match_in.template view<DeviceType>();
  auto execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  atomKK->sync(execution_space, X_MASK | MASK_MASK);
  d_x = atomKK->k_x.view<DeviceType>();
  d_mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagRegGridMatchAll>(0,nlocal),*this);
  copymode = 0;
  k_match_in.template modify<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void RegGridKokkos<DeviceType>::sync_grid_to_device()
{
  if (!grid3d || !griddata) return;

  allocate_grid_view();

  if (d_gridvals.extent(0) == 0) return;

  int col = (gridindex > 0) ? gridindex - 1 : 0;

  k_boxlo0 = domain->boxlo[0];
  k_boxlo1 = domain->boxlo[1];
  k_boxlo2 = domain->boxlo[2];
  k_dxinv = nx / domain->xprd;
  k_dyinv = ny / domain->yprd;
  k_dzinv = nz / domain->zprd;
  k_dx = domain->xprd / nx;
  k_dy = domain->yprd / ny;
  k_dz = domain->zprd / nz;
  grid_synced = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void RegGridKokkos<DeviceType>::operator()(TagRegGridMatchAll, const int &i) const {
  if (d_mask[i] & groupbit) {
    double x_tmp = d_x(i,0);
    double y_tmp = d_x(i,1);
    double z_tmp = d_x(i,2);
    d_match[i] = match_kokkos(x_tmp,y_tmp,z_tmp);
  }
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class RegGridKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class RegGridKokkos<LMPHostType>;
#endif
}
