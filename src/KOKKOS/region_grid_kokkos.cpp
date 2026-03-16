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
void RegGridKokkos<DeviceType>::match_all_kokkos(int groupbit_in, DAT::tdual_int_1d k_match_in)
{
  sync_grid_to_device();

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

  int nz_size = nzhi_out - nzlo_out + 1;
  int ny_size = nyhi_out - nylo_out + 1;
  int nx_size = nxhi_out - nxlo_out + 1;

  if (nz_size <= 0 || ny_size <= 0 || nx_size <= 0) return;

  if ((int)d_gridvals.extent(0) != nz_size ||
      (int)d_gridvals.extent(1) != ny_size ||
      (int)d_gridvals.extent(2) != nx_size) {
    d_gridvals = Kokkos::View<double***, DeviceType>(
        "RegGridKokkos:gridvals", nz_size, ny_size, nx_size);
  }

  auto h_gridvals = Kokkos::create_mirror_view(d_gridvals);

  int col = (gridindex > 0) ? gridindex - 1 : 0;

  if (ncol == 0) {
    auto vec3d = (double ***) griddata;
    for (int iz = nzlo_out; iz <= nzhi_out; iz++)
      for (int iy = nylo_out; iy <= nyhi_out; iy++)
        for (int ix = nxlo_out; ix <= nxhi_out; ix++)
          h_gridvals(iz - nzlo_out, iy - nylo_out, ix - nxlo_out) =
              vec3d[iz][iy][ix];
  } else {
    auto array3d = (double ****) griddata;
    for (int iz = nzlo_out; iz <= nzhi_out; iz++)
      for (int iy = nylo_out; iy <= nyhi_out; iy++)
        for (int ix = nxlo_out; ix <= nxhi_out; ix++)
          h_gridvals(iz - nzlo_out, iy - nylo_out, ix - nxlo_out) =
              array3d[iz][iy][ix][col];
  }

  Kokkos::deep_copy(d_gridvals, h_gridvals);

  k_boxlo0 = domain->boxlo[0];
  k_boxlo1 = domain->boxlo[1];
  k_boxlo2 = domain->boxlo[2];
  k_dxinv = nx / domain->xprd;
  k_dyinv = ny / domain->yprd;
  k_dzinv = nz / domain->zprd;
  k_dx = domain->xprd / nx;
  k_dy = domain->yprd / ny;
  k_dz = domain->zprd / nz;
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
