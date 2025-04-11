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
   Contributing author: Drew Rohskopf (SNL)
------------------------------------------------------------------------- */

#include "compute_gaussian_grid_local_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeGaussianGridLocalKokkos<DeviceType>::ComputeGaussianGridLocalKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeGaussianGridLocal(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  k_cutsq = tdual_fparams("ComputeSNAGridKokkos::cutsq",atom->ntypes+1,atom->ntypes+1);
  auto d_cutsq = k_cutsq.template view<DeviceType>();
  rnd_cutsq = d_cutsq;

  host_flag = (execution_space == Host);

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = 1; j <= atom->ntypes; j++){
      k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutsq[i][j]; //cutsq_tmp;
      k_cutsq.template modify<LMPHostType>();
    }
  }
  // Set up element lists
  int n = atom->ntypes;
  MemKK::realloc_kokkos(d_radelem,"ComputeSNAGridKokkos::radelem",n);
  MemKK::realloc_kokkos(d_sigmaelem,"ComputeSNAGridKokkos::sigmaelem",n+1);
  MemKK::realloc_kokkos(d_prefacelem,"ComputeSNAGridKokkos::prefacelem",n+1);
  MemKK::realloc_kokkos(d_argfacelem,"ComputeSNAGridKokkos::argfacelem",n+1);
  MemKK::realloc_kokkos(d_map,"ComputeSNAGridKokkos::map",n+1);
  auto h_radelem = Kokkos::create_mirror_view(d_radelem);
  auto h_sigmaelem = Kokkos::create_mirror_view(d_sigmaelem);
  auto h_prefacelem = Kokkos::create_mirror_view(d_prefacelem);
  auto h_argfacelem = Kokkos::create_mirror_view(d_argfacelem);
  auto h_map = Kokkos::create_mirror_view(d_map);
  // start from index 1 because of how compute sna/grid is
  for (int i = 1; i <= atom->ntypes; i++) {
    h_radelem(i-1) = radelem[i];
    h_sigmaelem(i-1) = sigmaelem[i];
    h_prefacelem(i-1) = prefacelem[i];
    h_argfacelem(i-1) = argfacelem[i];
  }
  Kokkos::deep_copy(d_radelem,h_radelem);
  Kokkos::deep_copy(d_sigmaelem,h_sigmaelem);
  Kokkos::deep_copy(d_prefacelem, h_prefacelem);
  Kokkos::deep_copy(d_argfacelem, h_argfacelem);
  Kokkos::deep_copy(d_map,h_map);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeGaussianGridLocalKokkos<DeviceType>::~ComputeGaussianGridLocalKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_cutsq,cutsq);
  memoryKK->destroy_kokkos(k_alocal,alocal);
  //gridlocal_allocated = 0;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeGaussianGridLocalKokkos<DeviceType>::setup()
{

  ComputeGridLocal::setup();

  // allocate arrays
  memoryKK->create_kokkos(k_alocal, alocal, size_local_rows, size_local_cols, "grid:alocal");
  array_local = alocal;
  d_alocal = k_alocal.template view<DeviceType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeGaussianGridLocalKokkos<DeviceType>::init()
{
  ComputeGaussianGridLocal::init();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeGaussianGridLocalKokkos<DeviceType>::compute_local()
{
  if (host_flag) {
    return;
  }

  invoked_local = update->ntimestep;

  copymode = 1;

  zlen = nzhi-nzlo+1;
  ylen = nyhi-nylo+1;
  xlen = nxhi-nxlo+1;
  total_range = (nzhi-nzlo+1)*(nyhi-nylo+1)*(nxhi-nxlo+1);

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  // max_neighs is defined here - think of more elaborate methods.
  max_neighs = 100;

  // Pair snap/kk uses grow_ij with some max number of neighs but compute sna/grid uses total
  // number of atoms.
  ntotal = atomKK->nlocal + atomKK->nghost;
  // Allocate view for number of neighbors per grid point
  MemKK::realloc_kokkos(d_ninside,"ComputeSNAGridKokkos:ninside",total_range);

  // "chunksize" variable is default 32768 in compute_sna_grid.cpp, and set by user
  // `total_range` is the number of grid points which may be larger than chunk size.
  // printf(">>> total_range: %d\n", total_range);
  chunksize = 32768; // 100*32768
  chunk_size = MIN(chunksize, total_range);
  chunk_offset = 0;

  int vector_length_default = 1;
  int team_size_default = 1;
  if (!host_flag)
    team_size_default = 1; // cost will increase with increasing team size //32;//max_neighs;

  if (triclinic){
    h0 = domain->h[0];
    h1 = domain->h[1];
    h2 = domain->h[2];
    h3 = domain->h[3];
    h4 = domain->h[4];
    h5 = domain->h[5];
    lo0 = domain->boxlo[0];
    lo1 = domain->boxlo[1];
    lo2 = domain->boxlo[2];
  }

  while (chunk_offset < total_range) { // chunk up loop to prevent running out of memory

    if (chunk_size > total_range - chunk_offset)
      chunk_size = total_range - chunk_offset;

    //Neigh
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagComputeGaussianGridLocalNeigh>(chunk_size,team_size,vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagComputeGaussianGridLocalNeigh> policy_neigh(chunk_size,team_size,vector_length);
      Kokkos::parallel_for("ComputeGaussianGridLocalNeigh",policy_neigh,*this);
    }

    // Proceed to the next chunk.
    chunk_offset += chunk_size;
  } // end while

  copymode = 0;

  k_alocal.template modify<DeviceType>();
  k_alocal.template sync<LMPHostType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeGaussianGridLocalKokkos<DeviceType>::operator() (TagComputeGaussianGridLocalNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagComputeGaussianGridLocalNeigh>::member_type& team) const
{
  const int ii = team.league_rank();

  if (ii >= chunk_size) return;

  // extract grid index
  int igrid = ii + chunk_offset;

  // convert to grid indices

  int iz = igrid/(xlen*ylen);
  int i2 = igrid - (iz*xlen*ylen);
  int iy = i2/xlen;
  int ix = i2 % xlen;
  iz += nzlo;
  iy += nylo;
  ix += nxlo;

  double xgrid[3];

  // index ii already captures the proper grid point
  //int igrid = iz * (nx * ny) + iy * nx + ix;

  // grid2x converts igrid to ix,iy,iz like we've done before
  // multiply grid integers by grid spacing delx, dely, delz
  //grid2x(igrid, xgrid);
  xgrid[0] = ix * delx;
  xgrid[1] = iy * dely;
  xgrid[2] = iz * delz;

  if (triclinic) {

    // Do a conversion on `xgrid` here like we do in the CPU version.

    // Can't do this:
    // domainKK->lamda2x(xgrid, xgrid);
    // Because calling a __host__ function("lamda2x") from a __host__ __device__ function("operator()") is not allowed

    // Using domainKK-> gives segfault, use domain-> instead since we're just accessing floats.
    xgrid[0] = h0*xgrid[0] + h5*xgrid[1] + h4*xgrid[2] + lo0;
    xgrid[1] = h1*xgrid[1] + h3*xgrid[2] + lo1;
    xgrid[2] = h2*xgrid[2] + lo2;
  }

  const F_FLOAT xtmp = xgrid[0];
  const F_FLOAT ytmp = xgrid[1];
  const F_FLOAT ztmp = xgrid[2];

  // Zeroing out the components, which are filled as a sum.
  for (int icol = size_local_cols_base; icol < size_local_cols; icol++){
    d_alocal(igrid, icol) = 0.0;
  }

  // Fill grid info columns
  d_alocal(igrid, 0) = ix;
  d_alocal(igrid, 1) = iy;
  d_alocal(igrid, 2) = iz;
  d_alocal(igrid, 3) = xtmp;
  d_alocal(igrid, 4) = ytmp;
  d_alocal(igrid, 5) = ztmp;

  // Looping over ntotal for now.
  for (int j = 0; j < ntotal; j++){
    const F_FLOAT dx = x(j,0) - xtmp;
    const F_FLOAT dy = x(j,1) - ytmp;
    const F_FLOAT dz = x(j,2) - ztmp;
    int jtype = type(j);
    const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;

    if (rsq < rnd_cutsq(jtype, jtype) ) {
      int icol = size_local_cols_base + jtype - 1;
      d_alocal(igrid, icol) += d_prefacelem(jtype-1) * exp(-rsq * d_argfacelem(jtype-1));
    }
  }
}

/* ----------------------------------------------------------------------
   check max team size
------------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void ComputeGaussianGridLocalKokkos<DeviceType>::check_team_size_for(int inum, int &team_size, int vector_length) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

namespace LAMMPS_NS {
template class ComputeGaussianGridLocalKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeGaussianGridLocalKokkos<LMPHostType>;
#endif
}
