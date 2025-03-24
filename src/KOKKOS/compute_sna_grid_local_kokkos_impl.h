// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing authors: Andrew Rohskopf (SNL)
------------------------------------------------------------------------- */

#include "compute_sna_grid_local_kokkos.h"
#include "pair_snap_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "domain.h"
#include "domain_kokkos.h"
#include "sna.h"
#include "update.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <iostream>

#define MAXLINE 1024
#define MAXWORD 3

namespace LAMMPS_NS {

// Constructor

template<class DeviceType, typename real_type, int vector_length>
ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::ComputeSNAGridLocalKokkos(LAMMPS *lmp, int narg, char **arg) : ComputeSNAGridLocal(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  domainKK = (DomainKokkos *) domain;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  k_cutsq = tdual_fparams("ComputeSNAGridLocalKokkos::cutsq",atom->ntypes+1,atom->ntypes+1);
  auto d_cutsq = k_cutsq.template view<DeviceType>();
  rnd_cutsq = d_cutsq;

  host_flag = (execution_space == Host);

  // TODO: Extract cutsq in double loop below, no need for cutsq_tmp

  cutsq_tmp = cutsq[1][1];

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = 1; j <= atom->ntypes; j++){
      k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutsq_tmp;
      k_cutsq.template modify<LMPHostType>();
    }
  }

   // Set up element lists
  MemKK::realloc_kokkos(d_radelem,"ComputeSNAGridLocalKokkos::radelem",nelements);
  MemKK::realloc_kokkos(d_wjelem,"ComputeSNAGridLocalKokkos:wjelem",nelements);
  MemKK::realloc_kokkos(d_sinnerelem,"ComputeSNAGridLocalKokkos:sinnerelem",nelements);
  MemKK::realloc_kokkos(d_dinnerelem,"ComputeSNAGridLocalKokkos:dinnerelem",nelements);
  // test
  MemKK::realloc_kokkos(d_test, "ComputeSNAGridLocalKokkos::test", nelements);

  int n = atom->ntypes;
  MemKK::realloc_kokkos(d_map,"ComputeSNAGridLocalKokkos::map",n+1);

  auto h_radelem = Kokkos::create_mirror_view(d_radelem);
  auto h_wjelem = Kokkos::create_mirror_view(d_wjelem);
  auto h_sinnerelem = Kokkos::create_mirror_view(d_sinnerelem);
  auto h_dinnerelem = Kokkos::create_mirror_view(d_dinnerelem);
  auto h_map = Kokkos::create_mirror_view(d_map);
  // test
  auto h_test = Kokkos::create_mirror_view(d_test);
  h_test(0) = 2.0;

  // start from index 1 because of how compute sna/grid is
  for (int i = 1; i <= atom->ntypes; i++) {
    h_radelem(i-1) = radelem[i];
    h_wjelem(i-1) = wjelem[i];
    if (switchinnerflag){
      h_sinnerelem(i) = sinnerelem[i];
      h_dinnerelem(i) = dinnerelem[i];
    }
  }

  // In pair snap some things like `map` get allocated regardless of chem flag.
  if (chemflag){
    for (int i = 1; i <= atom->ntypes; i++) {
      h_map(i) = map[i];
    }
  }

  Kokkos::deep_copy(d_radelem,h_radelem);
  Kokkos::deep_copy(d_wjelem,h_wjelem);
  if (switchinnerflag){
    Kokkos::deep_copy(d_sinnerelem,h_sinnerelem);
    Kokkos::deep_copy(d_dinnerelem,h_dinnerelem);
  }
  if (chemflag){
    Kokkos::deep_copy(d_map,h_map);
  }
  Kokkos::deep_copy(d_test,h_test);

  snaKK = SNAKokkos<DeviceType, real_type, vector_length>(*this);
  snaKK.grow_rij(0,0);
  snaKK.init();
}

// Destructor

template<class DeviceType, typename real_type, int vector_length>
ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::~ComputeSNAGridLocalKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_cutsq,cutsq);
  memoryKK->destroy_kokkos(k_alocal,alocal);
}

// Setup

template<class DeviceType, typename real_type, int vector_length>
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::setup()
{

  ComputeGridLocal::setup();

  // allocate arrays
  memoryKK->create_kokkos(k_alocal, alocal, size_local_rows, size_local_cols, "grid:alocal");
  array_local = alocal;
  d_alocal = k_alocal.template view<DeviceType>();
}

// Compute

template<class DeviceType, typename real_type, int vector_length>
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::compute_local()
{
  if (host_flag) {
    ComputeSNAGridLocal::compute_array();
    return;
  }

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
  MemKK::realloc_kokkos(d_ninside,"ComputeSNAGridLocalKokkos:ninside",total_range);

  // "chunksize" variable is default 32768 in compute_sna_grid.cpp, and set by user
  // `total_range` is the number of grid points which may be larger than chunk size.
  chunk_size = MIN(chunksize, total_range);
  chunk_offset = 0;
  //snaKK.grow_rij(chunk_size, ntotal);
  snaKK.grow_rij(chunk_size, max_neighs);

  //chunk_size = total_range;

  // Pre-compute ceil(chunk_size / vector_length) for code cleanliness
  const int chunk_size_div = (chunk_size + vector_length - 1) / vector_length;

  if (triclinic) {
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


    //ComputeNeigh
    {
      int scratch_size = scratch_size_helper<int>(team_size_compute_neigh * max_neighs); //ntotal);

      SnapAoSoATeamPolicy<DeviceType, team_size_compute_neigh, TagCSNAGridLocalComputeNeigh>
        policy_neigh(chunk_size, team_size_compute_neigh, vector_length);
      policy_neigh = policy_neigh.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
      Kokkos::parallel_for("ComputeNeigh",policy_neigh,*this);
    }

    //ComputeCayleyKlein
    {
      // tile_size_compute_ck is defined in `compute_sna_grid_kokkos.h`
      Snap3DRangePolicy<DeviceType, tile_size_compute_ck, TagCSNAGridLocalComputeCayleyKlein>
        policy_compute_ck({0,0,0}, {vector_length, max_neighs, chunk_size_div}, {vector_length, tile_size_compute_ck, 1});
      Kokkos::parallel_for("ComputeCayleyKlein", policy_compute_ck, *this);
    }

    //PreUi
    {
      auto policy_pre_ui = snap_get_policy<DeviceType, tile_size_pre_ui, TagCSNAGridLocalPreUi>(chunk_size_div, twojmax + 1);
      Kokkos::parallel_for("PreUi", policy_pre_ui, *this);
    }

    // ComputeUi w/ vector parallelism, shared memory, direct atomicAdd into ulisttot
    {
      // team_size_compute_ui is defined in `compute_sna_grid_kokkos.h`
      // scratch size: 32 atoms * (twojmax+1) cached values, no double buffer
      const int tile_size = vector_length * (twojmax + 1);
      const int scratch_size = scratch_size_helper<complex>(team_size_compute_ui * tile_size);

      if (chunk_size < parallel_thresh)
      {
        // Version with parallelism over j_bend

        // total number of teams needed: (natoms / 32) * (ntotal) * ("bend" locations)
        const int n_teams = chunk_size_div * max_neighs * (twojmax + 1);
        const int n_teams_div = (n_teams + team_size_compute_ui - 1) / team_size_compute_ui;

        SnapAoSoATeamPolicy<DeviceType, team_size_compute_ui, TagCSNAGridLocalComputeUiSmall>
          policy_ui(n_teams_div, team_size_compute_ui, vector_length);
        policy_ui = policy_ui.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
        Kokkos::parallel_for("ComputeUiSmall", policy_ui, *this);
      } else {
        // Version w/out parallelism over j_bend

        // total number of teams needed: (natoms / 32) * (ntotal)
        const int n_teams = chunk_size_div * max_neighs;
        const int n_teams_div = (n_teams + team_size_compute_ui - 1) / team_size_compute_ui;

        SnapAoSoATeamPolicy<DeviceType, team_size_compute_ui, TagCSNAGridLocalComputeUiLarge>
          policy_ui(n_teams_div, team_size_compute_ui, vector_length);
        policy_ui = policy_ui.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
        Kokkos::parallel_for("ComputeUiLarge", policy_ui, *this);
      }
    }

    //TransformUi: un-"fold" ulisttot, zero ylist
    {
      // Expand ulisttot_re,_im -> ulisttot
      // Zero out ylist
      auto policy_transform_ui = snap_get_policy<DeviceType, tile_size_transform_ui, TagCSNAGridLocalTransformUi>(chunk_size_div, snaKK.idxu_max);
      Kokkos::parallel_for("TransformUi", policy_transform_ui, *this);
    }

    //Compute bispectrum
    // team_size_[compute_zi, compute_bi, transform_bi] are defined in `pair_snap_kokkos.h`

    //ComputeZi and Bi
    if (nelements > 1) {
      auto policy_compute_zi = snap_get_policy<DeviceType, tile_size_compute_zi, TagCSNAGridLocalComputeZi<true>, min_blocks_compute_zi>(chunk_size_div, snaKK.idxz_max);
      Kokkos::parallel_for("ComputeZiChemsnap", policy_compute_zi, *this);

      auto policy_compute_bi = snap_get_policy<DeviceType, tile_size_compute_bi, TagCSNAGridLocalComputeBi<true>>(chunk_size_div, snaKK.idxb_max);
      Kokkos::parallel_for("ComputeBiChemsnap", policy_compute_bi, *this);
    } else {
      auto policy_compute_zi = snap_get_policy<DeviceType, tile_size_compute_zi, TagCSNAGridLocalComputeZi<false>, min_blocks_compute_zi>(chunk_size_div, snaKK.idxz_max);
      Kokkos::parallel_for("ComputeZi", policy_compute_zi, *this);

      auto policy_compute_bi = snap_get_policy<DeviceType, tile_size_compute_bi, TagCSNAGridLocalComputeBi<false>>(chunk_size_div, snaKK.idxb_max);
      Kokkos::parallel_for("ComputeBi", policy_compute_bi, *this);
    }

    // Fill the grid array with bispectrum values
    {
      typename Kokkos::RangePolicy<DeviceType,TagCSNAGridLocal2Fill> policy_fill(0,chunk_size);
      Kokkos::parallel_for(policy_fill, *this);
    }

    // Proceed to the next chunk.
    chunk_offset += chunk_size;

  } // end while

  copymode = 0;

  k_alocal.template modify<DeviceType>();
  k_alocal.template sync<LMPHostType>();
}

/* ----------------------------------------------------------------------
   Begin routines that are unique to the GPU codepath. These take advantage
   of AoSoA data layouts and scratch memory for recursive polynomials
------------------------------------------------------------------------- */

/*
 Simple team policy functor seeing how many layers deep we can go with the parallelism.
 */
template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeNeigh,const typename Kokkos::TeamPolicy<DeviceType,TagCSNAGridLocalComputeNeigh>::member_type& team) const {

  // This function follows similar procedure as ComputeNeigh of PairSNAPKokkos.
  // Main difference is that we don't use the neighbor class or neighbor variables here.
  // This is because the grid points are not atoms and therefore do not get assigned
  // neighbors in LAMMPS.
  // TODO: If we did make a neighborlist for each grid point, we could use current
  //       routines and avoid having to loop over all atoms (which limits us to
  //       natoms = max team size).

  // basic quantities associated with this team:
  // team_rank : rank of thread in this team
  // league_rank : rank of team in this league
  // team_size : number of threads in this team

  // extract loop index
  int ii = team.team_rank() + team.league_rank() * team.team_size();

  if (ii >= chunk_size) return;

  // extract grid index
  int igrid = ii + chunk_offset;

  // get a pointer to scratch memory
  // This is used to cache whether or not an atom is within the cutoff.
  // If it is, type_cache is assigned to the atom type.
  // If it's not, it's assigned to -1.
  //const int tile_size = ntotal; //max_neighs; // number of elements per thread
  //const int team_rank = team.team_rank();
  //const int scratch_shift = team_rank * tile_size; // offset into pointer for entire team
  //int* type_cache = (int*)team.team_shmem().get_shmem(team.team_size() * tile_size * sizeof(int), 0) + scratch_shift;

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

  // currently, all grid points are type 1
  // not clear what a better choice would be

  const int itype = 1;
  int ielem = 0;
  if (chemflag) ielem = d_map[itype];
  //const double radi = d_radelem[ielem];

  // Compute the number of neighbors, store rsq
  int ninside = 0;

  // Looping over ntotal for now.
  for (int j = 0; j < ntotal; j++){
    const F_FLOAT dx = x(j,0) - xtmp;
    const F_FLOAT dy = x(j,1) - ytmp;
    const F_FLOAT dz = x(j,2) - ztmp;
    int jtype = type(j);
    const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;

    // don't include atoms that share location with grid point
    if (rsq >= rnd_cutsq(itype,jtype) || rsq < 1e-20) {
      jtype = -1; // use -1 to signal it's outside the radius
    }

    if (jtype >= 0)
      ninside++;
  }

  d_ninside(ii) = ninside;

  // TODO: Adjust for multi-element, currently we set jelem = 0 regardless of type.
  int offset = 0;
  for (int j = 0; j < ntotal; j++){
    //const int jtype = type_cache[j];
    //if (jtype >= 0) {
    const F_FLOAT dx = x(j,0) - xtmp;
    const F_FLOAT dy = x(j,1) - ytmp;
    const F_FLOAT dz = x(j,2) - ztmp;
    const F_FLOAT rsq = dx*dx + dy*dy + dz*dz;
    int jtype = type(j);
    if (rsq < rnd_cutsq(itype,jtype) && rsq > 1e-20) {
      int jelem = 0;
      if (chemflag) jelem = d_map[jtype];
      snaKK.rij(ii,offset,0) = static_cast<real_type>(dx);
      snaKK.rij(ii,offset,1) = static_cast<real_type>(dy);
      snaKK.rij(ii,offset,2) = static_cast<real_type>(dz);
      // pair snap uses jelem here, but we use jtype, see compute_sna_grid.cpp
      // actually since the views here have values starting at 0, let's use jelem
      snaKK.wj(ii,offset) = static_cast<real_type>(d_wjelem[jelem]);
      snaKK.rcutij(ii,offset) = static_cast<real_type>((2.0 * d_radelem[jelem])*rcutfac);
      snaKK.inside(ii,offset) = j;
      if (switchinnerflag) {
        snaKK.sinnerij(ii,offset) = 0.5*(d_sinnerelem[ielem] + d_sinnerelem[jelem]);
        snaKK.dinnerij(ii,offset) = 0.5*(d_dinnerelem[ielem] + d_dinnerelem[jelem]);
      }
      if (chemflag)
        snaKK.element(ii,offset) = jelem;
      else
        snaKK.element(ii,offset) = 0;
      offset++;
    }
  }
}

/* ----------------------------------------------------------------------
  Pre-compute the Cayley-Klein parameters for reuse in later routines
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeCayleyKlein,const int iatom_mod, const int jnbor, const int iatom_div) const {

  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;

  const int ninside = d_ninside(iatom);
  if (jnbor >= ninside) return;

  snaKK.compute_cayley_klein(iatom, jnbor);
}

/* ----------------------------------------------------------------------
  Initialize the "ulisttot" structure with non-zero on-diagonal terms
  and zero terms elsewhere
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalPreUi, const int& iatom_mod, const int& j, const int& iatom_div) const {
  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;

  int itype = type(iatom);
  int ielem = d_map[itype];

  snaKK.pre_ui(iatom, j, ielem);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalPreUi, const int& iatom, const int& j) const {
  if (iatom >= chunk_size) return;

  int itype = type(iatom);
  int ielem = d_map[itype];

  snaKK.pre_ui(iatom, j, ielem);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalPreUi, const int& iatom) const {
  if (iatom >= chunk_size) return;

  const int itype = type(iatom);
  const int ielem = d_map[itype];

  for (int j = 0; j <= twojmax; j++)
    snaKK.pre_ui(iatom, j, ielem);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeUiSmall,const typename Kokkos::TeamPolicy<DeviceType,TagCSNAGridLocalComputeUiSmall>::member_type& team) const {

  // extract flattened atom_div / neighbor number / bend location
  int flattened_idx = team.team_rank() + team.league_rank() * team_size_compute_ui;

  // extract neighbor index, iatom_div
  int iatom_div = flattened_idx / (max_neighs * (twojmax + 1)); // removed "const" to work around GCC 7 bug
  const int jj_jbend = flattened_idx - iatom_div * (max_neighs * (twojmax + 1));
  const int jbend = jj_jbend / max_neighs;
  int jj = jj_jbend - jbend * max_neighs; // removed "const" to work around GCC 7 bug

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, vector_length),
    [&] (const int iatom_mod) {
    const int ii = iatom_mod + vector_length * iatom_div;
    if (ii >= chunk_size) return;

    const int ninside = d_ninside(ii);
    if (jj >= ninside) return;

    snaKK.compute_ui_small(team, iatom_mod, jbend, jj, iatom_div);
  });

}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeUiLarge,const typename Kokkos::TeamPolicy<DeviceType,TagCSNAGridLocalComputeUiLarge>::member_type& team) const {

  // extract flattened atom_div / neighbor number / bend location
  int flattened_idx = team.team_rank() + team.league_rank() * team_size_compute_ui;

  // extract neighbor index, iatom_div
  int iatom_div = flattened_idx / max_neighs; // removed "const" to work around GCC 7 bug
  int jj = flattened_idx - iatom_div * max_neighs;

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, vector_length),
    [&] (const int iatom_mod) {
    const int ii = iatom_mod + vector_length * iatom_div;
    if (ii >= chunk_size) return;

    const int ninside = d_ninside(ii);
    if (jj >= ninside) return;

    snaKK.compute_ui_large(team,iatom_mod, jj, iatom_div);
  });
}

/* ----------------------------------------------------------------------
  De-symmetrize ulisttot_re and _im and pack it into a unified ulisttot
  structure. Zero-initialize ylist. CPU and GPU.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalTransformUi, const int& iatom_mod, const int& idxu, const int& iatom_div) const {
  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;
  if (idxu >= snaKK.idxu_max) return;
  snaKK.transform_ui(iatom, idxu);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalTransformUi, const int& iatom, const int& idxu) const {
  if (iatom >= chunk_size) return;
  snaKK.transform_ui(iatom, idxu);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalTransformUi, const int& iatom) const {
  if (iatom >= chunk_size) return;
  for (int idxu = 0; idxu < snaKK.idxu_max; idxu++)
    snaKK.transform_ui(iatom, idxu);
}

/* ----------------------------------------------------------------------
  Compute all elements of the Z tensor and store them into the `zlist`
   view
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap> KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeZi<chemsnap>, const int& iatom_mod, const int& jjz, const int& iatom_div) const {
  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;
  if (jjz >= snaKK.idxz_max) return;
  snaKK.template compute_zi<chemsnap>(iatom, jjz);
}

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap> KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeZi<chemsnap>, const int& iatom, const int& jjz) const {
  if (iatom >= chunk_size) return;
  snaKK.template compute_zi<chemsnap>(iatom, jjz);
}

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap> KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeZi<chemsnap>, const int& iatom) const {
  if (iatom >= chunk_size) return;
  for (int jjz = 0; jjz < snaKK.idxz_max; jjz++)
    snaKK.template compute_zi<chemsnap>(iatom, jjz);
}

/* ----------------------------------------------------------------------
  Compute the energy triple products and store in the "blist" view
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap> KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeBi<chemsnap>, const int& iatom_mod, const int& jjb, const int& iatom_div) const {
  const int iatom = iatom_mod + iatom_div * vector_length;
  if (iatom >= chunk_size) return;
  if (jjb >= snaKK.idxb_max) return;
  snaKK.template compute_bi<chemsnap>(iatom, jjb);
}

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap> KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeBi<chemsnap>, const int& iatom, const int& jjb) const {
  if (iatom >= chunk_size) return;
  snaKK.template compute_bi<chemsnap>(iatom, jjb);
}

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap> KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocalComputeBi<chemsnap>, const int& iatom) const {
  if (iatom >= chunk_size) return;
  for (int jjb = 0; jjb < snaKK.idxb_max; jjb++)
    snaKK.template compute_bi<chemsnap>(iatom, jjb);
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::operator() (TagCSNAGridLocal2Fill, const int& ii) const {

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
  // int igrid = iz * (nx * ny) + iy * nx + ix;
  // printf("ii igrid: %d %d\n", ii, igrid);

  // grid2x converts igrid to ix,iy,iz like we've done before
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


  const auto idxb_max = snaKK.idxb_max;

  // linear contributions

  for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
    const auto idxb = icoeff % idxb_max;
    const auto idx_chem = icoeff / idxb_max;
    d_alocal(igrid,icoeff+6) = snaKK.blist(ii,idx_chem,idxb);
  }

}

/* ----------------------------------------------------------------------
   utility functions
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template<class TagStyle>
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::check_team_size_for(int inum, int &team_size) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

template<class DeviceType, typename real_type, int vector_length>
template<class TagStyle>
void ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::check_team_size_reduce(int inum, int &team_size) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelReduceTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

template<class DeviceType, typename real_type, int vector_length>
template<typename scratch_type>
int ComputeSNAGridLocalKokkos<DeviceType, real_type, vector_length>::scratch_size_helper(int values_per_team) {
  typedef Kokkos::View<scratch_type*, Kokkos::DefaultExecutionSpace::scratch_memory_space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > ScratchViewType;

  return ScratchViewType::shmem_size(values_per_team);
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   routines used by template reference classes
------------------------------------------------------------------------- */


template<class DeviceType>
ComputeSNAGridLocalKokkosDevice<DeviceType>::ComputeSNAGridLocalKokkosDevice(class LAMMPS *lmp, int narg, char **arg)
   : ComputeSNAGridLocalKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_DEVICE_VECLEN>(lmp, narg, arg) { ; }

template<class DeviceType>
void ComputeSNAGridLocalKokkosDevice<DeviceType>::compute_local()
{
  Base::compute_local();
}

#ifdef LMP_KOKKOS_GPU
template<class DeviceType>
ComputeSNAGridLocalKokkosHost<DeviceType>::ComputeSNAGridLocalKokkosHost(class LAMMPS *lmp, int narg, char **arg)
   : ComputeSNAGridLocalKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_HOST_VECLEN>(lmp, narg, arg) { ; }

template<class DeviceType>
void ComputeSNAGridLocalKokkosHost<DeviceType>::compute_local()
{
  Base::compute_local();
}
#endif

}
