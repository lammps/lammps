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

//
// Contributing author, Richard Meng, Queen's University at Kingston, 21.01.24, contact@richardzjm.com
//

#include "pair_mtp_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType> PairMTPKokkos<DeviceType>::PairMTPKokkos(LAMMPS(*lmp)) : PairMTP(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  host_flag = (execution_space == Host);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> PairMTPKokkos<DeviceType>::~PairMTPKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_eatom, eatom);
  memoryKK->destroy_kokkos(k_vatom, vatom);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template <class DeviceType> void PairMTPKokkos<DeviceType>::init_style()
{
  if (host_flag) {
    if (lmp->kokkos->nthreads > 1)
      error->all(FLERR, "Pair style mtp/kk can currently only run on a single CPU thread.");

    PairMTP::init_style();
    return;
  }

  if (force->newton_pair == 0) error->all(FLERR, "Pair style MTP requires newton pair on.");

  // neighbor list request for KOKKOS
  neighflag = lmp->kokkos->neighflag;

  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same_v<DeviceType, LMPHostType> &&
                           !std::is_same_v<DeviceType, LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType, LMPDeviceType>);
  if (neighflag == FULL) error->all(FLERR, "Must use half neighbor list style with pair mtp/kk.");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template <class DeviceType> double PairMTPKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairMTP::init_one(i, j);
  //Don't need to do anything with the cutoff because the MTP (and original MLIP package) only uses one cutoff for all species combos.
  return cutone;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template <class DeviceType> void PairMTPKokkos<DeviceType>::coeff(int narg, char **arg)
{ PairMTP::coeff(narg, arg); }

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template <class DeviceType> void PairMTPKokkos<DeviceType>::settings(int narg, char **arg)
{
  // We may need to process in chunks to deal with memory limitations
  // For now we expect the user to specify the chunk size

  if (narg != 3 || LAMMPS_NS::utils::lowercase(arg[1]) != "chunksize")
    error->all(FLERR,
               "Pair mtp/kk requires 3 arguments {{potential_file} \"chunksize\" {chunksize}}.");

  input_chunk_size = utils::inumeric(FLERR, arg[2], true, lmp);

  PairMTP ::settings(
      1, arg);    // This also calls read_file which parses and loads the necessary arrays in host

  // Prepare and check the alpha times waves
  PairMTPKokkos::prepare_waves();

  // ---------- Now we move arrays to device ----------
  // First we set up the index lists
  MemKK::realloc_kokkos(d_alpha_index_basic, "mtp/kk:alpha_index_basic", alpha_index_basic_count,
                        4);
  MemKK::realloc_kokkos(d_alpha_index_times, "mtp/kk:alpha_index_times", alpha_index_times_count,
                        4);
  MemKK::realloc_kokkos(d_alpha_moment_mapping, "mtp/kk:moment_mapping", alpha_scalar_count);

  // Setup the learned coefficients
  int radial_coeff_count = species_count * species_count * radial_basis_size * radial_func_count;
  MemKK::realloc_kokkos(d_radial_basis_coeffs, "mtp/kk:radial_coeffs", radial_coeff_count);
  MemKK::realloc_kokkos(d_species_coeffs, "mtp/kk:species_coeffs", species_count);
  MemKK::realloc_kokkos(d_linear_coeffs, "mtp/kk:linear_coeffs", alpha_scalar_count);

  // We need to init these as very small views to begin with because the user might specify a very large chunk_size which is much more than inum. We will resize these as needed in compute.
  MemKK::realloc_kokkos(d_valid_neighs, "mtp/kk:d_valid_neighs", 1, 1);
  MemKK::realloc_kokkos(d_num_valid_neighs, "mtp/kk:d_num_valid_neighs", 1);
  MemKK::realloc_kokkos(d_moment_jacobian, "mtp/kk:moment_jacobian", 1, 1, alpha_index_basic_count,
                        3);
  MemKK::realloc_kokkos(d_moment_tensor_vals, "mtp/kk:moment_tensor_vals", 1, alpha_moment_count);
  MemKK::realloc_kokkos(d_nbh_energy_ders_wrt_moments, "mtp/kk:nbh_energy_ders_wrt_moments", 1,
                        alpha_moment_count);

  //Declare host arrays
  auto h_alpha_index_basic = Kokkos::create_mirror_view(d_alpha_index_basic);
  auto h_alpha_index_times = Kokkos::create_mirror_view(d_alpha_index_times);
  auto h_alpha_moment_mapping = Kokkos::create_mirror_view(d_alpha_moment_mapping);
  auto h_radial_basis_coeffs = Kokkos::create_mirror_view(d_radial_basis_coeffs);
  auto h_species_coeffs = Kokkos::create_mirror_view(d_species_coeffs);
  auto h_linear_coeffs = Kokkos::create_mirror_view(d_linear_coeffs);

  //Populate the host arrays
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < alpha_index_basic_count; i++)
      h_alpha_index_basic(i, j) = alpha_index_basic[i][j];
    for (int i = 0; i < alpha_index_times_count; i++)
      h_alpha_index_times(i, j) = alpha_index_times[i][j];
  }
  for (int i = 0; i < alpha_scalar_count; i++) {
    h_alpha_moment_mapping(i) = alpha_moment_mapping[i];
    h_linear_coeffs(i) = linear_coeffs[i];
  }
  for (int i = 0; i < radial_coeff_count; i++) h_radial_basis_coeffs(i) = radial_basis_coeffs[i];
  for (int i = 0; i < species_count; i++) h_species_coeffs(i) = species_coeffs[i];

  // Peform the copy from host to device
  Kokkos::deep_copy(d_alpha_index_basic, h_alpha_index_basic);
  Kokkos::deep_copy(d_alpha_index_times, h_alpha_index_times);
  Kokkos::deep_copy(d_alpha_moment_mapping, h_alpha_moment_mapping);
  Kokkos::deep_copy(d_radial_basis_coeffs, h_radial_basis_coeffs);
  Kokkos::deep_copy(d_species_coeffs, h_species_coeffs);
  Kokkos::deep_copy(d_linear_coeffs, h_linear_coeffs);
  // No need to deep copy the working buffers.
}

/* ----------------------------------------------------------------------
   Finds the size of each alpha times waves.
------------------------------------------------------------------------- */
template <class DeviceType> void PairMTPKokkos<DeviceType>::prepare_waves()
{
  // The alpha times in the MLIP-3 format are already sorted by child node.
  std::vector<int> waves;
  int last_max_node = alpha_index_basic_count - 1;
  int last_max_edge = 0;

  for (int i = 0; i < alpha_index_times_count; i++) {
    if (alpha_index_times[i][0] > last_max_node || alpha_index_times[i][1] > last_max_node) {
      waves.push_back(i - last_max_edge);
      last_max_node = alpha_index_times[i - 1][3];
      last_max_edge = i;
    }
  }
  waves.push_back(alpha_index_times_count - last_max_edge);

  MemKK::realloc_kokkos(d_waves, "mtp/kk:d_waves", waves.size());
  auto h_waves = Kokkos::create_mirror_view(d_waves);
  for (int i = 0; i < waves.size(); i++) { h_waves(i) = waves[i]; }
  Kokkos::deep_copy(d_waves, h_waves);
}

// Finds the maximum number of neighbours in all neigbhourhoods. This enables use to set the size (2nd index) of the jacobian. (Copied from other potentials)
template <class DeviceType> struct FindMaxNumNeighs {
  typedef DeviceType device_type;
  NeighListKokkos<DeviceType> k_list;

  FindMaxNumNeighs(NeighListKokkos<DeviceType> *nl) : k_list(*nl) {}
  ~FindMaxNumNeighs() { k_list.copymode = 1; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int &ii, int &max_neighs) const
  {
    const int i = k_list.d_ilist[ii];
    const int num_neighs = k_list.d_numneigh(i);
    if (max_neighs < num_neighs) max_neighs = num_neighs;
  }
};

// Finds the maximum number of valid MTP neighbours in all neigbhourhoods.
template <class DeviceType> struct FindMaxValidNeighs {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  const KK_FLOAT max_cutoff_sq;
  Kokkos::View<int *, DeviceType> d_num_valid_neighs;
  Kokkos::View<int **, DeviceType> d_valid_neighs;

  FindMaxValidNeighs(typename AT::t_int_1d_randomread d_ilist,
                     typename AT::t_int_1d_randomread d_numneigh,
                     typename AT::t_neighbors_2d d_neighbors,
                     typename AT::t_kkfloat_1d_3_lr_randomread x, KK_FLOAT max_cutoff_sq,
                     Kokkos::View<int *, DeviceType> d_num_valid_neighs,
                     Kokkos::View<int **, DeviceType> d_valid_neighs) :
      d_ilist(d_ilist), d_numneigh(d_numneigh), d_neighbors(d_neighbors), x(x),
      max_cutoff_sq(max_cutoff_sq), d_num_valid_neighs(d_num_valid_neighs),
      d_valid_neighs(d_valid_neighs)
  {
  }
  ~FindMaxValidNeighs() {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<DeviceType>::member_type &team,
                  int &max_valid_neighs) const
  {
    const int ii = team.league_rank();
    const int i = d_ilist[ii];
    const int num_neighs = d_numneigh(i);

    const KK_FLOAT xi[3] = {x(i, 0), x(i, 1), x(i, 2)};

    Kokkos::parallel_scan(Kokkos::TeamThreadRange(team, num_neighs),
                          [&](const int jj, int &prefix, const bool final) {
                            const int j = d_neighbors(i, jj) & NEIGHMASK;

                            const KK_FLOAT r0 = x(j, 0) - xi[0];
                            const KK_FLOAT r1 = x(j, 1) - xi[1];
                            const KK_FLOAT r2 = x(j, 2) - xi[2];
                            const KK_FLOAT rsq = Kokkos::fma(r0, r0, Kokkos::fma(r1, r1, r2 * r2));

                            const int is_valid = (rsq < max_cutoff_sq) ? 1 : 0;
                            const int pos = prefix;
                            prefix += is_valid;

                            if (final) {
                              if (is_valid) { d_valid_neighs(pos, ii) = j; }

                              // The last iteration’s final prefix is the total number of valid neighbors.
                              if (jj == num_neighs - 1) {
                                d_num_valid_neighs(ii) = prefix;
                                if (max_valid_neighs < prefix) max_valid_neighs = prefix;
                              }
                            }
                          });
  }
};

/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */

template <class DeviceType> void PairMTPKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  // If we are running on host we just use the base implementation
  if (host_flag) {
    atomKK->sync(Host, X_MASK | F_MASK | TYPE_MASK);
    PairMTP::compute(eflag_in, vflag_in);
    atomKK->modified(Host, F_MASK);
    return;
  }

  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag, vflag, 0);

  // reallocate per-atom arrays if necessary
  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  copymode = 1;
  int newton_pair = force->newton_pair;
  if (newton_pair == false) error->all(FLERR, "PairMTPKokkos requires 'newton on'.");

  // Now, ensure the atom data is synced
  atomKK->sync(execution_space, X_MASK | F_MASK | TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  NeighListKokkos<DeviceType> *k_list = static_cast<NeighListKokkos<DeviceType> *>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;
  inum = list->inum;

  need_dup = lmp->kokkos->need_dup<DeviceType>();
  // clang-format off
  if (need_dup) {
    dup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
    // clang-format on
  }

  //Precalc the max neighs.
  max_neighs = 0;
  Kokkos::parallel_reduce("PairMTPKokkos::find_max_neighs", inum,
                          FindMaxNumNeighs<DeviceType>(k_list), Kokkos::Max<int>(max_neighs));
  // std::cout << max_neighs << std::endl;

  if ((int) d_num_valid_neighs.extent(0) < inum) {
    Kokkos::realloc(Kokkos::WithoutInitializing, d_num_valid_neighs, inum);
  }
  if ((int) d_valid_neighs.extent(1) < inum || (int) d_valid_neighs.extent(0) < max_neighs) {
    Kokkos::realloc(Kokkos::WithoutInitializing, d_valid_neighs, max_neighs, inum);
  }

  // Precalculate the number of valid MTP neighs and stream compact them
  int max_valid_neighs = 0;
  {
    const int team_size = 64;
    Kokkos::TeamPolicy<DeviceType> policy_valid_neighs(inum, team_size);
    Kokkos::parallel_reduce("PairMTPKokkos::find_max_valid_neighs", policy_valid_neighs,
                            FindMaxValidNeighs<DeviceType>(d_ilist, d_numneigh, d_neighbors, x,
                                                           max_cutoff_sq, d_num_valid_neighs,
                                                           d_valid_neighs),
                            Kokkos::Max<int>(max_valid_neighs));
  }
  // Handling batching
  chunk_size =    // chunk_size is the working chunk size and may change per compute pass
      MIN(input_chunk_size,
          inum);    // chunksize is the maximum atoms per pass as defined by the user
  chunk_offset = 0;

  // Team sizes. We specify 32 for 1 warp per thread block.
  int team_size_default = 1;
  int vector_length_default = 1;
  if (!host_flag) team_size_default = 64;

  // Resize the arrays to the chunksize if needed. Do not initialize values, we do so in the loop.
  if ((int) d_moment_tensor_vals.extent(0) < chunk_size) {
    Kokkos::realloc(Kokkos::WithoutInitializing, d_moment_tensor_vals, chunk_size,
                    alpha_moment_count);
    Kokkos::realloc(Kokkos::WithoutInitializing, d_nbh_energy_ders_wrt_moments, chunk_size,
                    alpha_moment_count);
  }
  // Resize the jacobian and within _cutoff if max_neighs is too large. Do not initalize; first access is write.
  if ((int) d_moment_jacobian.extent(1) < chunk_size ||
      (int) d_moment_jacobian.extent(0) < max_valid_neighs) {
    Kokkos::realloc(Kokkos::WithoutInitializing, d_moment_jacobian, max_valid_neighs, chunk_size,
                    alpha_index_basic_count, 3);
  }

  EV_FLOAT ev;

  // ========== Begin Main Computation ==========
  while (chunk_offset < inum) {    // batching to prevent OOM on device
    EV_FLOAT ev_tmp;
    if (chunk_size > inum - chunk_offset) chunk_size = inum - chunk_offset;
    // ========== Init working views as 0  ==========
    {
      typename Kokkos::MDRangePolicy<Kokkos::Rank<2>, DeviceType, TagPairMTPInitMomentValsDers>
          policy_moment_init({0, 0}, {alpha_moment_count, chunk_size});
      Kokkos::parallel_for("InitMomentValDers", policy_moment_init, *this);
    }

    // ========== Calculate the basic alphas (Per outer-atom parallelizaton) ==========
    {
      int team_size = team_size_default;
      if (!host_flag && max_valid_neighs < 32) team_size = 32;
      int vector_length = vector_length_default;
      check_team_size_for<TagPairMTPComputeAlphaBasic>(chunk_size * max_valid_neighs, team_size,
                                                       vector_length);
      int radial_scratch_count = 2 * (radial_func_count + radial_basis_size);
      int dist_coords_scratch_count = 4 * max_alpha_index_basic;

      // Reduce the scratch size to the max number of neighbors
      int scratch_size = scratch_size_helper<KK_FLOAT>(
          min(team_size, max_valid_neighs) * (radial_scratch_count + dist_coords_scratch_count));
      Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaBasic> policy_basic_alpha(chunk_size,
                                                                                     team_size);
      policy_basic_alpha = policy_basic_alpha.set_scratch_size(0, Kokkos::PerTeam(scratch_size));
      Kokkos::parallel_for("ComputeAlphaBasic", policy_basic_alpha, *this);
    }

    // ========== Calculate the non-elementary alphas  ==========
    {
      int team_size = team_size_default;
      // Best team size depends on the max number of blocks per SM. 64 is good for CC8, and CC > 9+.
      Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaTimes> policy_basic_alpha(chunk_size,
                                                                                     team_size);
      Kokkos::parallel_for("ComputeAlphaTimes", policy_basic_alpha, *this);
    }

    // ========== Set the scalar nbh ders wrt moments ==========
    {
      typename Kokkos::MDRangePolicy<Kokkos::Rank<2>, DeviceType, TagPairMTPSetScalarNbhDers>
          policy_nbh_init({0, 0}, {alpha_scalar_count, chunk_size});
      Kokkos::parallel_for("SetScalarNbhDers", policy_nbh_init, *this);
    }

    // ========== Calc the nbh ders wrt moments ==========
    {
      int team_size = team_size_default;
      // Best team size depends on the max number of blocks per SM. 64 is good for CC8, and CC > 9+.
      Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeNbhDers> policy_basic_alpha(chunk_size,
                                                                                  team_size);
      Kokkos::parallel_for("ComputeNbhDers", policy_basic_alpha, *this);
    }

    // ========== Compute force (and dot product with alphas to get energy if needed) ==========
    {
      int team_size = team_size_default;
      if (!host_flag && max_neighs < 32) team_size = 32;

      if (neighflag == HALF) {
        Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeForce<HALF, 1>> policy_force(chunk_size,
                                                                                     team_size);
        Kokkos::parallel_reduce(policy_force, *this, ev_tmp);
      } else if (neighflag == HALFTHREAD) {
        Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeForce<HALFTHREAD, 1>> policy_force(
            chunk_size, team_size);
        Kokkos::parallel_reduce(policy_force, *this, ev_tmp);
      }
    }

    ev += ev_tmp;
    chunk_offset += chunk_size;    // Manage halt condition
  }    // end batching while loop

  // ========== End Main Computation ==========

  if (need_dup) Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup) Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  atomKK->modified(execution_space, F_MASK);

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f = decltype(dup_f)();
    dup_vatom = decltype(dup_vatom)();
  }
}

// ========== Kernels ==========

// Inits the working arrays: moment and ders, jacobian not needed.
template <class DeviceType>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(TagPairMTPInitMomentValsDers,
                                                                  const int &k, const int &ii) const
{
  d_moment_tensor_vals(ii, k) = 0;
  d_nbh_energy_ders_wrt_moments(ii, k) = 0;
}

// Calculates the basic alphas using fused operations where possible
template <class DeviceType>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(
    TagPairMTPComputeAlphaBasic,
    const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaBasic>::member_type &team)
    const
{
  // Extract the atom number
  int ii = team.league_rank();
  int thread = team.team_rank();

  // Get information about the central atom
  const int i = d_ilist[ii + chunk_offset];
  const KK_FLOAT xi[3] = {x(i, 0), x(i, 1), x(i, 2)};
  const int itype = type[i] - 1;    // switch to zero indexing
  const int jnum = d_num_valid_neighs(ii + chunk_offset);
  const int array_size = Kokkos::min(team.team_size(), jnum);

  shared_double_2d s_radial_vals(team.team_scratch(0), array_size, radial_func_count);
  shared_double_2d s_radial_ders(team.team_scratch(0), array_size, radial_func_count);
  shared_double_2d s_dist_powers(team.team_scratch(0), array_size, max_alpha_index_basic);
  shared_double_3d s_coord_powers(team.team_scratch(0), array_size, max_alpha_index_basic);
  shared_double_2d s_radial_basis_vals(team.team_scratch(0), array_size, radial_basis_size);
  shared_double_2d s_radial_basis_ders(team.team_scratch(0), array_size, radial_basis_size);

  // Now we calculate the alpha basics.
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, jnum), [&](const int jj) {
    const int j = d_valid_neighs(jj, ii + chunk_offset);
    const int jtype = type[j] - 1;    // switch to zero indexing
    const KK_FLOAT r[3] = {x(j, 0) - xi[0], x(j, 1) - xi[1], x(j, 2) - xi[2]};
    const KK_FLOAT rsq = Kokkos::fma(r[0], r[0], Kokkos::fma(r[1], r[1], r[2] * r[2]));
    const KK_FLOAT dist = sqrt(rsq);

    s_dist_powers(thread, 0) = s_coord_powers(thread, 0, 0) = s_coord_powers(thread, 0, 1) =
        s_coord_powers(thread, 0, 2) = 1;    // Set the constants

    // Precompute the coord and distance power
    for (int k = 1; k < max_alpha_index_basic; k++) {
      s_dist_powers(thread, k) = s_dist_powers(thread, k - 1) * dist;
      for (int a = 0; a < 3; a++)
        s_coord_powers(thread, k, a) = s_coord_powers(thread, k - 1, a) * r[a];
    }

    // Calculate the radial basis and store in shared memory
    KK_FLOAT mult = 2.0 / (max_cutoff - min_cutoff);
    KK_FLOAT ksi = Kokkos::fma(2.0, dist, -(min_cutoff + max_cutoff)) / (max_cutoff - min_cutoff);

    KK_FLOAT temp = dist - max_cutoff;
    s_radial_basis_vals(thread, 0) = scaling * temp * temp;
    s_radial_basis_vals(thread, 1) = scaling * ksi * temp * temp;
    for (int k = 2; k < radial_basis_size; k++) {
      s_radial_basis_vals(thread, k) = Kokkos::fma(2.0 * ksi, s_radial_basis_vals(thread, k - 1),
                                                   -s_radial_basis_vals(thread, k - 2));
    }

    // Do the same with the derivatives
    s_radial_basis_ders(thread, 0) = scaling * 2.0 * temp;
    s_radial_basis_ders(thread, 1) = scaling * Kokkos::fma(mult, temp * temp, 2.0 * ksi * temp);
    for (int k = 2; k < radial_basis_size; k++) {
      KK_FLOAT tmp = Kokkos::fma(mult, s_radial_basis_vals(thread, k - 1),
                                 ksi * s_radial_basis_ders(thread, k - 1));
      s_radial_basis_ders(thread, k) = Kokkos::fma(2.0, tmp, -s_radial_basis_ders(thread, k - 2));
    }

    // Precompute the mu vals and ders
    int pair_offset = itype * species_count + jtype;
    for (int mu = 0; mu < radial_func_count; mu++) {
      KK_FLOAT val = 0;
      KK_FLOAT der = 0;
      int offset = (pair_offset * radial_basis_size * radial_func_count) + mu * radial_basis_size;

      for (int ri = 0; ri < radial_basis_size; ri++) {
        val = Kokkos::fma(d_radial_basis_coeffs(offset + ri), s_radial_basis_vals(thread, ri), val);
        der = Kokkos::fma(d_radial_basis_coeffs(offset + ri), s_radial_basis_ders(thread, ri), der);
      }

      s_radial_vals(thread, mu) = val;
      s_radial_ders(thread, mu) = der;
    }

    //Now, we loop through all the basic alphas
    for (int k = 0; k < alpha_index_basic_count; k++) {

      int mu = d_alpha_index_basic(k, 0);
      int a0 = d_alpha_index_basic(k, 1);
      int a1 = d_alpha_index_basic(k, 2);
      int a2 = d_alpha_index_basic(k, 3);

      KK_FLOAT val = s_radial_vals(thread, mu);
      KK_FLOAT der = s_radial_ders(thread, mu);

      // Normalize by the rank of alpha's coresponding tensor
      int norm_rank = a0 + a1 + a2;
      KK_FLOAT norm_fac = 1.0 / s_dist_powers(thread, norm_rank);
      val *= norm_fac;
      der = Kokkos::fma(norm_fac, der, -norm_rank * val / dist);

      KK_FLOAT pow0 = s_coord_powers(thread, a0, 0);
      KK_FLOAT pow1 = s_coord_powers(thread, a1, 1);
      KK_FLOAT pow2 = s_coord_powers(thread, a2, 2);
      KK_FLOAT pow = pow0 * pow1 * pow2;
      Kokkos::atomic_add(&d_moment_tensor_vals(ii, k), val * pow);

      // Get the component's derivatives too
      KK_FLOAT temp_jac[3] = {pow * r[0], pow * r[1], pow * r[2]};

      pow *= der / dist;
      temp_jac[0] = pow * r[0];
      temp_jac[1] = pow * r[1];
      temp_jac[2] = pow * r[2];

      if (a0 != 0)
        temp_jac[0] =
            Kokkos::fma(val * a0, s_coord_powers(thread, a0 - 1, 0) * pow1 * pow2, temp_jac[0]);
      if (a1 != 0)
        temp_jac[1] =
            Kokkos::fma(val * a1, pow0 * s_coord_powers(thread, a1 - 1, 1) * pow2, temp_jac[1]);
      if (a2 != 0)
        temp_jac[2] =
            Kokkos::fma(val * a2, pow0 * pow1 * s_coord_powers(thread, a2 - 1, 2), temp_jac[2]);

      d_moment_jacobian(jj, ii, k, 0) = temp_jac[0];
      d_moment_jacobian(jj, ii, k, 1) = temp_jac[1];
      d_moment_jacobian(jj, ii, k, 2) = temp_jac[2];
    }
  });
}

// Calculates the non-elementary alpha from the basic alphas
template <class DeviceType>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(
    TagPairMTPComputeAlphaTimes,
    const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaTimes>::member_type &team)
    const
{
  const int ii = team.league_rank();

  int offset = 0;
  // Traverse all edges in the alpha times compute graph. We need to do this in waves to ensure dependencies.
  for (int i = 0; i < 3; i++) {
    const int wave_size = d_waves[i];
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, wave_size), [&](const int kk) {
      const int k = offset + kk;    // Offset for the wave
      const int a0 = d_alpha_index_times(k, 0);
      const int a1 = d_alpha_index_times(k, 1);
      const int mult = d_alpha_index_times(k, 2);
      const int a3 = d_alpha_index_times(k, 3);

      const KK_FLOAT val0 = d_moment_tensor_vals(ii, a0);
      const KK_FLOAT val1 = d_moment_tensor_vals(ii, a1);

      Kokkos::atomic_add(&d_moment_tensor_vals(ii, a3), mult * val0 * val1);
    });
    offset += wave_size;
    team.team_barrier();    // Wait for the wave to finish
  }
}

// Sets the nbh energy ders as the linear coeffs
template <class DeviceType>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(TagPairMTPSetScalarNbhDers,
                                                                  const int &k, const int &ii) const
{ d_nbh_energy_ders_wrt_moments(ii, d_alpha_moment_mapping(k)) = d_linear_coeffs(k); }

template <class DeviceType>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(
    TagPairMTPComputeNbhDers,
    const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeNbhDers>::member_type &team)
    const
{

  const int ii = team.league_rank();

  int offset = alpha_index_times_count;
  // Traverse all edges in the alpha times compute graph. We need to do this in reverse waves to ensure dependencies.
  for (int i = 2; i >= 0; i--) {
    const int wave_size = d_waves[i];
    offset -= wave_size;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, wave_size), [&](const int kk) {
      const int k = kk + offset;    // Offset for the wave
      const int a0 = d_alpha_index_times(k, 0);
      const int a1 = d_alpha_index_times(k, 1);
      const int mult = d_alpha_index_times(k, 2);
      const int a3 = d_alpha_index_times(k, 3);

      const KK_FLOAT val0 = d_moment_tensor_vals(ii, a0);
      const KK_FLOAT val1 = d_moment_tensor_vals(ii, a1);
      const KK_FLOAT val3 = d_nbh_energy_ders_wrt_moments(ii, a3);

      Kokkos::atomic_add(&d_nbh_energy_ders_wrt_moments(ii, a1), val3 * mult * val0);
      Kokkos::atomic_add(&d_nbh_energy_ders_wrt_moments(ii, a0), val3 * mult * val1);
    });
    team.team_barrier();    // Wait for the wave to finish
  }
}

// Computes forces from jac and nbh ders
template <class DeviceType>
template <int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(
    const TagPairMTPComputeForce<NEIGHFLAG, EVFLAG> &,    // Tag parameter added here
    const typename Kokkos::TeamPolicy<DeviceType,
                                      TagPairMTPComputeForce<NEIGHFLAG, EVFLAG>>::member_type &team,
    EV_FLOAT &ev) const
{
  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial
  auto v_f =
      ScatterViewHelper<NeedDup_v<NEIGHFLAG, DeviceType>, decltype(dup_f), decltype(ndup_f)>::get(
          dup_f, ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG, DeviceType>>();

  const int ii = team.league_rank();
  const int i = d_ilist[ii + chunk_offset];
  const int jnum = d_num_valid_neighs(ii + chunk_offset);
  bool need_energies = EVFLAG && eflag_either;

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, jnum), [&](const int jj) {
    const int j = d_valid_neighs(jj, ii + chunk_offset);
    KK_ACC_FLOAT temp_force[3] = {0, 0, 0};
    for (int k = 0; k < alpha_index_basic_count; k++) {
      for (int a = 0; a < 3; a++) {
        //Calculate forces
        temp_force[a] += d_nbh_energy_ders_wrt_moments(ii, k) * d_moment_jacobian(jj, ii, k, a);
      }
    }

    a_f(i, 0) += temp_force[0];
    a_f(i, 1) += temp_force[1];
    a_f(i, 2) += temp_force[2];

    a_f(j, 0) -= temp_force[0];
    a_f(j, 1) -= temp_force[1];
    a_f(j, 2) -= temp_force[2];

    if (need_energies) {
      KK_FLOAT r[3] = {x(j, 0) - x(j, 0), x(j, 1) - x(j, 1), x(j, 2) - x(j, 2)};
      v_tally_xyz<NEIGHFLAG>(ev, i, j, temp_force[0], temp_force[1], temp_force[2], r[0], r[1],
                             r[2]);
    }
  });

  if (need_energies) {
    const int itype = type(i) - 1;    // zero indexing
    KK_FLOAT nbh_energy = 0;

    // Reduction to find the dot product of the linear coeffs and the moment tensor vals
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, alpha_scalar_count),
        [&](const int k, KK_FLOAT &sum) {
          sum += d_linear_coeffs(k) * d_moment_tensor_vals(ii, d_alpha_moment_mapping(k));
        },
        nbh_energy);

    // A single team member updates the global array
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      nbh_energy += d_species_coeffs[itype];    // Essentially the reference energy
      if (eflag_global) ev.evdwl += nbh_energy;
      if (eflag_atom) d_eatom[i] = nbh_energy;
    });
  }
}

// =========== Helper Functions (Also used in other Kokkos potentials)===========
template <class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION void
PairMTPKokkos<DeviceType>::v_tally_xyz(EV_FLOAT &ev, const int &i, const int &j, const KK_FLOAT &fx,
                                       const KK_FLOAT &fy, const KK_FLOAT &fz, const KK_FLOAT &delx,
                                       const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  // The vatom array is duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG, DeviceType>, decltype(dup_vatom),
                                   decltype(ndup_vatom)>::get(dup_vatom, ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG, DeviceType>>();

  const KK_FLOAT v0 = delx * fx;
  const KK_FLOAT v1 = dely * fy;
  const KK_FLOAT v2 = delz * fz;
  const KK_FLOAT v3 = delx * fy;
  const KK_FLOAT v4 = delx * fz;
  const KK_FLOAT v5 = dely * fz;

  if (vflag_global) {
    ev.v[0] += v0;
    ev.v[1] += v1;
    ev.v[2] += v2;
    ev.v[3] += v3;
    ev.v[4] += v4;
    ev.v[5] += v5;
  }

  if (vflag_atom) {
    a_vatom(i, 0) += 0.5 * v0;
    a_vatom(i, 1) += 0.5 * v1;
    a_vatom(i, 2) += 0.5 * v2;
    a_vatom(i, 3) += 0.5 * v3;
    a_vatom(i, 4) += 0.5 * v4;
    a_vatom(i, 5) += 0.5 * v5;
    a_vatom(j, 0) += 0.5 * v0;
    a_vatom(j, 1) += 0.5 * v1;
    a_vatom(j, 2) += 0.5 * v2;
    a_vatom(j, 3) += 0.5 * v3;
    a_vatom(j, 4) += 0.5 * v4;
    a_vatom(j, 5) += 0.5 * v5;
  }
}

template <class DeviceType>
template <class TagStyle>
void PairMTPKokkos<DeviceType>::check_team_size_for(int inum, int &team_size, int vector_length)
{
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType, TagStyle>(inum, Kokkos::AUTO)
                      .team_size_max(*this, Kokkos::ParallelForTag());

  if (team_size * vector_length > team_size_max) team_size = team_size_max / vector_length;
}

template <class DeviceType>
template <typename scratch_type>
int PairMTPKokkos<DeviceType>::scratch_size_helper(int values_per_team)
{
  typedef Kokkos::View<scratch_type *, Kokkos::DefaultExecutionSpace::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      ScratchViewType;

  return ScratchViewType::shmem_size(values_per_team);
}    // namespace LAMMPS_NS

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairMTPKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairMTPKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS