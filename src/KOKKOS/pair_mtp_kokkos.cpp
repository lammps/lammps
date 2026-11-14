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
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"

#include <algorithm>
#include <vector>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template <class DeviceType> PairMTPKokkos<DeviceType>::PairMTPKokkos(LAMMPS *lmp) : PairMTP(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  host_flag = (execution_space == Host);

  chunk_size = 0;
  chunk_offset = 0;
  inum = 0;
  max_valid_neighs = 0;
  num_waves = 0;
  wave_begin = wave_end = 0;
  node_partitions = 1;

  inv_cutoff_range = 0.0;
  cutoff_sum = 0.0;
  radial_mult = 0.0;

  ts_basic = ts_times = ts_nbh = ts_nbh_long = 0;
  ts_force[0][0] = ts_force[0][1] = ts_force[1][0] = ts_force[1][1] = 0;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> PairMTPKokkos<DeviceType>::~PairMTPKokkos()
{
  if (copymode) return;
  if (host_flag) return;

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
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

template <class DeviceType> void PairMTPKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairMTP::coeff(narg, arg);
  PairMTPKokkos::prepare_waves();

  // The occupancy probes below depend on the model dimensions, so a second
  // pair_coeff must re-run them instead of reusing the previous model's sizes.
  ts_basic = ts_times = ts_nbh = ts_nbh_long = 0;
  ts_force[0][0] = ts_force[0][1] = ts_force[1][0] = ts_force[1][1] = 0;

  inv_cutoff_range = 1.0 / (max_cutoff - min_cutoff);
  cutoff_sum = min_cutoff + max_cutoff;
  radial_mult = 2.0 / (max_cutoff - min_cutoff);

  // ---------- Now we move arrays to device ----------
  // First we set up the index lists
  MemKK::realloc_kokkos(d_alpha_index_basic, "mtp/kk:alpha_index_basic", alpha_index_basic_count,
                        4);
  MemKK::realloc_kokkos(d_alpha_index_times, "mtp/kk:alpha_index_times", alpha_index_times_count,
                        4);
  MemKK::realloc_kokkos(d_alpha_moment_mapping, "mtp/kk:moment_mapping", alpha_scalar_count);
  MemKK::realloc_kokkos(d_map, "mtp/kk:mapping", atom->ntypes + 1);

  // Setup the learned coefficients
  MemKK::realloc_kokkos(d_radial_basis_coeffs, "mtp/kk:radial_coeffs", radial_coeff_count);
  MemKK::realloc_kokkos(d_species_coeffs, "mtp/kk:species_coeffs", species_count);
  MemKK::realloc_kokkos(d_linear_coeffs, "mtp/kk:linear_coeffs", alpha_scalar_count);
  MemKK::realloc_kokkos(d_moment_coeffs, "mtp/kk:moment_coeffs", alpha_moment_count);

  // We will grow these as needed in compute.
  MemKK::realloc_kokkos(d_valid_neighs, "mtp/kk:d_valid_neighs", 1, 1);
  MemKK::realloc_kokkos(d_num_valid_neighs, "mtp/kk:d_num_valid_neighs", 1);
  MemKK::realloc_kokkos(d_radial_vals, "mtp/kk:cached_radial_vals", 1, radial_func_count, 1);
  MemKK::realloc_kokkos(d_radial_ders, "mtp/kk:cached_radial_ders", 1, radial_func_count, 1);
  MemKK::realloc_kokkos(d_inv_dist, "mtp/kk:cached_inv_dist", 1, 1);
  MemKK::realloc_kokkos(d_moment_tensor_vals, "mtp/kk:moment_tensor_vals", 1, alpha_moment_count);
  MemKK::realloc_kokkos(d_nbh_energy_ders_wrt_moments, "mtp/kk:nbh_energy_ders_wrt_moments", 1,
                        alpha_moment_count);

  //Declare host arrays
  auto h_alpha_index_basic = Kokkos::create_mirror_view(d_alpha_index_basic);
  auto h_alpha_index_times = Kokkos::create_mirror_view(d_alpha_index_times);
  auto h_alpha_moment_mapping = Kokkos::create_mirror_view(d_alpha_moment_mapping);
  auto h_map = Kokkos::create_mirror_view(d_map);
  auto h_radial_basis_coeffs = Kokkos::create_mirror_view(d_radial_basis_coeffs);
  auto h_species_coeffs = Kokkos::create_mirror_view(d_species_coeffs);
  auto h_linear_coeffs = Kokkos::create_mirror_view(d_linear_coeffs);
  auto h_moment_coeffs = Kokkos::create_mirror_view(d_moment_coeffs);

  //Populate the host arrays
  for (int j = 0; j < 4; j++) {
    for (int i = 0; i < alpha_index_basic_count; i++)
      h_alpha_index_basic(i, j) = alpha_index_basic[i][j];
    for (int i = 0; i < alpha_index_times_count; i++)
      h_alpha_index_times(i, j) = alpha_index_times[i][j];
  }
  for (int i = 0; i < alpha_moment_count; i++) h_moment_coeffs(i) = 0;
  for (int i = 0; i < alpha_scalar_count; i++) {
    h_alpha_moment_mapping(i) = alpha_moment_mapping[i];
    h_linear_coeffs(i) = linear_coeffs[i];
    h_moment_coeffs(alpha_moment_mapping[i]) += linear_coeffs[i];
  }

  h_map[0] = -1;
  for (int i = 1; i < atom->ntypes + 1; i++) h_map[i] = map[i];
  for (int i = 0; i < radial_coeff_count; i++) h_radial_basis_coeffs(i) = radial_basis_coeffs[i];
  for (int i = 0; i < species_count; i++) h_species_coeffs(i) = species_coeffs[i];

  // Peform the copy from host to device
  Kokkos::deep_copy(d_alpha_index_basic, h_alpha_index_basic);
  Kokkos::deep_copy(d_alpha_index_times, h_alpha_index_times);
  Kokkos::deep_copy(d_alpha_moment_mapping, h_alpha_moment_mapping);
  Kokkos::deep_copy(d_map, h_map);
  Kokkos::deep_copy(d_radial_basis_coeffs, h_radial_basis_coeffs);
  Kokkos::deep_copy(d_species_coeffs, h_species_coeffs);
  Kokkos::deep_copy(d_linear_coeffs, h_linear_coeffs);
  Kokkos::deep_copy(d_moment_coeffs, h_moment_coeffs);
}

/* ----------------------------------------------------------------------
   Groups rules and nodes by dependency.
------------------------------------------------------------------------- */
template <class DeviceType> void PairMTPKokkos<DeviceType>::prepare_waves()
{
  const int n = alpha_moment_count;
  std::vector<int> forward_offsets(n + 1, 0), reverse_offsets(n + 1, 0);
  for (int k = 0; k < alpha_index_times_count; k++) {
    const int a0 = alpha_index_times[k][0];
    const int a1 = alpha_index_times[k][1];
    const int a3 = alpha_index_times[k][3];
    if (a0 < 0 || a0 >= n || a1 < 0 || a1 >= n || a3 < alpha_index_basic_count || a3 >= n)
      error->all(FLERR, "Invalid MTP contraction indices.");
    forward_offsets[a3 + 1]++;
    reverse_offsets[a0 + 1]++;
    reverse_offsets[a1 + 1]++;
  }
  for (int i = 0; i < n; i++) {
    forward_offsets[i + 1] += forward_offsets[i];
    reverse_offsets[i + 1] += reverse_offsets[i];
  }

  std::vector<int> forward_rules(alpha_index_times_count);
  std::vector<int> reverse_rules(2 * alpha_index_times_count);
  auto forward_pos = forward_offsets;
  auto reverse_pos = reverse_offsets;
  for (int k = 0; k < alpha_index_times_count; k++) {
    forward_rules[forward_pos[alpha_index_times[k][3]]++] = k;
    // Each operand occurrence contributes a derivative.
    reverse_rules[reverse_pos[alpha_index_times[k][0]]++] = 2 * k;
    reverse_rules[reverse_pos[alpha_index_times[k][1]]++] = 2 * k + 1;
  }

  std::vector<int> pending(n), depth(n, 0), ready;
  ready.reserve(n);
  for (int i = 0; i < n; i++) {
    pending[i] = 2 * (forward_offsets[i + 1] - forward_offsets[i]);
    if (pending[i] == 0) ready.push_back(i);
  }
  num_waves = 1;
  for (size_t q = 0; q < ready.size(); q++) {
    const int i = ready[q];
    for (int e = reverse_offsets[i]; e < reverse_offsets[i + 1]; e++) {
      const int child = alpha_index_times[reverse_rules[e] / 2][3];
      depth[child] = MAX(depth[child], depth[i] + 1);
      num_waves = MAX(num_waves, depth[child] + 1);
      if (--pending[child] == 0) ready.push_back(child);
    }
  }
  if ((int) ready.size() != n) error->all(FLERR, "Cyclic MTP contraction graph.");

  std::vector<int> waves(num_waves + 1, 0), nodes(n);
  for (int i = 0; i < n; i++) waves[depth[i] + 1]++;
  for (int i = 0; i < num_waves; i++) waves[i + 1] += waves[i];
  auto wave_pos = waves;
  for (int i = 0; i < n; i++) nodes[wave_pos[depth[i]]++] = i;

  std::vector<int> long_waves(num_waves + 1, 0), long_nodes;
  for (int wave = 0; wave < num_waves; wave++) {
    for (int q = waves[wave]; q < waves[wave + 1]; q++) {
      const int node = nodes[q];
      if (reverse_offsets[node + 1] - reverse_offsets[node] > REVERSE_LONG_THRESHOLD)
        long_nodes.push_back(node);
    }
    long_waves[wave + 1] = long_nodes.size();
  }

  std::vector<int> reverse_split(n);
  for (int i = 0; i < n; i++) {
    const auto split =
        std::stable_partition(reverse_rules.begin() + reverse_offsets[i],
                              reverse_rules.begin() + reverse_offsets[i + 1], [&](const int use) {
                                const int child = alpha_index_times[use / 2][3];
                                return reverse_offsets[child] != reverse_offsets[child + 1];
                              });
    reverse_split[i] = split - reverse_rules.begin();
  }

  auto copy_indices = [](auto &dst, const std::vector<int> &src, const char *label) {
    MemKK::realloc_kokkos(dst, label, src.size());
    auto h_dst = Kokkos::create_mirror_view(dst);
    for (size_t i = 0; i < src.size(); i++) h_dst(i) = src[i];
    Kokkos::deep_copy(dst, h_dst);
  };
  copy_indices(h_waves, waves, "mtp/kk:h_waves");
  copy_indices(d_wave_nodes, nodes, "mtp/kk:wave_nodes");
  copy_indices(d_forward_offsets, forward_offsets, "mtp/kk:forward_offsets");
  copy_indices(d_forward_rules, forward_rules, "mtp/kk:forward_rules");
  copy_indices(d_reverse_offsets, reverse_offsets, "mtp/kk:reverse_offsets");
  copy_indices(d_reverse_split, reverse_split, "mtp/kk:reverse_split");
  copy_indices(h_long_waves, long_waves, "mtp/kk:long_waves");
  copy_indices(d_long_nodes, long_nodes, "mtp/kk:long_nodes");
  MemKK::realloc_kokkos(d_reverse_terms, "mtp/kk:reverse_terms", reverse_rules.size());
  auto h_reverse_terms = Kokkos::create_mirror_view(d_reverse_terms);
  for (size_t e = 0; e < reverse_rules.size(); e++) {
    const int use = reverse_rules[e];
    const int k = use / 2;
    h_reverse_terms(e, 0) = alpha_index_times[k][3];
    h_reverse_terms(e, 1) = alpha_index_times[k][1 - use % 2];
    h_reverse_terms(e, 2) = alpha_index_times[k][2];
  }
  Kokkos::deep_copy(d_reverse_terms, h_reverse_terms);
}

// Finds the maximum number of valid MTP neighbours in all neigbhourhoods.
template <class DeviceType> struct FindMaxValidNeighsMTP {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  const KK_FLOAT max_cutoff_sq;
  Kokkos::View<int *, DeviceType> d_num_valid_neighs;
  Kokkos::View<int **, DeviceType> d_valid_neighs;
  const int capacity;

  FindMaxValidNeighsMTP(typename AT::t_int_1d_randomread d_ilist,
                        typename AT::t_int_1d_randomread d_numneigh,
                        typename AT::t_neighbors_2d d_neighbors,
                        typename AT::t_kkfloat_1d_3_lr_randomread x, KK_FLOAT max_cutoff_sq,
                        Kokkos::View<int *, DeviceType> d_num_valid_neighs,
                        Kokkos::View<int **, DeviceType> d_valid_neighs, int capacity) :
      d_ilist(d_ilist), d_numneigh(d_numneigh), d_neighbors(d_neighbors), x(x),
      max_cutoff_sq(max_cutoff_sq), d_num_valid_neighs(d_num_valid_neighs),
      d_valid_neighs(d_valid_neighs), capacity(capacity)
  {
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<DeviceType>::member_type &team,
                  int &max_valid_neighs) const
  {
    const int ii = team.league_rank();
    const int i = d_ilist[ii];
    const int num_neighs = d_numneigh(i);

    if (num_neighs == 0) {
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        d_num_valid_neighs(ii) = 0;
      });
      if (max_valid_neighs < 0) max_valid_neighs = 0;
      return;
    }

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
                              if (is_valid && pos < capacity) { d_valid_neighs(pos, ii) = j; }

                              // The last iteration's final prefix is the total number of valid neighbors.
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
  ev_init(eflag, vflag, 0);

  typename DeviceType::execution_space graph_space;

  // Reuse per-atom storage, but clear each requested tally on every compute.
  // Clearing the full views also resets ghost entries used with Newton on.
  if (eflag_atom) {
    if ((int) d_eatom.extent(0) < maxeatom) {
      memoryKK->destroy_kokkos(k_eatom, eatom);
      memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    }
    d_eatom = k_eatom.view<DeviceType>();
    Kokkos::deep_copy(graph_space, d_eatom, KK_ACC_FLOAT(0));
  }
  if (vflag_atom) {
    if ((int) d_vatom.extent(0) < maxvatom) {
      memoryKK->destroy_kokkos(k_vatom, vatom);
      memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
    }
    d_vatom = k_vatom.view<DeviceType>();
    Kokkos::deep_copy(graph_space, d_vatom, KK_ACC_FLOAT(0));
  }

  copymode = 1;

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
    if (vflag_atom)
      dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    if (vflag_atom)
      ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
    // clang-format on
  }

  if ((int) d_num_valid_neighs.extent(0) < inum) {
    Kokkos::realloc(Kokkos::WithoutInitializing, d_num_valid_neighs, inum);
  }
  if ((int) d_valid_neighs.extent(1) < inum) {
    Kokkos::realloc(Kokkos::WithoutInitializing, d_valid_neighs,
                    MAX((int) d_valid_neighs.extent(0), 1), inum);
  }

  if (!ts_basic) {
    using BasicPolicy = Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaBasic>;
    const int basic_limit =
        BasicPolicy(inum, Kokkos::AUTO).team_size_max(*this, Kokkos::ParallelForTag());
    ts_basic = MIN(basic_limit, MAX_TEAM_SIZE_BASIC);
    if (ts_basic < 1) error->all(FLERR, "Insufficient device resources for MTP basic alpha.");
  }
  const int team_size_default = ts_basic;

  // One pass counts, compacts and reduces the max.
  // Reruns if capactiy from last step was too small
  for (int attempt = 0; attempt < 2; attempt++) {
    max_valid_neighs = 0;
    Kokkos::TeamPolicy<DeviceType> policy_valid_neighs(inum, team_size_default);
    Kokkos::parallel_reduce("PairMTPKokkos::find_valid_neighs", policy_valid_neighs,
                            FindMaxValidNeighsMTP<DeviceType>(
                                d_ilist, d_numneigh, d_neighbors, x, max_cutoff_sq,
                                d_num_valid_neighs, d_valid_neighs, (int) d_valid_neighs.extent(0)),
                            Kokkos::Max<int>(max_valid_neighs));
    if (max_valid_neighs <= (int) d_valid_neighs.extent(0)) break;
    const int grown = max_valid_neighs + max_valid_neighs / 8 + 1;
    Kokkos::realloc(
        Kokkos::WithoutInitializing, d_valid_neighs,
        ((grown + NEIGH_CAPACITY_ALIGN - 1) / NEIGH_CAPACITY_ALIGN) * NEIGH_CAPACITY_ALIGN, inum);
  }

  // Handling batching
  chunk_size = MIN(input_chunk_size,
                   inum);    // chunksize is the maximum atoms per pass as defined by the user
  chunk_offset = 0;

  // Resize the arrays to the chunksize if needed. Do not initialize.
  const int capacity_tiles = (chunk_size + ATOM_TILE_SIZE - 1) / ATOM_TILE_SIZE;
  if ((int) d_moment_tensor_vals.extent(0) < capacity_tiles) {
    Kokkos::realloc(Kokkos::WithoutInitializing, d_moment_tensor_vals, capacity_tiles,
                    alpha_moment_count);
    Kokkos::realloc(Kokkos::WithoutInitializing, d_nbh_energy_ders_wrt_moments, capacity_tiles,
                    alpha_moment_count);
  }

  // Preserve capacity when either the neighbor count or chunk size grows.
  if ((int) d_radial_vals.extent(0) < max_valid_neighs ||
      (int) d_radial_vals.extent(2) < chunk_size) {
    const int cache_neighs = MAX(max_valid_neighs, (int) d_radial_vals.extent(0));
    const int cache_atoms = MAX(chunk_size, (int) d_radial_vals.extent(2));
    Kokkos::realloc(Kokkos::WithoutInitializing, d_radial_vals, cache_neighs, radial_func_count,
                    cache_atoms);
    Kokkos::realloc(Kokkos::WithoutInitializing, d_radial_ders, cache_neighs, radial_func_count,
                    cache_atoms);
    Kokkos::realloc(Kokkos::WithoutInitializing, d_inv_dist, cache_neighs, cache_atoms);
  }

  EV_FLOAT ev;

  using TimesPolicy = Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaTimes>;
  using NbhDersPolicy = Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeNbhDers>;
  using LongPolicy = Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeNbhDersLong>;
  const int long_scratch = shared_kk_float_1d::shmem_size(MAX_TEAM_SIZE_GRAPH * ATOM_TILE_SIZE);
  if (d_long_nodes.extent(0) &&
      (size_t) long_scratch >
          LongPolicy(graph_space, 1, Kokkos::AUTO, ATOM_TILE_SIZE).scratch_size_max(0))
    error->all(FLERR, "Insufficient scratch memory for MTP long reverse accumulation.");

  // Occupancy Queries.
  if (!ts_times) {
    const int limit_times = TimesPolicy(graph_space, 1, Kokkos::AUTO, ATOM_TILE_SIZE)
                                .team_size_max(*this, Kokkos::ParallelForTag());
    const int limit_nbh = NbhDersPolicy(graph_space, 1, Kokkos::AUTO, ATOM_TILE_SIZE)
                              .team_size_max(*this, Kokkos::ParallelForTag());
    ts_times = MIN(limit_times, MAX_TEAM_SIZE_GRAPH);
    ts_nbh = MIN(limit_nbh, MAX_TEAM_SIZE_GRAPH);
    ts_nbh_long = 1;
    if (d_long_nodes.extent(0)) {
      LongPolicy long_limits(graph_space, 1, Kokkos::AUTO, ATOM_TILE_SIZE);
      long_limits.set_scratch_size(0, Kokkos::PerTeam(long_scratch));
      const int limit_long = long_limits.team_size_max(*this, Kokkos::ParallelForTag());
      ts_nbh_long = MIN(limit_long, MAX_TEAM_SIZE_GRAPH);
    }
    if (ts_times < 1 || ts_nbh < 1 || ts_nbh_long < 1)
      error->all(FLERR, "Insufficient device resources for MTP graph kernels.");
  }

  const int radial_scratch_count = radial_func_count;
  const int coords_scratch_count = 3 * max_alpha_index_basic;
  const int basis_scratch_count = 2 * radial_basis_size;

  const int basic_scratch_size = scratch_size_helper<KK_FLOAT>(
      max_valid_neighs * (radial_scratch_count + coords_scratch_count) +
      Kokkos::min(team_size_default, max_valid_neighs) * basis_scratch_count);

  int force_team_size = team_size_default;
  if (max_valid_neighs < MAX_TEAM_SIZE_BASIC)
    force_team_size = MIN(force_team_size,
                          ((max_valid_neighs + NEIGH_CAPACITY_ALIGN - 1) / NEIGH_CAPACITY_ALIGN) *
                              NEIGH_CAPACITY_ALIGN);
  force_team_size = MAX(force_team_size, 1);

  const int force_scratch_size =
      scratch_size_helper<KK_FLOAT>(2 * radial_func_count + 3 * max_alpha_index_basic);
  const int force_team_scratch_size = shared_kk_float_1d::shmem_size(alpha_index_basic_count);

  // ========== Begin Main Computation ==========
  while (chunk_offset < inum) {    // batching to prevent OOM on device
    EV_FLOAT ev_tmp;
    if (chunk_size > inum - chunk_offset) chunk_size = inum - chunk_offset;
    const int atom_tiles = (chunk_size + ATOM_TILE_SIZE - 1) / ATOM_TILE_SIZE;
    const int graph_partitions =
        MIN(GRAPH_PARTITION_MAX, MAX(1, (GRAPH_PARTITION_TARGET + atom_tiles - 1) / atom_tiles));

    // ========== Calculate the basic alphas (Per outer-atom parallelizaton) ==========
    {
      Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaBasic> policy_basic_alpha(
          graph_space, chunk_size, team_size_default);
      if ((size_t) basic_scratch_size > policy_basic_alpha.scratch_size_max(0))
        error->all(FLERR, "Insufficient scratch memory for MTP basic alpha computation.");
      policy_basic_alpha =
          policy_basic_alpha.set_scratch_size(0, Kokkos::PerTeam(basic_scratch_size));
      Kokkos::parallel_for("ComputeAlphaBasic", policy_basic_alpha, *this);
    }

    // ========== Calculate the composite moment values  ==========
    {
      const int team_size = ts_times;
      for (int wave = 0; wave < num_waves; wave++) {
        if (!eflag_either && wave == num_waves - 1) break;
        wave_begin = MAX(alpha_index_basic_count, h_waves(wave));
        wave_end = h_waves(wave + 1);
        if (wave_begin == wave_end) continue;
        node_partitions =
            MIN(graph_partitions, (wave_end - wave_begin + team_size - 1) / team_size);
        TimesPolicy policy(graph_space, atom_tiles * node_partitions, team_size, ATOM_TILE_SIZE);
        Kokkos::parallel_for("ComputeAlphaTimes", policy, *this);
      }
    }

    // ========== Calc the nbh ders wrt moments ==========
    {
      const int team_size = ts_nbh;
      const int long_team_size = ts_nbh_long;
      for (int wave = num_waves - 2; wave >= 0; wave--) {
        wave_begin = h_waves(wave);
        wave_end = h_waves(wave + 1);
        if (wave_begin == wave_end) continue;
        node_partitions =
            MIN(graph_partitions, (wave_end - wave_begin + team_size - 1) / team_size);
        NbhDersPolicy policy(graph_space, atom_tiles * node_partitions, team_size, ATOM_TILE_SIZE);
        Kokkos::parallel_for("ComputeNbhDers", policy, *this);
        wave_begin = h_long_waves(wave);
        wave_end = h_long_waves(wave + 1);
        if (wave_begin < wave_end) {
          LongPolicy long_policy(graph_space, atom_tiles * (wave_end - wave_begin), long_team_size,
                                 ATOM_TILE_SIZE);
          long_policy.set_scratch_size(0, Kokkos::PerTeam(long_scratch));
          Kokkos::parallel_for("ComputeNbhDersLong", long_policy, *this);
        }
      }
    }

    // ========== Compute force (and dot product with alphas to get energy if needed) ==========
    {
      if (neighflag == HALF) {
        if (evflag)
          ev_tmp = compute_force<HALF, 1>(graph_space, force_team_size, force_team_scratch_size,
                                          force_scratch_size);
        else
          ev_tmp = compute_force<HALF, 0>(graph_space, force_team_size, force_team_scratch_size,
                                          force_scratch_size);
      } else if (neighflag == HALFTHREAD) {
        if (evflag)
          ev_tmp = compute_force<HALFTHREAD, 1>(graph_space, force_team_size,
                                                force_team_scratch_size, force_scratch_size);
        else
          ev_tmp = compute_force<HALFTHREAD, 0>(graph_space, force_team_size,
                                                force_team_scratch_size, force_scratch_size);
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

template <class DeviceType>
template <int NEIGHFLAG, int EVFLAG>
EV_FLOAT PairMTPKokkos<DeviceType>::compute_force(const typename DeviceType::execution_space &space,
                                                  int team_size, int team_scratch,
                                                  int thread_scratch)
{
  using ForcePolicy = Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeForce<NEIGHFLAG, EVFLAG>>;
  int &cached = ts_force[NEIGHFLAG == HALFTHREAD][EVFLAG];

  if (!cached) {
    const int probe_scratch = team_scratch + team_size * thread_scratch;
    ForcePolicy probe(space, MAX(chunk_size, 1), Kokkos::AUTO);
    if ((size_t) probe_scratch > probe.scratch_size_max(0))
      error->all(FLERR, "Insufficient scratch memory for MTP force recomputation.");
    probe.set_scratch_size(0, Kokkos::PerTeam(probe_scratch));
    if constexpr (EVFLAG)
      cached = probe.team_size_max(*this, Kokkos::ParallelReduceTag());
    else
      cached = probe.team_size_max(*this, Kokkos::ParallelForTag());
    if (cached < 1) error->all(FLERR, "Insufficient device resources for MTP force recomputation.");
  }
  if (cached < team_size) team_size = cached;

  int league = chunk_size;
  if (league > FORCE_MAX_BLOCKS) league = FORCE_MAX_BLOCKS;

  ForcePolicy policy =
      ForcePolicy(space, league, team_size)
          .set_scratch_size(0, Kokkos::PerTeam(team_scratch + team_size * thread_scratch));

  EV_FLOAT ev = {};
  if constexpr (EVFLAG)
    Kokkos::parallel_reduce("ComputeForceEV", policy, *this, ev);
  else
    Kokkos::parallel_for("ComputeForce", policy, *this);
  return ev;
}

// ========== Kernels ==========

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
  const int itype = d_map(type(i));
  const int jnum = d_num_valid_neighs(ii + chunk_offset);

  const int array_size = Kokkos::min(team.team_size(), jnum);
  const int power_stride = max_alpha_index_basic;
  shared_kk_float_2d s_radial_vals(team.team_scratch(0), radial_func_count, jnum);
  shared_kk_float_3d s_coord_powers(team.team_scratch(0), power_stride, jnum);
  shared_kk_float_2d s_radial_basis_vals(team.team_scratch(0), array_size, radial_basis_size);
  shared_kk_float_2d s_radial_basis_ders(team.team_scratch(0), array_size, radial_basis_size);

  // First we calc every neighbour's radial functions and powers.
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, jnum), [&](const int jj) {
    const int j = d_valid_neighs(jj, ii + chunk_offset);
    const int jtype = d_map(type(j));
    const KK_FLOAT r[3] = {x(j, 0) - xi[0], x(j, 1) - xi[1], x(j, 2) - xi[2]};
    const KK_FLOAT rsq = Kokkos::fma(r[0], r[0], Kokkos::fma(r[1], r[1], r[2] * r[2]));
    const KK_FLOAT dist = Kokkos::sqrt(rsq);
    const KK_FLOAT inv_dist = 1.0 / dist;
    const KK_FLOAT u[3] = {r[0] * inv_dist, r[1] * inv_dist, r[2] * inv_dist};
    d_inv_dist(jj, ii) = inv_dist;

    s_coord_powers(0, jj, 0) = s_coord_powers(0, jj, 1) = s_coord_powers(0, jj, 2) =
        1;    // Set the constants

    // Powers of the unit vector already carry the rank normalization
    for (int k = 1; k < max_alpha_index_basic; k++) {
      for (int a = 0; a < 3; a++) s_coord_powers(k, jj, a) = s_coord_powers(k - 1, jj, a) * u[a];
    }

    // Calculate the radial basis and store in shared memory.
    const KK_FLOAT mult = radial_mult;
    const KK_FLOAT ksi = Kokkos::fma(2.0, dist, -cutoff_sum) * inv_cutoff_range;

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

    // Precompute values
    int pair_offset = itype * species_count + jtype;
    for (int mu = 0; mu < radial_func_count; mu++) {
      KK_FLOAT val = 0;
      KK_FLOAT der = 0;
      int offset = (pair_offset * radial_basis_size * radial_func_count) + mu * radial_basis_size;

      for (int ri = 0; ri < radial_basis_size; ri++) {
        val = Kokkos::fma(d_radial_basis_coeffs(offset + ri), s_radial_basis_vals(thread, ri), val);
        der = Kokkos::fma(d_radial_basis_coeffs(offset + ri), s_radial_basis_ders(thread, ri), der);
      }

      s_radial_vals(mu, jj) = val;
      d_radial_vals(jj, mu, ii) = val;
      d_radial_ders(jj, mu, ii) = der;
    }
  });

  team.team_barrier();

  // Each thread sums over an alpha basic
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, alpha_index_basic_count), [&](const int k) {
    const int mu = d_alpha_index_basic(k, 0);
    const int a0 = d_alpha_index_basic(k, 1);
    const int a1 = d_alpha_index_basic(k, 2);
    const int a2 = d_alpha_index_basic(k, 3);

    KK_FLOAT moment_val = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const KK_FLOAT pow =
          s_coord_powers(a0, jj, 0) * s_coord_powers(a1, jj, 1) * s_coord_powers(a2, jj, 2);
      moment_val = Kokkos::fma(s_radial_vals(mu, jj), pow, moment_val);
    }

    d_moment_tensor_vals(ii / ATOM_TILE_SIZE, k, ii % ATOM_TILE_SIZE) = moment_val;
    if (d_reverse_offsets(k) == d_reverse_offsets(k + 1))
      d_nbh_energy_ders_wrt_moments(ii / ATOM_TILE_SIZE, k, ii % ATOM_TILE_SIZE) =
          d_moment_coeffs(k);
  });
}

// Calculates the non-elementary alpha from the basic alphas
template <class DeviceType>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(
    TagPairMTPComputeAlphaTimes,
    const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeAlphaTimes>::member_type &team)
    const
{
  const int tile = team.league_rank() / node_partitions;
  const int part = team.league_rank() % node_partitions;
  const int count = (wave_end - wave_begin + node_partitions - 1 - part) / node_partitions;

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, count), [&](const int q) {
    const int node = d_wave_nodes(wave_begin + part + q * node_partitions);
    if (!eflag_either && d_reverse_offsets(node) == d_reverse_offsets(node + 1)) return;
    const int begin = d_forward_offsets(node);
    const int end = d_forward_offsets(node + 1);
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, ATOM_TILE_SIZE), [&](const int lane) {
      if (tile * ATOM_TILE_SIZE + lane >= chunk_size) return;
      KK_FLOAT sum = 0;
      for (int e = begin; e < end; e++) {
        const int k = d_forward_rules(e);
        const int a0 = d_alpha_index_times(k, 0);
        const int a1 = d_alpha_index_times(k, 1);
        const int mult = d_alpha_index_times(k, 2);
        sum += mult * d_moment_tensor_vals(tile, a0, lane) * d_moment_tensor_vals(tile, a1, lane);
      }
      d_moment_tensor_vals(tile, node, lane) = sum;
    });
  });
}

// Calculates the nbh ders (backwards pass)
template <class DeviceType>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(
    TagPairMTPComputeNbhDers,
    const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeNbhDers>::member_type &team)
    const
{
  const int tile = team.league_rank() / node_partitions;
  const int part = team.league_rank() % node_partitions;
  const int count = (wave_end - wave_begin + node_partitions - 1 - part) / node_partitions;

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, count), [&](const int q) {
    const int node = d_wave_nodes(wave_begin + part + q * node_partitions);
    const int begin = d_reverse_offsets(node);
    const int end = d_reverse_offsets(node + 1);
    if (begin == end || end - begin > REVERSE_LONG_THRESHOLD) return;
    const int split = d_reverse_split(node);
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, ATOM_TILE_SIZE), [&](const int lane) {
      if (tile * ATOM_TILE_SIZE + lane >= chunk_size) return;
      KK_FLOAT sum = d_moment_coeffs(node);
      for (int e = begin; e < split; e++) {
        const int child = d_reverse_terms(e, 0);
        const int partner = d_reverse_terms(e, 1);
        const int mult = d_reverse_terms(e, 2);
        sum += d_nbh_energy_ders_wrt_moments(tile, child, lane) * mult *
            d_moment_tensor_vals(tile, partner, lane);
      }
      for (int e = split; e < end; e++) {
        const int child = d_reverse_terms(e, 0);
        const int partner = d_reverse_terms(e, 1);
        const int mult = d_reverse_terms(e, 2);
        sum += d_moment_coeffs(child) * mult * d_moment_tensor_vals(tile, partner, lane);
      }
      d_nbh_energy_ders_wrt_moments(tile, node, lane) = sum;
    });
  });
}

// Long lists share their terms across the team.
template <class DeviceType>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(
    TagPairMTPComputeNbhDersLong,
    const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeNbhDersLong>::member_type &team)
    const
{
  const int count = wave_end - wave_begin;
  const int tile = team.league_rank() / count;
  const int node = d_long_nodes(wave_begin + team.league_rank() % count);
  const int thread = team.team_rank();
  const int begin = d_reverse_offsets(node);
  const int end = d_reverse_offsets(node + 1);
  const int split = d_reverse_split(node);
  shared_kk_float_1d partial(team.team_scratch(0), team.team_size() * ATOM_TILE_SIZE);

  Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, ATOM_TILE_SIZE), [&](const int lane) {
    if (tile * ATOM_TILE_SIZE + lane >= chunk_size) return;
    KK_FLOAT sum = 0;
    for (int e = begin + thread; e < split; e += team.team_size()) {
      const int child = d_reverse_terms(e, 0);
      const int partner = d_reverse_terms(e, 1);
      const int mult = d_reverse_terms(e, 2);
      sum += d_nbh_energy_ders_wrt_moments(tile, child, lane) * mult *
          d_moment_tensor_vals(tile, partner, lane);
    }
    for (int e = split + thread; e < end; e += team.team_size()) {
      const int child = d_reverse_terms(e, 0);
      const int partner = d_reverse_terms(e, 1);
      const int mult = d_reverse_terms(e, 2);
      sum += d_moment_coeffs(child) * mult * d_moment_tensor_vals(tile, partner, lane);
    }
    partial(thread * ATOM_TILE_SIZE + lane) = sum;
  });
  team.team_barrier();

  if (thread == 0) {
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, ATOM_TILE_SIZE), [&](const int lane) {
      if (tile * ATOM_TILE_SIZE + lane >= chunk_size) return;
      KK_FLOAT sum = d_moment_coeffs(node);
      for (int t = 0; t < team.team_size(); t++) sum += partial(t * ATOM_TILE_SIZE + lane);
      d_nbh_energy_ders_wrt_moments(tile, node, lane) = sum;
    });
  }
}

// Computes forces from  radial functions and ders
template <class DeviceType>
template <int NEIGHFLAG>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(
    const TagPairMTPComputeForce<NEIGHFLAG, 0> &tag,
    const typename Kokkos::TeamPolicy<DeviceType, TagPairMTPComputeForce<NEIGHFLAG, 0>>::member_type
        &team) const
{
  EV_FLOAT ev = {};
  operator()(tag, team, ev);
}

template <class DeviceType>
template <int NEIGHFLAG, int EVFLAG>
KOKKOS_INLINE_FUNCTION void PairMTPKokkos<DeviceType>::operator()(
    const TagPairMTPComputeForce<NEIGHFLAG, EVFLAG> &,
    const typename Kokkos::TeamPolicy<DeviceType,
                                      TagPairMTPComputeForce<NEIGHFLAG, EVFLAG>>::member_type &team,
    EV_FLOAT &ev) const
{
  // The f array is duplicated for OpenMP, atomic for GPU, and neither for Serial
  auto v_f =
      ScatterViewHelper<NeedDup_v<NEIGHFLAG, DeviceType>, decltype(dup_f), decltype(ndup_f)>::get(
          dup_f, ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG, DeviceType>>();

  const int thread = team.team_rank();

  shared_kk_float_1d s_basic_adj(team.team_scratch(0), alpha_index_basic_count);
  shared_kk_float_2d s_radial_vals(team.team_scratch(0), team.team_size(), radial_func_count);
  shared_kk_float_2d s_radial_ders(team.team_scratch(0), team.team_size(), radial_func_count);
  shared_kk_float_3d s_coord_powers(team.team_scratch(0), team.team_size(), max_alpha_index_basic);

  const bool need_energies = EVFLAG && eflag_either;
  const bool need_virial = EVFLAG && vflag_either;

  for (int ii = team.league_rank(); ii < chunk_size; ii += team.league_size()) {
    const int i = d_ilist[ii + chunk_offset];
    const int jnum = d_num_valid_neighs(ii + chunk_offset);
    const KK_FLOAT xi[3] = {x(i, 0), x(i, 1), x(i, 2)};

    // Reuse the central atom's adjoints across neighbours.
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, alpha_index_basic_count), [&](const int k) {
      s_basic_adj(k) = d_nbh_energy_ders_wrt_moments(ii / ATOM_TILE_SIZE, k, ii % ATOM_TILE_SIZE);
    });
    team.team_barrier();

    KK_ACC_FLOAT fix = 0;
    KK_ACC_FLOAT fiy = 0;
    KK_ACC_FLOAT fiz = 0;

    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, jnum),
        [&](const int jj, KK_ACC_FLOAT &fx, KK_ACC_FLOAT &fy, KK_ACC_FLOAT &fz) {
          const int j = d_valid_neighs(jj, ii + chunk_offset);
          const KK_FLOAT r[3] = {x(j, 0) - xi[0], x(j, 1) - xi[1], x(j, 2) - xi[2]};
          const KK_FLOAT inv_dist = d_inv_dist(jj, ii);
          const KK_FLOAT u[3] = {r[0] * inv_dist, r[1] * inv_dist, r[2] * inv_dist};

          s_coord_powers(thread, 0, 0) = s_coord_powers(thread, 0, 1) =
              s_coord_powers(thread, 0, 2) = 1;

          // Precompute the unit vector powers
          for (int k = 1; k < max_alpha_index_basic; k++) {
            for (int a = 0; a < 3; a++)
              s_coord_powers(thread, k, a) = s_coord_powers(thread, k - 1, a) * u[a];
          }

          // Load radial functions
          for (int mu = 0; mu < radial_func_count; mu++) {
            s_radial_vals(thread, mu) = d_radial_vals(jj, mu, ii);
            s_radial_ders(thread, mu) = d_radial_ders(jj, mu, ii);
          }

          KK_ACC_FLOAT temp_force[3] = {0, 0, 0};
          KK_ACC_FLOAT radial_force = 0;

          // Recompute each derivative
          for (int k = 0; k < alpha_index_basic_count; k++) {

            int mu = d_alpha_index_basic(k, 0);
            int a0 = d_alpha_index_basic(k, 1);
            int a1 = d_alpha_index_basic(k, 2);
            int a2 = d_alpha_index_basic(k, 3);

            KK_FLOAT val = s_radial_vals(thread, mu);
            KK_FLOAT der = s_radial_ders(thread, mu);

            // Normalize by the rank
            int norm_rank = a0 + a1 + a2;
            der = Kokkos::fma(-norm_rank * inv_dist, val, der);

            KK_FLOAT pow0 = s_coord_powers(thread, a0, 0);
            KK_FLOAT pow1 = s_coord_powers(thread, a1, 1);
            KK_FLOAT pow2 = s_coord_powers(thread, a2, 2);
            KK_FLOAT pow = pow0 * pow1 * pow2;

            // Get the component's derivatives too
            KK_FLOAT adj = s_basic_adj(k);
            radial_force += adj * der * pow;
            val *= adj * inv_dist;

            if (a0 != 0)
              temp_force[0] += val * a0 * (s_coord_powers(thread, a0 - 1, 0) * pow1 * pow2);
            if (a1 != 0)
              temp_force[1] += val * a1 * (pow0 * s_coord_powers(thread, a1 - 1, 1) * pow2);
            if (a2 != 0)
              temp_force[2] += val * a2 * (pow0 * pow1 * s_coord_powers(thread, a2 - 1, 2));
          }

          // Apply the radial direction once after summing the moments.
          for (int a = 0; a < 3; a++) temp_force[a] += radial_force * u[a];

          fx += temp_force[0];
          fy += temp_force[1];
          fz += temp_force[2];
          a_f(j, 0) -= temp_force[0];
          a_f(j, 1) -= temp_force[1];
          a_f(j, 2) -= temp_force[2];

          if (need_virial)
            v_tally_xyz<NEIGHFLAG>(ev, i, j, temp_force[0], temp_force[1], temp_force[2], -r[0],
                                   -r[1], -r[2]);
        },
        fix, fiy, fiz);

    // A single team member updates the central atom
    Kokkos::single(Kokkos::PerTeam(team), [&]() {
      a_f(i, 0) += fix;
      a_f(i, 1) += fiy;
      a_f(i, 2) += fiz;
    });

    if (need_energies) {
      team.team_barrier();
      const int itype = d_map(type(i));
      KK_FLOAT nbh_energy = 0;

      // Reduction to find the dot product of the linear coeffs and the moment tensor vals
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team, alpha_scalar_count),
          [&](const int k, KK_FLOAT &sum) {
            sum += d_linear_coeffs(k) *
                d_moment_tensor_vals(ii / ATOM_TILE_SIZE, d_alpha_moment_mapping(k),
                                     ii % ATOM_TILE_SIZE);
          },
          nbh_energy);

      // A single team member updates the global array
      Kokkos::single(Kokkos::PerTeam(team), [&]() {
        nbh_energy += d_species_coeffs[itype];    // Essentially the reference energy
        if (eflag_global) ev.evdwl += nbh_energy;
        if (eflag_atom) d_eatom[i] = nbh_energy;
      });
    }

    // s_basic_adj is overwritten by the next atom's prefetch.
    team.team_barrier();
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
    // The vatom array is duplicated for OpenMP, atomic for GPU, and neither for Serial
    auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG, DeviceType>, decltype(dup_vatom),
                                     decltype(ndup_vatom)>::get(dup_vatom, ndup_vatom);
    auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG, DeviceType>>();

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
template <typename scratch_type>
int PairMTPKokkos<DeviceType>::scratch_size_helper(int values_per_team)
{
  typedef Kokkos::View<scratch_type *, Kokkos::DefaultExecutionSpace::scratch_memory_space,
                       Kokkos::MemoryTraits<Kokkos::Unmanaged>>
      ScratchViewType;

  return ScratchViewType::shmem_size(values_per_team);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairMTPKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairMTPKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
