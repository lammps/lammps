/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Evan Weinberg (NVIDIA)
------------------------------------------------------------------------- */

#include "kokkos_type.h"

namespace LAMMPS_NS {

// This helper class implements optimized out-of-place transposes for Rank-2 views.
// In the case where both views have the same layout, it uses Kokkos' default `deep_copy`.
// In the case where the views have different layouts (LayoutLeft/LayoutRight), it implements
// a transpose through scratch memory staging
template <class DeviceType_, class t_view_dst_, class t_view_src_>
struct TransposeHelperKokkos {

  using DeviceType = DeviceType_;

  using t_view_dst = t_view_dst_;
  using t_view_src = t_view_src_;

  static_assert(std::is_same<typename t_view_dst::value_type, typename t_view_src::value_type>::value, "Value types do not match");
  static_assert(t_view_dst::Rank == 2, "Destination view rank != 2");
  static_assert(t_view_src::Rank == 2, "Source view rank != 2");

  using dst_layout = typename t_view_dst::traits::array_layout;
  using src_layout = typename t_view_src::traits::array_layout;

  typedef ArrayTypes<DeviceType> AT;

  using t_view_value = typename t_view_dst::value_type;

  // 32x32 tiles, will update so each thread does multiple loads
  static constexpr int vector_length = 32;
  static constexpr int bank_pad = 1;
  static constexpr int elem_size = sizeof(t_view_value);

  static constexpr int threads_per_team = 4;

  t_view_dst d_dst;
  t_view_src d_src;

  bool src_is_layout_right;

  // extents divided by vector length, rounded up
  int extent_tiles[2];

  // 1 if extent is divisible by vector length, 0 otherwise
  int extent_is_multiple_vector_length[2];

  // number of teams
  int n_teams;

  // amount of shared memory per thread
  int shared_mem_per_thread;

  TransposeHelperKokkos(t_view_dst d_dst_, t_view_src d_src_)
   : d_dst(d_dst_), d_src(d_src_) {

    assert(d_dst.extent(0) == d_src.extent(0) && d_dst.extent(1) == d_dst.extent(1));

    if (std::is_same<dst_layout, src_layout>::value) {
      Kokkos::deep_copy(d_dst, d_src);
    } else {

      src_is_layout_right = std::is_same<src_layout, Kokkos::LayoutRight>::value;

      extent_tiles[0] = (d_dst.extent(0) + vector_length - 1) / vector_length;
      extent_tiles[1] = (d_dst.extent(1) + vector_length - 1) / vector_length;

      extent_is_multiple_vector_length[0] = (extent_tiles[0] * vector_length == d_dst.extent(0)) ? 1 : 0;
      extent_is_multiple_vector_length[1] = (extent_tiles[1] * vector_length == d_dst.extent(1)) ? 1 : 0;

      n_teams = (extent_tiles[0] * extent_tiles[1] + threads_per_team - 1) / threads_per_team;

      shared_mem_per_thread = vector_length * (vector_length + bank_pad) * elem_size;

      Kokkos::TeamPolicy<DeviceType> transpose_policy(n_teams, threads_per_team, vector_length);
      transpose_policy = transpose_policy.set_scratch_size(0, Kokkos::PerThread(shared_mem_per_thread));

      Kokkos::parallel_for(transpose_policy, *this);
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team_member) const {

    t_view_value* buffer = (t_view_value*)(team_member.team_shmem().get_shmem(shared_mem_per_thread * threads_per_team, 0)) + (shared_mem_per_thread / elem_size) * team_member.team_rank();

    // extract flattened tile
    const int flattened_idx = team_member.team_rank() + team_member.league_rank() * threads_per_team;

    // get range x, range y tile
    int extent_tile_id[2];
    if (src_is_layout_right) {
      // keep extent 1 tiles close together b/c loading from layout right
      extent_tile_id[0] = flattened_idx / extent_tiles[1];
      extent_tile_id[1] = flattened_idx - extent_tile_id[0] * extent_tiles[1];
    } else {
      // keep extent 0 tiles close together b/c loading from layout left
      extent_tile_id[1] = flattened_idx / extent_tiles[0];
      extent_tile_id[0] = flattened_idx - extent_tile_id[1] * extent_tiles[0];
    }

    int elem[2];
    elem[0] = extent_tile_id[0] * vector_length;
    elem[1] = extent_tile_id[1] * vector_length;

    if (elem[0] >= d_dst.extent(0) ||
      elem[1] >= d_dst.extent(1)) return;

    // determine if a row/column is a full `vector_length` in size or not
    bool perfect_pad[2];
    perfect_pad[0] = (extent_is_multiple_vector_length[0] == 1 || extent_tile_id[0] + 1 < extent_tiles[0]);
    perfect_pad[1] = (extent_is_multiple_vector_length[1] == 1 || extent_tile_id[1] + 1 < extent_tiles[1]);

    // load phase
    if (src_is_layout_right) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team_member, vector_length),
        [&] (const int j) {

        if (elem[1] + j < d_src.extent(1)) {
          if (perfect_pad[0]) {
            for (int i = 0; i < vector_length; i++)
              buffer[i * (vector_length + bank_pad) + j] = d_src(elem[0] + i, elem[1] + j);
          } else {
            for (int i = 0; i < (d_src.extent(0) - elem[0]); i++)
              buffer[i * (vector_length + bank_pad) + j] = d_src(elem[0] + i, elem[1] + j);
          }
        }
      });

    } else {
      // src is layout left
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team_member, vector_length),
        [&] (const int i) {

        if (elem[0] + i < d_src.extent(0)) {
          if (perfect_pad[1]) {
            for (int j = 0; j < vector_length; j++)
              buffer[i * (vector_length + bank_pad) + j] = d_src(elem[0] + i, elem[1] + j);
          } else {
            for (int j = 0; j < (d_src.extent(1) - elem[1]); j++)
              buffer[i * (vector_length + bank_pad) + j] = d_src(elem[0] + i, elem[1] + j);
          }
        }
      });
    }

    // No need for an extra sync b/c, as confirmed by asking on the Kokkos Slack, there
    // is an implicit sync at the end of a ThreadVectorRange as per the Kokkos
    // programming model.

    // save phase
    if (src_is_layout_right) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team_member, vector_length),
        [&] (const int i) {

        if (elem[0] + i < d_dst.extent(0)) {
          if (perfect_pad[1]) {
            for (int j = 0; j < vector_length; j++)
              d_dst(elem[0] + i, elem[1] + j) = buffer[i * (vector_length + bank_pad) + j];
          } else {
            for (int j = 0; j < (d_dst.extent(1) - elem[1]); j++)
              d_dst(elem[0] + i, elem[1] + j) = buffer[i * (vector_length + bank_pad) + j];
          }
        }
      });
    } else {

      // src is layout left
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team_member, vector_length),
        [&] (const int j) {

        if (elem[1] + j < d_dst.extent(1)) {
          if (perfect_pad[0]) {
            for (int i = 0; i < vector_length; i++)
              d_dst(elem[0] + i, elem[1] + j) = buffer[i * (vector_length + bank_pad) + j];
          } else {
            for (int i = 0; i < (d_dst.extent(0) - elem[0]); i++)
              d_dst(elem[0] + i, elem[1] + j) = buffer[i * (vector_length + bank_pad) + j];
          }
        }
      });
    }

  }
};

}
