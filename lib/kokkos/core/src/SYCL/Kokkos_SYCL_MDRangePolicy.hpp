// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SYCL_MDRANGEPOLICY_HPP_
#define KOKKOS_SYCL_MDRANGEPOLICY_HPP_

#include <KokkosExp_MDRangePolicy.hpp>

namespace Kokkos {

template <>
struct default_outer_direction<Kokkos::SYCL> {
  using type                     = Iterate;
  static constexpr Iterate value = Iterate::Left;
};

template <>
struct default_inner_direction<Kokkos::SYCL> {
  using type                     = Iterate;
  static constexpr Iterate value = Iterate::Left;
};

namespace Impl {

// Settings for MDRangePolicy
template <>
inline TileSizeProperties get_tile_size_properties<Kokkos::SYCL>(
    const Kokkos::SYCL& space) {
  TileSizeProperties properties;
  properties.max_threads =
      space.impl_internal_space_instance()->m_maxWorkgroupSize;
  properties.max_total_tile_size = properties.max_threads;
  auto device                    = space.sycl_queue().get_device();
  auto max_work_item_sizes =
      device.get_info<sycl::info::device::max_work_item_sizes<3>>();
  properties.max_threads_dimensions[0] = max_work_item_sizes[0];
  properties.max_threads_dimensions[1] = max_work_item_sizes[1];
  properties.max_threads_dimensions[2] = max_work_item_sizes[2];
  return properties;
}

template <>
struct TileSizeRecommended<Kokkos::SYCL> {
  template <typename Policy>
  static auto get(Policy const&) {
    constexpr auto InnerDirection = Policy::inner_direction;
    constexpr int Rank            = Policy::rank;

    using tile_type = typename Policy::tile_type;

    if constexpr (InnerDirection == Iterate::Left) {
      if constexpr (Rank == 2) {
        return tile_type{32, 8};
      } else if constexpr (Rank == 3) {
        return tile_type{32, 2, 4};
      } else if constexpr (Rank == 4) {
        return tile_type{16, 2, 2, 2};
      } else if constexpr (Rank == 5) {
        return tile_type{16, 2, 2, 2, 2};
      } else if constexpr (Rank == 6) {
        return tile_type{16, 2, 2, 2, 2, 1};
      }
      tile_type tile_sizes{};
      for (int i = 0; i < Rank; ++i) {
        tile_sizes[i] = 2;
      }
      tile_sizes[0] = 16;
      return tile_sizes;
    } else {
      if constexpr (Rank == 2) {
        return tile_type{8, 32};
      } else if constexpr (Rank == 3) {
        return tile_type{4, 2, 32};
      } else if constexpr (Rank == 4) {
        return tile_type{2, 2, 2, 16};
      } else if constexpr (Rank == 5) {
        return tile_type{2, 2, 2, 2, 16};
      } else if constexpr (Rank == 6) {
        return tile_type{1, 2, 2, 2, 2, 16};
      }
      tile_type tile_sizes{};
      for (int i = 0; i < Rank; ++i) {
        tile_sizes[i] = 2;
      }
      tile_sizes[Rank - 1] = 16;
      return tile_sizes;
    }
  }
};

// Settings for TeamMDRangePolicy
template <typename Rank, TeamMDRangeThreadAndVector ThreadAndVector>
struct ThreadAndVectorNestLevel<Rank, Kokkos::SYCL, ThreadAndVector>
    : AcceleratorBasedNestLevel<Rank, ThreadAndVector> {};

}  // namespace Impl
}  // Namespace Kokkos
#endif
