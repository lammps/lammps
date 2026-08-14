// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_MDRANGEPOLICY_HPP_
#define KOKKOS_HIP_MDRANGEPOLICY_HPP_

#include <KokkosExp_MDRangePolicy.hpp>

namespace Kokkos {

template <>
struct default_outer_direction<HIP> {
  using type                     = Iterate;
  static constexpr Iterate value = Iterate::Left;
};

template <>
struct default_inner_direction<HIP> {
  using type                     = Iterate;
  static constexpr Iterate value = Iterate::Left;
};

namespace Impl {

template <>
struct TileSizeRecommended<Kokkos::HIP> {
  template <typename Policy>
  static auto get(Policy const&) {
    constexpr auto InnerDirection = Policy::inner_direction;
    constexpr int Rank            = Policy::rank;

    using tile_type = typename Policy::tile_type;

    if constexpr (InnerDirection == Iterate::Left) {
      if constexpr (Rank == 2) {
        return tile_type{64, 4};
      } else if constexpr (Rank == 3) {
        return tile_type{32, 2, 4};
      } else if constexpr (Rank == 4) {
        return tile_type{16, 4, 2, 2};
      } else if constexpr (Rank == 5) {
        return tile_type{16, 4, 2, 2, 1};
      } else if constexpr (Rank == 6) {
        return tile_type{8, 4, 2, 2, 2, 1};
      }
      tile_type tile_sizes{};
      for (int i = 0; i < Rank; ++i) {
        tile_sizes[i] = 2;
      }
      tile_sizes[0] = 16;
      return tile_sizes;
    } else {
      if constexpr (Rank == 2) {
        return tile_type{4, 64};
      } else if constexpr (Rank == 3) {
        return tile_type{4, 2, 32};
      } else if constexpr (Rank == 4) {
        return tile_type{2, 2, 4, 16};
      } else if constexpr (Rank == 5) {
        return tile_type{1, 2, 2, 4, 16};
      } else if constexpr (Rank == 6) {
        return tile_type{1, 2, 2, 2, 4, 8};
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

// Settings for MDRangePolicy
template <>
inline TileSizeProperties get_tile_size_properties<HIP>(const HIP& space) {
  TileSizeProperties properties;
  const auto& device_prop              = space.hip_device_prop();
  properties.max_threads               = device_prop.maxThreadsPerBlock;
  properties.max_total_tile_size       = HIPTraits::MaxThreadsPerBlock;
  properties.max_threads_dimensions[0] = device_prop.maxThreadsDim[0];
  properties.max_threads_dimensions[1] = device_prop.maxThreadsDim[1];
  properties.max_threads_dimensions[2] = device_prop.maxThreadsDim[2];
  return properties;
}

// Settings for TeamMDRangePolicy
template <typename Rank, TeamMDRangeThreadAndVector ThreadAndVector>
struct ThreadAndVectorNestLevel<Rank, HIP, ThreadAndVector>
    : AcceleratorBasedNestLevel<Rank, ThreadAndVector> {};

}  // Namespace Impl
}  // Namespace Kokkos
#endif
