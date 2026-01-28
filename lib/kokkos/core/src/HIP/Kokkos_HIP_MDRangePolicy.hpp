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

// Settings for MDRangePolicy
template <>
inline TileSizeProperties get_tile_size_properties<HIP>(const HIP& space) {
  TileSizeProperties properties;
  const auto& device_prop              = space.hip_device_prop();
  properties.max_threads               = device_prop.maxThreadsPerBlock;
  properties.default_largest_tile_size = 16;
  properties.default_tile_size         = 4;
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
