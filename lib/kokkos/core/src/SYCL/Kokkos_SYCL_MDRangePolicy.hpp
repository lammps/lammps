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
  properties.default_largest_tile_size = 16;
  properties.default_tile_size         = 2;
  properties.max_total_tile_size       = properties.max_threads;

  auto device = space.sycl_queue().get_device();
  auto max_work_item_sizes =
      device.get_info<sycl::info::device::max_work_item_sizes<3>>();
  properties.max_threads_dimensions[0] = max_work_item_sizes[0];
  properties.max_threads_dimensions[1] = max_work_item_sizes[1];
  properties.max_threads_dimensions[2] = max_work_item_sizes[2];
  return properties;
}

// Settings for TeamMDRangePolicy
template <typename Rank, TeamMDRangeThreadAndVector ThreadAndVector>
struct ThreadAndVectorNestLevel<Rank, Kokkos::SYCL, ThreadAndVector>
    : AcceleratorBasedNestLevel<Rank, ThreadAndVector> {};

}  // namespace Impl
}  // Namespace Kokkos
#endif
