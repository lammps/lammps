#ifndef KOKKOS_HIP_MDRANGEPOLICY_HPP_
#define KOKKOS_HIP_MDRANGEPOLICY_HPP_

#include <KokkosExp_MDRangePolicy.hpp>

namespace Kokkos {

template <>
struct default_outer_direction<Kokkos::Experimental::HIP> {
  using type                     = Iterate;
  static constexpr Iterate value = Iterate::Left;
};

template <>
struct default_inner_direction<Kokkos::Experimental::HIP> {
  using type                     = Iterate;
  static constexpr Iterate value = Iterate::Left;
};

namespace Impl {

// Settings for MDRangePolicy
template <>
inline TileSizeProperties get_tile_size_properties<Kokkos::Experimental::HIP>(
    const Kokkos::Experimental::HIP& space) {
  TileSizeProperties properties;
  properties.max_threads =
      space.impl_internal_space_instance()->m_maxThreadsPerSM;
  properties.default_largest_tile_size = 16;
  properties.default_tile_size         = 4;
  properties.max_total_tile_size =
      Kokkos::Experimental::Impl::HIPTraits::MaxThreadsPerBlock;
  return properties;
}

}  // Namespace Impl
}  // Namespace Kokkos
#endif
