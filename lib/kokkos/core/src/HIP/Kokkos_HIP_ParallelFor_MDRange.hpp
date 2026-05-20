// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_PARALLEL_FOR_MDRANGE_HPP
#define KOKKOS_HIP_PARALLEL_FOR_MDRANGE_HPP

#include <Kokkos_Parallel.hpp>

#include <HIP/Kokkos_HIP_BlockSize_Deduction.hpp>
#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <KokkosExp_MDRangePolicy.hpp>
#include <impl/KokkosExp_IterateTileGPU.hpp>

namespace Kokkos {
namespace Impl {

// ParallelFor
template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>, HIP> {
 public:
  using Policy       = Kokkos::MDRangePolicy<Traits...>;
  using functor_type = FunctorType;

 private:
  using array_index_type = typename Policy::array_index_type;
  using index_type       = typename Policy::index_type;
  using LaunchBounds     = typename Policy::launch_bounds;
  using MaxGridSize      = Kokkos::Array<index_type, 3>;
  using array_type       = typename Policy::point_type;

  const FunctorType m_functor;
  const Policy m_policy;
  const MaxGridSize m_max_grid_size;

  array_type m_lower;
  array_type m_upper;
  array_type m_extent;  // tile_size * num_tiles

 public:
  ParallelFor() = delete;

  inline __device__ void operator()() const {
    Kokkos::Impl::DeviceIterate<Policy::rank, array_index_type, index_type,
                                FunctorType, Policy::inner_direction,
                                typename Policy::work_tag>(m_lower, m_upper,
                                                           m_extent, m_functor)
        .exec_range();
  }

  inline void execute() const {
    using ClosureType = ParallelFor<FunctorType, Policy, HIP>;
    if (m_policy.m_num_tiles == 0) return;

    const auto [grid, block] =
        Kokkos::Impl::compute_device_launch_params(m_policy, m_max_grid_size);

    hip_parallel_launch<ClosureType, LaunchBounds>(
        *this, grid, block, 0, m_policy.space().impl_internal_space_instance(),
        false);
  }  // end execute

  ParallelFor(FunctorType const& arg_functor, Policy const& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_max_grid_size({
            static_cast<index_type>(
                m_policy.space().hip_device_prop().maxGridSize[0]),
            static_cast<index_type>(
                m_policy.space().hip_device_prop().maxGridSize[1]),
            static_cast<index_type>(
                m_policy.space().hip_device_prop().maxGridSize[2]),
        }) {
    // Initialize begins and ends based on layout
    // Swap the fastest indexes to x dimension
    for (array_index_type i = 0; i < Policy::rank; ++i) {
      if constexpr (Policy::inner_direction == Iterate::Left) {
        m_lower[i]  = m_policy.m_lower[i];
        m_upper[i]  = m_policy.m_upper[i];
        m_extent[i] = m_policy.m_tile[i] * m_policy.m_tile_end[i];
      } else {
        m_lower[i]  = m_policy.m_lower[Policy::rank - 1 - i];
        m_upper[i]  = m_policy.m_upper[Policy::rank - 1 - i];
        m_extent[i] = m_policy.m_tile[Policy::rank - 1 - i] *
                      m_policy.m_tile_end[Policy::rank - 1 - i];
      }
    }
  }

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy&, const Functor&) {
    using closure_type =
        ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>, HIP>;
    unsigned block_size = hip_get_max_blocksize<closure_type, LaunchBounds>();
    if (block_size == 0)
      Kokkos::Impl::throw_runtime_exception(
          std::string("Kokkos::Impl::ParallelFor< HIP > could not find a valid "
                      "tile size."));
    return block_size;
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
