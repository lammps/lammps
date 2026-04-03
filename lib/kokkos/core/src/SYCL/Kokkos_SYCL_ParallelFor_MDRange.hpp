// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SYCL_PARALLEL_FOR_MDRANGE_HPP_
#define KOKKOS_SYCL_PARALLEL_FOR_MDRANGE_HPP_

#include <limits>

#include <sycl/sycl.hpp>

#include <impl/KokkosExp_IterateTileGPU.hpp>

#ifdef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
#include <vector>
#endif

template <class FunctorType, class... Traits>
class Kokkos::Impl::ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                                Kokkos::SYCL> {
 public:
  using Policy = Kokkos::MDRangePolicy<Traits...>;

 private:
  using array_index_type = typename Policy::array_index_type;
  using index_type       = typename Policy::index_type;
  using WorkTag          = typename Policy::work_tag;
  using MaxGridSize      = Kokkos::Array<index_type, 3>;
  using array_type       = typename Policy::point_type;

  const FunctorType m_functor;
  const Policy m_policy;
  const MaxGridSize m_max_grid_size;

  array_type m_lower;
  array_type m_upper;
  array_type m_extent;  // tile_size * num_tiles

  template <typename FunctorWrapper>
  sycl::event sycl_direct_launch(const FunctorWrapper& functor_wrapper,
                                 const sycl::event& memcpy_event) const {
    // Convenience references
    const Kokkos::SYCL& space = m_policy.space();
    sycl::queue& q            = space.sycl_queue();

    if (m_policy.m_num_tiles == 0) return {};

    const auto lower_bound = m_lower;
    const auto upper_bound = m_upper;
    const auto extent      = m_extent;

    desul::ensure_sycl_lock_arrays_on_device(q);

    const auto range =
        Kokkos::Impl::compute_device_launch_params(m_policy, m_max_grid_size);

    auto cgh_lambda = [&, range](sycl::handler& cgh) {
      const sycl::range<3> global_range = range.get_global_range();
      const sycl::range<3> local_range  = range.get_local_range();
      const sycl::nd_range sycl_swapped_range{
          sycl::range<3>{global_range[2], global_range[1], global_range[0]},
          sycl::range<3>{local_range[2], local_range[1], local_range[0]}};

#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
      cgh.depends_on(memcpy_event);
#else
      (void)memcpy_event;
#endif
      cgh.parallel_for(sycl_swapped_range, [lower_bound, upper_bound, extent,
                                            functor_wrapper](
                                               sycl::nd_item<3> item) {
        // swap back for correct index calculations in DeviceIterateTile
        const index_type local_x    = item.get_local_id(2);
        const index_type local_y    = item.get_local_id(1);
        const index_type local_z    = item.get_local_id(0);
        const index_type n_local_x  = item.get_local_range(2);
        const index_type n_local_y  = item.get_local_range(1);
        const index_type n_local_z  = item.get_local_range(0);
        const index_type global_x   = item.get_group(2);
        const index_type global_y   = item.get_group(1);
        const index_type global_z   = item.get_group(0);
        const index_type n_global_x = item.get_group_range(2);
        const index_type n_global_y = item.get_group_range(1);
        const index_type n_global_z = item.get_group_range(0);

        Kokkos::Impl::DeviceIterate<Policy::rank, array_index_type, index_type,
                                    FunctorType, Policy::inner_direction,
                                    typename Policy::work_tag>(
            lower_bound, upper_bound, extent, functor_wrapper.get_functor(),
            {n_global_x, n_global_y, n_global_z},
            {n_local_x, n_local_y, n_local_z}, {global_x, global_y, global_z},
            {local_x, local_y, local_z})
            .exec_range();
      });
    };

#ifdef KOKKOS_IMPL_SYCL_GRAPH_SUPPORT
    if constexpr (Policy::is_graph_kernel::value) {
      sycl_attach_kernel_to_node(*this, cgh_lambda);
      return {};
    } else
#endif
    {
      auto parallel_for_event = q.submit(cgh_lambda);

#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
      q.ext_oneapi_submit_barrier(std::vector<sycl::event>{parallel_for_event});
#endif
      return parallel_for_event;
    }
  }

  static MaxGridSize get_max_grid_size(const Policy& policy) {
    // the SYCL specs do not allow to get the maximum grid size (maximum
    // ND-range size, maximum number of work groups)
    // TODO update this when the specs change
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
    // we use an Intel extension if possible
    auto max_grid_size =
        policy.space()
            .sycl_queue()
            .get_device()
            .template get_info<sycl::ext::oneapi::experimental::info::device::
                                   max_work_groups<3>>();

    // note that SYCL represents a (x, y, z) range with the the right-most term
    // as the one varying the fastest, so the order must be reversed for Kokkos
    // see:
    // https://registry.khronos.org/SYCL/specs/sycl-2020/html/sycl-2020.html#sec:multi-dim-linearization
    return {static_cast<index_type>(max_grid_size[2]),
            static_cast<index_type>(max_grid_size[1]),
            static_cast<index_type>(max_grid_size[0])};
#else
    // otherwise, we consider that the max is infinite
    return {std::numeric_limits<index_type>::max(),
            std::numeric_limits<index_type>::max(),
            std::numeric_limits<index_type>::max()};
#endif
  }

 public:
  using functor_type = FunctorType;

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy& policy, const Functor&) {
    return policy.space().impl_internal_space_instance()->m_maxWorkgroupSize;
  }

  void execute() const {
    auto space_instance = m_policy.space().impl_internal_space_instance();
    Kokkos::Impl::SYCLInternal::IndirectKernelMem& indirectKernelMem =
        space_instance->get_indirect_kernel_mem();

    auto functor_wrapper =
        Impl::make_sycl_function_wrapper(m_functor, indirectKernelMem);
    sycl::event event =
        sycl_direct_launch(functor_wrapper, functor_wrapper.get_copy_event());
    functor_wrapper.register_event(event);
  }

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_max_grid_size(get_max_grid_size(arg_policy)) {
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
};

#endif  // KOKKOS_SYCL_PARALLEL_FOR_MDRANGE_HPP_
