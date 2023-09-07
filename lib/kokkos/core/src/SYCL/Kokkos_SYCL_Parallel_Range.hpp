//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_SYCL_PARALLEL_RANGE_HPP_
#define KOKKOS_SYCL_PARALLEL_RANGE_HPP_

#include <impl/KokkosExp_IterateTileGPU.hpp>

#include <vector>

namespace Kokkos::Impl {
template <typename FunctorWrapper, typename Policy>
struct FunctorWrapperRangePolicyParallelFor {
  using WorkTag = typename Policy::work_tag;

  void operator()(sycl::item<1> item) const {
    const typename Policy::index_type id = item.get_linear_id() + m_begin;
    if constexpr (std::is_void_v<WorkTag>)
      m_functor_wrapper.get_functor()(id);
    else
      m_functor_wrapper.get_functor()(WorkTag(), id);
  }

  typename Policy::index_type m_begin;
  FunctorWrapper m_functor_wrapper;
};

// Same as above but for a user-provided workgroup size
template <typename FunctorWrapper, typename Policy>
struct FunctorWrapperRangePolicyParallelForCustom {
  using WorkTag = typename Policy::work_tag;

  void operator()(sycl::item<1> item) const {
    const typename Policy::index_type id = item.get_linear_id();
    if (id < m_work_size) {
      const auto shifted_id = id + m_begin;
      if constexpr (std::is_void_v<WorkTag>)
        m_functor_wrapper.get_functor()(shifted_id);
      else
        m_functor_wrapper.get_functor()(WorkTag(), shifted_id);
    }
  }

  typename Policy::index_type m_begin;
  FunctorWrapper m_functor_wrapper;
  typename Policy::index_type m_work_size;
};
}  // namespace Kokkos::Impl

template <class FunctorType, class... Traits>
class Kokkos::Impl::ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>,
                                Kokkos::Experimental::SYCL> {
 public:
  using Policy = Kokkos::RangePolicy<Traits...>;

 private:
  using Member  = typename Policy::member_type;
  using WorkTag = typename Policy::work_tag;

  const FunctorType m_functor;
  const Policy m_policy;

  template <typename Functor>
  static sycl::event sycl_direct_launch(const Policy& policy,
                                        const Functor& functor,
                                        const sycl::event& memcpy_event) {
    // Convenience references
    const Kokkos::Experimental::SYCL& space = policy.space();
    sycl::queue& q                          = space.sycl_queue();

    auto parallel_for_event = q.submit([&](sycl::handler& cgh) {
      cgh.depends_on(memcpy_event);
      if (policy.chunk_size() <= 1) {
        FunctorWrapperRangePolicyParallelFor<Functor, Policy> f{policy.begin(),
                                                                functor};
        sycl::range<1> range(policy.end() - policy.begin());
        cgh.parallel_for<FunctorWrapperRangePolicyParallelFor<Functor, Policy>>(
            range, f);
      } else {
        // Use the chunk size as workgroup size. We need to make sure that the
        // range the kernel is launched with is a multiple of the workgroup
        // size. Hence, we need to restrict the execution of the functor in the
        // kernel to the actual range.
        const auto actual_range = policy.end() - policy.begin();
        const auto wgroup_size  = policy.chunk_size();
        const auto launch_range =
            (actual_range + wgroup_size - 1) / wgroup_size * wgroup_size;
        FunctorWrapperRangePolicyParallelForCustom<Functor, Policy> f{
            policy.begin(), functor, actual_range};
        sycl::nd_range<1> range(launch_range, wgroup_size);
        cgh.parallel_for<
            FunctorWrapperRangePolicyParallelForCustom<Functor, Policy>>(range,
                                                                         f);
      }
    });
    q.ext_oneapi_submit_barrier(std::vector<sycl::event>{parallel_for_event});

    return parallel_for_event;
  }

 public:
  using functor_type = FunctorType;

  void execute() const {
    if (m_policy.begin() == m_policy.end()) return;

    Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem&
        indirectKernelMem = m_policy.space()
                                .impl_internal_space_instance()
                                ->get_indirect_kernel_mem();

    auto functor_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_functor, indirectKernelMem);
    sycl::event event = sycl_direct_launch(m_policy, functor_wrapper,
                                           functor_wrapper.get_copy_event());
    functor_wrapper.register_event(event);
  }

  ParallelFor(const ParallelFor&) = delete;
  ParallelFor(ParallelFor&&)      = delete;
  ParallelFor& operator=(const ParallelFor&) = delete;
  ParallelFor& operator=(ParallelFor&&) = delete;
  ~ParallelFor()                        = default;

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

// ParallelFor
template <class FunctorType, class... Traits>
class Kokkos::Impl::ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                                Kokkos::Experimental::SYCL> {
 public:
  using Policy = Kokkos::MDRangePolicy<Traits...>;

 private:
  using array_index_type = typename Policy::array_index_type;
  using index_type       = typename Policy::index_type;
  using WorkTag          = typename Policy::work_tag;

  const FunctorType m_functor;
  // MDRangePolicy is not trivially copyable. Hence, replicate the data we
  // really need in DeviceIterateTile in a trivially copyable struct.
  const struct BarePolicy {
    using index_type = typename Policy::index_type;

    BarePolicy(const Policy& policy)
        : m_lower(policy.m_lower),
          m_upper(policy.m_upper),
          m_tile(policy.m_tile),
          m_tile_end(policy.m_tile_end),
          m_num_tiles(policy.m_num_tiles) {}

    const typename Policy::point_type m_lower;
    const typename Policy::point_type m_upper;
    const typename Policy::tile_type m_tile;
    const typename Policy::point_type m_tile_end;
    const typename Policy::index_type m_num_tiles;
    static constexpr Iterate inner_direction = Policy::inner_direction;
  } m_policy;
  const Kokkos::Experimental::SYCL& m_space;

  sycl::nd_range<3> compute_ranges() const {
    const auto& m_tile     = m_policy.m_tile;
    const auto& m_tile_end = m_policy.m_tile_end;

    if constexpr (Policy::rank == 2) {
      sycl::range<3> local_sizes(m_tile[0], m_tile[1], 1);
      sycl::range<3> global_sizes(m_tile_end[0] * m_tile[0],
                                  m_tile_end[1] * m_tile[1], 1);
      return {global_sizes, local_sizes};
    }
    if constexpr (Policy::rank == 3) {
      sycl::range<3> local_sizes(m_tile[0], m_tile[1], m_tile[2]);
      sycl::range<3> global_sizes(m_tile_end[0] * m_tile[0],
                                  m_tile_end[1] * m_tile[1],
                                  m_tile_end[2] * m_tile[2]);
      return {global_sizes, local_sizes};
    }
    if constexpr (Policy::rank == 4) {
      // id0,id1 encoded within first index; id2 to second index; id3 to third
      // index
      sycl::range<3> local_sizes(m_tile[0] * m_tile[1], m_tile[2], m_tile[3]);
      sycl::range<3> global_sizes(
          m_tile_end[0] * m_tile[0] * m_tile_end[1] * m_tile[1],
          m_tile_end[2] * m_tile[2], m_tile_end[3] * m_tile[3]);
      return {global_sizes, local_sizes};
    }
    if constexpr (Policy::rank == 5) {
      // id0,id1 encoded within first index; id2,id3 to second index; id4 to
      // third index
      sycl::range<3> local_sizes(m_tile[0] * m_tile[1], m_tile[2] * m_tile[3],
                                 m_tile[4]);
      sycl::range<3> global_sizes(
          m_tile_end[0] * m_tile[0] * m_tile_end[1] * m_tile[1],
          m_tile_end[2] * m_tile[2] * m_tile_end[3] * m_tile[3],
          m_tile_end[4] * m_tile[4]);
      return {global_sizes, local_sizes};
    }
    if constexpr (Policy::rank == 6) {
      // id0,id1 encoded within first index; id2,id3 to second index; id4,id5 to
      // third index
      sycl::range<3> local_sizes(m_tile[0] * m_tile[1], m_tile[2] * m_tile[3],
                                 m_tile[4] * m_tile[5]);
      sycl::range<3> global_sizes(
          m_tile_end[0] * m_tile[0] * m_tile_end[1] * m_tile[1],
          m_tile_end[2] * m_tile[2] * m_tile_end[3] * m_tile[3],
          m_tile_end[4] * m_tile[4] * m_tile_end[5] * m_tile[5]);
      return {global_sizes, local_sizes};
    }
    static_assert(Policy::rank > 1 && Policy::rank < 7,
                  "Kokkos::MDRange Error: Exceeded rank bounds with SYCL\n");
  }

  template <typename FunctorWrapper>
  sycl::event sycl_direct_launch(const FunctorWrapper& functor_wrapper,
                                 const sycl::event& memcpy_event) const {
    // Convenience references
    sycl::queue& q = m_space.sycl_queue();

    if (m_policy.m_num_tiles == 0) return {};

    const BarePolicy bare_policy(m_policy);

    auto parallel_for_event = q.submit([&](sycl::handler& cgh) {
      const auto range                  = compute_ranges();
      const sycl::range<3> global_range = range.get_global_range();
      const sycl::range<3> local_range  = range.get_local_range();
      const sycl::nd_range sycl_swapped_range{
          sycl::range<3>{global_range[2], global_range[1], global_range[0]},
          sycl::range<3>{local_range[2], local_range[1], local_range[0]}};

      cgh.depends_on(memcpy_event);
      cgh.parallel_for(sycl_swapped_range, [functor_wrapper, bare_policy](
                                               sycl::nd_item<3> item) {
        // swap back for correct index calculations in DeviceIterateTile
        const index_type local_x    = item.get_local_id(2);
        const index_type local_y    = item.get_local_id(1);
        const index_type local_z    = item.get_local_id(0);
        const index_type global_x   = item.get_group(2);
        const index_type global_y   = item.get_group(1);
        const index_type global_z   = item.get_group(0);
        const index_type n_global_x = item.get_group_range(2);
        const index_type n_global_y = item.get_group_range(1);
        const index_type n_global_z = item.get_group_range(0);

        Kokkos::Impl::DeviceIterateTile<Policy::rank, BarePolicy, FunctorType,
                                        typename Policy::work_tag>(
            bare_policy, functor_wrapper.get_functor(),
            {n_global_x, n_global_y, n_global_z},
            {global_x, global_y, global_z}, {local_x, local_y, local_z})
            .exec_range();
      });
    });
    q.ext_oneapi_submit_barrier(std::vector<sycl::event>{parallel_for_event});

    return parallel_for_event;
  }

 public:
  using functor_type = FunctorType;

  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy& policy, const Functor&) {
    return policy.space().impl_internal_space_instance()->m_maxWorkgroupSize;
  }

  void execute() const {
    Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem&
        indirectKernelMem =
            m_space.impl_internal_space_instance()->get_indirect_kernel_mem();

    auto functor_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_functor, indirectKernelMem);
    sycl::event event =
        sycl_direct_launch(functor_wrapper, functor_wrapper.get_copy_event());
    functor_wrapper.register_event(event);
  }

  ParallelFor(const ParallelFor&) = delete;
  ParallelFor(ParallelFor&&)      = delete;
  ParallelFor& operator=(const ParallelFor&) = delete;
  ParallelFor& operator=(ParallelFor&&) = delete;
  ~ParallelFor()                        = default;

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_space(arg_policy.space()) {}
};

#endif  // KOKKOS_SYCL_PARALLEL_RANGE_HPP_
