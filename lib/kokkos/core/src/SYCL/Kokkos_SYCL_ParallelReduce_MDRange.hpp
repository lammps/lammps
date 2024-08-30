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

#ifndef KOKKOS_SYCL_PARALLEL_REDUCE_MDRANGE_HPP
#define KOKKOS_SYCL_PARALLEL_REDUCE_MDRANGE_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Parallel_Reduce.hpp>
#include <SYCL/Kokkos_SYCL_WorkgroupReduction.hpp>
#include <Kokkos_BitManipulation.hpp>

#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
#include <vector>
#endif

template <class CombinedFunctorReducerType, class... Traits>
class Kokkos::Impl::ParallelReduce<CombinedFunctorReducerType,
                                   Kokkos::MDRangePolicy<Traits...>,
                                   Kokkos::Experimental::SYCL> {
 public:
  using Policy      = Kokkos::MDRangePolicy<Traits...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

 private:
  using value_type     = typename ReducerType::value_type;
  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;

  using WorkTag = typename Policy::work_tag;

  // MDRangePolicy is not trivially copyable. Hence, replicate the data we
  // really need in DeviceIterateTile in a trivially copyable struct.
  struct BarePolicy {
    using index_type = typename Policy::index_type;

    BarePolicy(const Policy& policy)
        : m_lower(policy.m_lower),
          m_upper(policy.m_upper),
          m_tile(policy.m_tile),
          m_tile_end(policy.m_tile_end),
          m_num_tiles(policy.m_num_tiles),
          m_prod_tile_dims(policy.m_prod_tile_dims) {}

    const typename Policy::point_type m_lower;
    const typename Policy::point_type m_upper;
    const typename Policy::tile_type m_tile;
    const typename Policy::point_type m_tile_end;
    const typename Policy::index_type m_num_tiles;
    const typename Policy::index_type m_prod_tile_dims;
    static constexpr Iterate inner_direction = Policy::inner_direction;
    static constexpr int rank                = Policy::rank;
  };

 public:
  // V - View
  template <typename View>
  ParallelReduce(const CombinedFunctorReducerType& f, const Policy& p,
                 const View& v)
      : m_functor_reducer(f),
        m_policy(p),
        m_space(p.space()),
        m_result_ptr(v.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::Experimental::SYCLDeviceUSMSpace,
                              typename View::memory_space>::accessible),
        m_scratch_buffers_lock(
            m_space.impl_internal_space_instance()->m_mutexScratchSpace) {}

 private:
  template <typename CombinedFunctorReducerWrapper>
  sycl::event sycl_direct_launch(
      const CombinedFunctorReducerWrapper& functor_reducer_wrapper,
      const sycl::event& memcpy_event) const {
    // Convenience references
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *m_space.impl_internal_space_instance();
    sycl::queue& q = m_space.sycl_queue();

    const typename Policy::index_type n_tiles = m_policy.m_num_tiles;
    const unsigned int value_count =
        m_functor_reducer.get_reducer().value_count();
    sycl::device_ptr<value_type> results_ptr;
    auto host_result_ptr =
        (m_result_ptr && !m_result_ptr_device_accessible)
            ? static_cast<sycl::host_ptr<value_type>>(
                  instance.scratch_host(sizeof(value_type) * value_count))
            : nullptr;

    sycl::event last_reduction_event;

    desul::ensure_sycl_lock_arrays_on_device(q);

    // If n_tiles==0 we only call init() and final() working with the global
    // scratch memory but don't copy back to m_result_ptr yet.
    if (n_tiles == 0) {
      auto parallel_reduce_event = q.submit([&](sycl::handler& cgh) {
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
        cgh.depends_on(memcpy_event);
#else
        (void)memcpy_event;
#endif
        results_ptr = static_cast<sycl::device_ptr<value_type>>(
            instance.scratch_space(sizeof(value_type) * value_count));
        auto device_accessible_result_ptr =
            m_result_ptr_device_accessible
                ? static_cast<sycl::global_ptr<value_type>>(m_result_ptr)
                : static_cast<sycl::global_ptr<value_type>>(host_result_ptr);
        cgh.single_task([=]() {
          const CombinedFunctorReducerType& functor_reducer =
              functor_reducer_wrapper.get_functor();
          const ReducerType& reducer = functor_reducer.get_reducer();
          reducer.init(results_ptr);
          reducer.final(results_ptr);
          if (device_accessible_result_ptr)
            reducer.copy(device_accessible_result_ptr.get(), results_ptr.get());
        });
      });
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
      q.ext_oneapi_submit_barrier(
          std::vector<sycl::event>{parallel_reduce_event});
#endif
      last_reduction_event = parallel_reduce_event;
    } else {
      // Otherwise (when n_tiles is not zero), we perform a reduction on the
      // values in all workgroups separately, write the workgroup results back
      // to global memory and recurse until only one workgroup does the
      // reduction and thus gets the final value.
      const int wgroup_size = Kokkos::bit_ceil(
          static_cast<unsigned int>(m_policy.m_prod_tile_dims));

      // FIXME_SYCL Find a better way to determine a good limit for the
      // maximum number of work groups, also see
      // https://github.com/intel/llvm/blob/756ba2616111235bba073e481b7f1c8004b34ee6/sycl/source/detail/reduction.cpp#L51-L62
      size_t max_work_groups =
          2 * q.get_device().get_info<sycl::info::device::max_compute_units>();
      int values_per_thread = 1;
      size_t n_wgroups      = n_tiles;
      while (n_wgroups > max_work_groups) {
        values_per_thread *= 2;
        n_wgroups = (n_tiles + values_per_thread - 1) / values_per_thread;
      }

      results_ptr = static_cast<sycl::device_ptr<value_type>>(
          instance.scratch_space(sizeof(value_type) * value_count * n_wgroups));
      auto device_accessible_result_ptr =
          m_result_ptr_device_accessible
              ? static_cast<sycl::global_ptr<value_type>>(m_result_ptr)
              : static_cast<sycl::global_ptr<value_type>>(host_result_ptr);
      auto scratch_flags = static_cast<sycl::device_ptr<unsigned int>>(
          instance.scratch_flags(sizeof(unsigned int)));

      auto parallel_reduce_event = q.submit([&](sycl::handler& cgh) {
        sycl::local_accessor<value_type> local_mem(
            sycl::range<1>(wgroup_size) * value_count, cgh);
        sycl::local_accessor<unsigned int> num_teams_done(1, cgh);

        const BarePolicy bare_policy = m_policy;

#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
        cgh.depends_on(memcpy_event);
#else
        (void)memcpy_event;
#endif

        // REMEMBER swap local x<->y to be conforming with Cuda/HIP
        // implementation
        cgh.parallel_for(
            sycl::nd_range<1>{n_wgroups * wgroup_size, wgroup_size},
            [=](sycl::nd_item<1> item) {
              const int local_id = item.get_local_linear_id();
              const CombinedFunctorReducerType& functor_reducer =
                  functor_reducer_wrapper.get_functor();
              const FunctorType& functor = functor_reducer.get_functor();
              const ReducerType& reducer = functor_reducer.get_reducer();

              // In the first iteration, we call functor to initialize the local
              // memory. Otherwise, the local memory is initialized with the
              // results from the previous iteration that are stored in global
              // memory.
              using index_type = typename Policy::index_type;

              // SWAPPED here to be conforming with CUDA implementation
              const index_type local_x    = 0;
              const index_type local_y    = item.get_local_id(0);
              const index_type local_z    = 0;
              const index_type global_y   = 0;
              const index_type global_z   = 0;
              const index_type n_global_x = n_tiles;
              const index_type n_global_y = 1;
              const index_type n_global_z = 1;

              if constexpr (!SYCLReduction::use_shuffle_based_algorithm<
                                ReducerType>) {
                reference_type update =
                    reducer.init(&local_mem[local_id * value_count]);

                for (index_type global_x = item.get_group(0);
                     global_x < n_tiles; global_x += item.get_group_range(0))
                  Kokkos::Impl::Reduce::DeviceIterateTile<
                      Policy::rank, BarePolicy, FunctorType,
                      typename Policy::work_tag, reference_type>(
                      bare_policy, functor, update,
                      {n_global_x, n_global_y, n_global_z},
                      {global_x, global_y, global_z},
                      {local_x, local_y, local_z})
                      .exec_range();
                item.barrier(sycl::access::fence_space::local_space);

                SYCLReduction::workgroup_reduction<>(
                    item, local_mem, results_ptr, device_accessible_result_ptr,
                    value_count, reducer, false, wgroup_size);

                if (local_id == 0) {
                  sycl::atomic_ref<unsigned, sycl::memory_order::acq_rel,
                                   sycl::memory_scope::device,
                                   sycl::access::address_space::global_space>
                      scratch_flags_ref(*scratch_flags);
                  num_teams_done[0] = ++scratch_flags_ref;
                }
                item.barrier(sycl::access::fence_space::local_space);
                if (num_teams_done[0] == n_wgroups) {
                  if (local_id == 0) *scratch_flags = 0;
                  if (local_id >= static_cast<int>(n_wgroups))
                    reducer.init(&local_mem[local_id * value_count]);
                  else {
                    reducer.copy(&local_mem[local_id * value_count],
                                 &results_ptr[local_id * value_count]);
                    for (unsigned int id = local_id + wgroup_size;
                         id < n_wgroups; id += wgroup_size) {
                      reducer.join(&local_mem[local_id * value_count],
                                   &results_ptr[id * value_count]);
                    }
                  }

                  SYCLReduction::workgroup_reduction<>(
                      item, local_mem, results_ptr,
                      device_accessible_result_ptr, value_count, reducer, true,
                      std::min<int>(n_wgroups, wgroup_size));
                }
              } else {
                value_type local_value;
                reference_type update = reducer.init(&local_value);

                for (index_type global_x = item.get_group(0);
                     global_x < n_tiles; global_x += item.get_group_range(0))
                  Kokkos::Impl::Reduce::DeviceIterateTile<
                      Policy::rank, BarePolicy, FunctorType,
                      typename Policy::work_tag, reference_type>(
                      bare_policy, functor, update,
                      {n_global_x, n_global_y, n_global_z},
                      {global_x, global_y, global_z},
                      {local_x, local_y, local_z})
                      .exec_range();

                SYCLReduction::workgroup_reduction<>(
                    item, local_mem, local_value, results_ptr,
                    device_accessible_result_ptr, reducer, false, wgroup_size);

                if (local_id == 0) {
                  sycl::atomic_ref<unsigned, sycl::memory_order::acq_rel,
                                   sycl::memory_scope::device,
                                   sycl::access::address_space::global_space>
                      scratch_flags_ref(*scratch_flags);
                  num_teams_done[0] = ++scratch_flags_ref;
                }
                item.barrier(sycl::access::fence_space::local_space);
                if (num_teams_done[0] == n_wgroups) {
                  if (local_id == 0) *scratch_flags = 0;
                  if (local_id >= static_cast<int>(n_wgroups))
                    reducer.init(&local_value);
                  else {
                    local_value = results_ptr[local_id];
                    for (unsigned int id = local_id + wgroup_size;
                         id < n_wgroups; id += wgroup_size) {
                      reducer.join(&local_value, &results_ptr[id]);
                    }
                  }

                  SYCLReduction::workgroup_reduction<>(
                      item, local_mem, local_value, results_ptr,
                      device_accessible_result_ptr, reducer, true,
                      std::min<int>(n_wgroups, wgroup_size));
                }
              }
            });
      });
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
      q.ext_oneapi_submit_barrier(
          std::vector<sycl::event>{parallel_reduce_event});
#endif
      last_reduction_event = parallel_reduce_event;
    }

    // At this point, the reduced value is written to the entry in results_ptr
    // and all that is left is to copy it back to the given result pointer if
    // necessary.
    // Using DeepCopy instead of fence+memcpy turned out to be up to 2x slower.
    if (host_result_ptr) {
      m_space.fence(
          "Kokkos::Impl::ParallelReduce<SYCL, MDRangePolicy>::execute: result "
          "not device-accessible");
      std::memcpy(m_result_ptr, host_result_ptr,
                  sizeof(value_type) * value_count);
    }

    return last_reduction_event;
  }

 public:
  template <typename Policy, typename Functor>
  static int max_tile_size_product(const Policy& policy, const Functor&) {
    return policy.space().impl_internal_space_instance()->m_maxWorkgroupSize;
  }

  void execute() const {
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *m_space.impl_internal_space_instance();
    using IndirectKernelMem =
        Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem;
    IndirectKernelMem& indirectKernelMem = instance.get_indirect_kernel_mem();

    auto functor_reducer_wrapper =
        Experimental::Impl::make_sycl_function_wrapper(m_functor_reducer,
                                                       indirectKernelMem);

    sycl::event event = sycl_direct_launch(
        functor_reducer_wrapper, functor_reducer_wrapper.get_copy_event());
    functor_reducer_wrapper.register_event(event);
  }

 private:
  const CombinedFunctorReducerType m_functor_reducer;
  const BarePolicy m_policy;
  const Kokkos::Experimental::SYCL& m_space;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;

  // Only let one ParallelReduce instance at a time use the host scratch memory.
  // The constructor acquires the mutex which is released in the destructor.
  std::scoped_lock<std::mutex> m_scratch_buffers_lock;
};

#endif /* KOKKOS_SYCL_PARALLEL_REDUCE_MDRANGE_HPP */
