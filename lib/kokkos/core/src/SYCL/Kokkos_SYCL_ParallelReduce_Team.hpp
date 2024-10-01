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

#ifndef KOKKOS_SYCL_PARALLEL_REDUCE_TEAM_HPP
#define KOKKOS_SYCL_PARALLEL_REDUCE_TEAM_HPP

#include <Kokkos_Parallel.hpp>

#include <SYCL/Kokkos_SYCL_Team.hpp>
#include <SYCL/Kokkos_SYCL_TeamPolicy.hpp>
#include <SYCL/Kokkos_SYCL_WorkgroupReduction.hpp>

#include <sstream>
#include <vector>

template <class CombinedFunctorReducerType, class... Properties>
class Kokkos::Impl::ParallelReduce<CombinedFunctorReducerType,
                                   Kokkos::TeamPolicy<Properties...>,
                                   Kokkos::Experimental::SYCL> {
 public:
  using Policy      = TeamPolicy<Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

 private:
  using member_type   = typename Policy::member_type;
  using WorkTag       = typename Policy::work_tag;
  using launch_bounds = typename Policy::launch_bounds;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;
  using value_type     = typename ReducerType::value_type;

 public:
  using functor_type = FunctorType;
  using size_type    = Kokkos::Experimental::SYCL::size_type;

 private:
  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const bool m_result_ptr_device_accessible;
  size_type m_shmem_begin;
  size_type m_shmem_size;
  size_t m_scratch_size[2];
  const size_type m_league_size;
  int m_team_size;
  const size_type m_vector_size;

  template <typename CombinedFunctorReducerWrapper>
  sycl::event sycl_direct_launch(
      const sycl_device_ptr<char> global_scratch_ptr,
      const CombinedFunctorReducerWrapper& functor_reducer_wrapper,
      const sycl::event& memcpy_event) const {
    // Convenience references
    const Kokkos::Experimental::SYCL& space = m_policy.space();
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *space.impl_internal_space_instance();
    sycl::queue& q = space.sycl_queue();

    const unsigned int value_count =
        m_functor_reducer.get_reducer().value_count();
    std::size_t size = std::size_t(m_league_size) * m_team_size * m_vector_size;
    value_type* results_ptr = nullptr;
    auto host_result_ptr =
        (m_result_ptr && !m_result_ptr_device_accessible)
            ? static_cast<sycl_host_ptr<value_type>>(
                  instance.scratch_host(sizeof(value_type) * value_count))
            : nullptr;

    sycl::event last_reduction_event;

    desul::ensure_sycl_lock_arrays_on_device(q);

    // If size<=1 we only call init(), the functor and possibly final once
    // working with the global scratch memory but don't copy back to
    // m_result_ptr yet.
    if (size <= 1) {
      results_ptr =
          static_cast<sycl_device_ptr<value_type>>(instance.scratch_space(
              sizeof(value_type) * std::max(value_count, 1u)));
      auto device_accessible_result_ptr =
          m_result_ptr_device_accessible
              ? static_cast<sycl::global_ptr<value_type>>(m_result_ptr)
              : static_cast<sycl::global_ptr<value_type>>(host_result_ptr);

      auto cgh_lambda = [&](sycl::handler& cgh) {
        // FIXME_SYCL accessors seem to need a size greater than zero at least
        // for host queues
        sycl::local_accessor<char, 1> team_scratch_memory_L0(
            sycl::range<1>(
                std::max(m_scratch_size[0] + m_shmem_begin, size_t(1))),
            cgh);

        // Avoid capturing *this since it might not be trivially copyable
        const auto shmem_begin       = m_shmem_begin;
        const size_t scratch_size[2] = {m_scratch_size[0], m_scratch_size[1]};

#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
        cgh.depends_on(memcpy_event);
#else
        (void)memcpy_event;
#endif
        cgh.parallel_for(
            sycl::nd_range<2>(sycl::range<2>(1, 1), sycl::range<2>(1, 1)),
            [=](sycl::nd_item<2> item) {
              const CombinedFunctorReducerType& functor_reducer =
                  functor_reducer_wrapper.get_functor();
              const FunctorType& functor = functor_reducer.get_functor();
              const ReducerType& reducer = functor_reducer.get_reducer();

              reference_type update = reducer.init(results_ptr);
              if (size == 1) {
                const member_type team_member(
                    KOKKOS_IMPL_SYCL_GET_MULTI_PTR(team_scratch_memory_L0),
                    shmem_begin, scratch_size[0], global_scratch_ptr,
                    scratch_size[1], item, item.get_group_linear_id(),
                    item.get_group_range(1));
                if constexpr (std::is_void_v<WorkTag>)
                  functor(team_member, update);
                else
                  functor(WorkTag(), team_member, update);
              }
              reducer.final(results_ptr);
              if (device_accessible_result_ptr)
                reducer.copy(device_accessible_result_ptr, &results_ptr[0]);
            });
      };
#ifdef SYCL_EXT_ONEAPI_GRAPH
      if constexpr (Policy::is_graph_kernel::value) {
        sycl_attach_kernel_to_node(*this, cgh_lambda);
      } else
#endif
      {
        last_reduction_event = q.submit(cgh_lambda);
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
        q.ext_oneapi_submit_barrier(
            std::vector<sycl::event>{last_reduction_event});
#endif
      }
    } else {
      // Otherwise, (if the total range has more than one element) we perform a
      // reduction on the values in all workgroups separately, write the
      // workgroup results back to global memory and recurse until only one
      // workgroup does the reduction and thus gets the final value.
      auto cgh_lambda = [&](sycl::handler& cgh) {
        auto scratch_flags = static_cast<sycl_device_ptr<unsigned int>>(
            instance.scratch_flags(sizeof(unsigned int)));

        // FIXME_SYCL accessors seem to need a size greater than zero at least
        // for host queues
        sycl::local_accessor<char, 1> team_scratch_memory_L0(
            sycl::range<1>(
                std::max(m_scratch_size[0] + m_shmem_begin, size_t(1))),
            cgh);

        // Avoid capturing *this since it might not be trivially copyable
        const auto shmem_begin       = m_shmem_begin;
        const auto league_size       = m_league_size;
        const size_t scratch_size[2] = {m_scratch_size[0], m_scratch_size[1]};
        sycl::local_accessor<unsigned int> num_teams_done(1, cgh);

        auto team_reduction_factory =
            [&](sycl::local_accessor<value_type, 1> local_mem,
                sycl_device_ptr<value_type> results_ptr) {
              auto device_accessible_result_ptr =
                  m_result_ptr_device_accessible
                      ? static_cast<sycl::global_ptr<value_type>>(m_result_ptr)
                      : static_cast<sycl::global_ptr<value_type>>(
                            host_result_ptr);
              auto lambda = [=](sycl::nd_item<2> item) {
                auto n_wgroups = item.get_group_range()[1];
                int wgroup_size =
                    item.get_local_range()[0] * item.get_local_range()[1];
                auto group_id = item.get_group_linear_id();
                auto size     = n_wgroups * wgroup_size;

                const auto local_id = item.get_local_linear_id();
                const CombinedFunctorReducerType& functor_reducer =
                    functor_reducer_wrapper.get_functor();
                const FunctorType& functor = functor_reducer.get_functor();
                const ReducerType& reducer = functor_reducer.get_reducer();

                if constexpr (!SYCLReduction::use_shuffle_based_algorithm<
                                  ReducerType>) {
                  reference_type update =
                      reducer.init(&local_mem[local_id * value_count]);
                  for (int league_rank = group_id; league_rank < league_size;
                       league_rank += n_wgroups) {
                    const member_type team_member(
                        KOKKOS_IMPL_SYCL_GET_MULTI_PTR(team_scratch_memory_L0),
                        shmem_begin, scratch_size[0],
                        global_scratch_ptr +
                            item.get_group(1) * scratch_size[1],
                        scratch_size[1], item, league_rank, league_size);
                    if constexpr (std::is_void_v<WorkTag>)
                      functor(team_member, update);
                    else
                      functor(WorkTag(), team_member, update);
                  }
                  item.barrier(sycl::access::fence_space::local_space);

                  SYCLReduction::workgroup_reduction<>(
                      item, local_mem, results_ptr,
                      device_accessible_result_ptr, value_count, reducer, false,
                      std::min<std::size_t>(size,
                                            item.get_local_range()[0] *
                                                item.get_local_range()[1]));

                  if (local_id == 0) {
                    sycl::atomic_ref<unsigned, sycl::memory_order::acq_rel,
                                     sycl::memory_scope::device,
                                     sycl::access::address_space::global_space>
                        scratch_flags_ref(*scratch_flags);
                    num_teams_done[0] = ++scratch_flags_ref;
                  }
                  sycl::group_barrier(item.get_group());
                  if (num_teams_done[0] == n_wgroups) {
                    if (local_id == 0) *scratch_flags = 0;
                    if (local_id >= n_wgroups)
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
                        device_accessible_result_ptr, value_count, reducer,
                        true,
                        std::min(n_wgroups, item.get_local_range()[0] *
                                                item.get_local_range()[1]));
                  }
                } else {
                  value_type local_value;
                  reference_type update = reducer.init(&local_value);
                  for (int league_rank = group_id; league_rank < league_size;
                       league_rank += n_wgroups) {
                    const member_type team_member(
                        KOKKOS_IMPL_SYCL_GET_MULTI_PTR(team_scratch_memory_L0),
                        shmem_begin, scratch_size[0],
                        global_scratch_ptr +
                            item.get_group(1) * scratch_size[1],
                        scratch_size[1], item, league_rank, league_size);
                    if constexpr (std::is_void_v<WorkTag>)
                      functor(team_member, update);
                    else
                      functor(WorkTag(), team_member, update);
                  }

                  SYCLReduction::workgroup_reduction<>(
                      item, local_mem, local_value, results_ptr,
                      device_accessible_result_ptr, reducer, false,
                      std::min<std::size_t>(size,
                                            item.get_local_range()[0] *
                                                item.get_local_range()[1]));

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
                    if (local_id >= n_wgroups)
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
                        std::min(n_wgroups, item.get_local_range()[0] *
                                                item.get_local_range()[1]));
                  }
                }
              };
              return lambda;
            };

        auto dummy_reduction_lambda = team_reduction_factory({1, cgh}, nullptr);

        static sycl::kernel kernel = [&] {
          sycl::kernel_id functor_kernel_id =
              sycl::get_kernel_id<decltype(dummy_reduction_lambda)>();
          auto kernel_bundle =
              sycl::get_kernel_bundle<sycl::bundle_state::executable>(
                  q.get_context(), std::vector{functor_kernel_id});
          return kernel_bundle.get_kernel(functor_kernel_id);
        }();
        auto max_sg_size = kernel.get_info<
            sycl::info::kernel_device_specific::max_sub_group_size>(
            q.get_device());
        auto final_vector_size = std::min<int>(m_vector_size, max_sg_size);
        // FIXME_SYCL For some reason, explicitly enforcing the kernel bundle to
        // be used gives a runtime error.

        //     cgh.use_kernel_bundle(kernel_bundle);

        auto wgroup_size = m_team_size * final_vector_size;
        std::size_t size = std::size_t(m_league_size) * wgroup_size;
        sycl::local_accessor<value_type, 1> local_mem(
            sycl::range<1>(wgroup_size) * std::max(value_count, 1u), cgh);

        const auto init_size =
            std::max<std::size_t>((size + wgroup_size - 1) / wgroup_size, 1);
        results_ptr =
            static_cast<sycl_device_ptr<value_type>>(instance.scratch_space(
                sizeof(value_type) * std::max(value_count, 1u) * init_size));

        size_t max_work_groups =
            2 *
            q.get_device().get_info<sycl::info::device::max_compute_units>();
        int values_per_thread = 1;
        size_t n_wgroups      = m_league_size;
        while (n_wgroups > max_work_groups) {
          values_per_thread *= 2;
          n_wgroups =
              ((size_t(m_league_size) * wgroup_size + values_per_thread - 1) /
                   values_per_thread +
               wgroup_size - 1) /
              wgroup_size;
        }

        auto reduction_lambda = team_reduction_factory(local_mem, results_ptr);

#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
        cgh.depends_on(memcpy_event);
#endif

        cgh.parallel_for(
            sycl::nd_range<2>(
                sycl::range<2>(m_team_size, n_wgroups * m_vector_size),
                sycl::range<2>(m_team_size, m_vector_size)),
            reduction_lambda);
      };
#ifdef SYCL_EXT_ONEAPI_GRAPH
      if constexpr (Policy::is_graph_kernel::value) {
        sycl_attach_kernel_to_node(*this, cgh_lambda);
      } else
#endif
      {
        last_reduction_event = q.submit(cgh_lambda);
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
        q.ext_oneapi_submit_barrier(
            std::vector<sycl::event>{last_reduction_event});
#endif
      }
    }

    // At this point, the reduced value is written to the entry in results_ptr
    // and all that is left is to copy it back to the given result pointer if
    // necessary.
    // Using DeepCopy instead of fence+memcpy turned out to be up to 2x slower.
    if (host_result_ptr) {
      if constexpr (Policy::is_graph_kernel::value)
        Kokkos::abort(
            "parallel_reduce not implemented for graph kernels if result is "
            "not device-accessible!");

      space.fence(
          "Kokkos::Impl::ParallelReduce<SYCL, TeamPolicy>::execute: result not "
          "device-accessible");
      std::memcpy(m_result_ptr, host_result_ptr,
                  sizeof(*m_result_ptr) * value_count);
    }

    return last_reduction_event;
  }

 public:
  inline void execute() {
    Kokkos::Experimental::Impl::SYCLInternal& instance =
        *m_policy.space().impl_internal_space_instance();

    // Only let one instance at a time resize the instance's scratch memory
    // allocations.
    std::scoped_lock<std::mutex> scratch_buffers_lock(
        instance.m_mutexScratchSpace);
    std::scoped_lock<std::mutex> team_scratch_lock(
        instance.m_team_scratch_mutex);

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.
    int scratch_pool_id = instance.acquire_team_scratch_space();
    const sycl_device_ptr<char> global_scratch_ptr =
        static_cast<sycl_device_ptr<char>>(instance.resize_team_scratch_space(
            scratch_pool_id,
            static_cast<ptrdiff_t>(m_scratch_size[1]) * m_league_size));

    using IndirectKernelMem =
        Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem;
    IndirectKernelMem& indirectKernelMem = instance.get_indirect_kernel_mem();

    auto functor_reducer_wrapper =
        Experimental::Impl::make_sycl_function_wrapper(m_functor_reducer,
                                                       indirectKernelMem);

    sycl::event event =
        sycl_direct_launch(global_scratch_ptr, functor_reducer_wrapper,
                           functor_reducer_wrapper.get_copy_event());
    functor_reducer_wrapper.register_event(event);
    instance.register_team_scratch_event(scratch_pool_id, event);
  }

  template <class ViewType>
  ParallelReduce(CombinedFunctorReducerType const& arg_functor_reducer,
                 Policy const& arg_policy, ViewType const& arg_result)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result.data()),
        m_result_ptr_device_accessible(
            MemorySpaceAccess<Kokkos::Experimental::SYCLDeviceUSMSpace,
                              typename ViewType::memory_space>::accessible),
        m_league_size(arg_policy.league_size()),
        m_team_size(arg_policy.team_size()),
        m_vector_size(arg_policy.impl_vector_length()) {
    // FIXME_SYCL optimize
    if (m_team_size < 0)
      m_team_size = m_policy.team_size_recommended(
          m_functor_reducer.get_functor(), m_functor_reducer.get_reducer(),
          ParallelReduceTag{});
    // Must be a power of two greater than two, get the one not bigger than the
    // requested one.
    if ((m_team_size & m_team_size - 1) || m_team_size < 2) {
      int temp_team_size = 2;
      while ((temp_team_size << 1) < m_team_size) temp_team_size <<= 1;
      m_team_size = temp_team_size;
    }

    m_shmem_begin     = (sizeof(double) * (m_team_size + 2));
    m_shmem_size      = (m_policy.scratch_size(0, m_team_size) +
                    FunctorTeamShmemSize<FunctorType>::value(
                        m_functor_reducer.get_functor(), m_team_size));
    m_scratch_size[0] = m_shmem_size;
    m_scratch_size[1] = m_policy.scratch_size(1, m_team_size);

    const Kokkos::Experimental::Impl::SYCLInternal& instance =
        *m_policy.space().impl_internal_space_instance();
    if (static_cast<int>(instance.m_maxShmemPerBlock) <
        m_shmem_size - m_shmem_begin) {
      std::stringstream out;
      out << "Kokkos::Impl::ParallelFor<SYCL> insufficient shared memory! "
             "Requested "
          << m_shmem_size - m_shmem_begin << " bytes but maximum is "
          << instance.m_maxShmemPerBlock << '\n';
      Kokkos::Impl::throw_runtime_exception(out.str());
    }

    if (m_team_size > m_policy.team_size_max(m_functor_reducer.get_functor(),
                                             m_functor_reducer.get_reducer(),
                                             ParallelReduceTag{}))
      Kokkos::Impl::throw_runtime_exception(
          "Kokkos::Impl::ParallelFor<SYCL> requested too large team size.");
  }
};

#endif
