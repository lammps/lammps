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

#ifndef KOKKOS_SYCL_PARALLEL_FOR_TEAM_HPP
#define KOKKOS_SYCL_PARALLEL_FOR_TEAM_HPP

#include <Kokkos_Parallel.hpp>

#include <SYCL/Kokkos_SYCL_Team.hpp>
#include <SYCL/Kokkos_SYCL_TeamPolicy.hpp>

#include <sstream>
#include <vector>

template <typename FunctorType, typename... Properties>
class Kokkos::Impl::ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                                Kokkos::Experimental::SYCL> {
 public:
  using Policy       = TeamPolicy<Properties...>;
  using functor_type = FunctorType;
  using size_type    = ::Kokkos::Experimental::SYCL::size_type;

 private:
  using member_type   = typename Policy::member_type;
  using work_tag      = typename Policy::work_tag;
  using launch_bounds = typename Policy::launch_bounds;

  FunctorType const m_functor;
  Policy const m_policy;
  size_type const m_league_size;
  int m_team_size;
  size_type const m_vector_size;
  int m_shmem_begin;
  int m_shmem_size;
  size_t m_scratch_size[2];

  template <typename FunctorWrapper>
  sycl::event sycl_direct_launch(const sycl_device_ptr<char> global_scratch_ptr,
                                 const FunctorWrapper& functor_wrapper,
                                 const sycl::event& memcpy_event) const {
    // Convenience references
    const Kokkos::Experimental::SYCL& space = m_policy.space();
    sycl::queue& q                          = space.sycl_queue();

    desul::ensure_sycl_lock_arrays_on_device(q);

    auto cgh_lambda = [&](sycl::handler& cgh) {
      // FIXME_SYCL accessors seem to need a size greater than zero at least for
      // host queues
      sycl::local_accessor<char, 1> team_scratch_memory_L0(
          sycl::range<1>(
              std::max(m_scratch_size[0] + m_shmem_begin, size_t(1))),
          cgh);

      // Avoid capturing *this since it might not be trivially copyable
      const auto shmem_begin       = m_shmem_begin;
      const size_t scratch_size[2] = {m_scratch_size[0], m_scratch_size[1]};

      auto lambda = [=](sycl::nd_item<2> item) {
        const member_type team_member(
            KOKKOS_IMPL_SYCL_GET_MULTI_PTR(team_scratch_memory_L0), shmem_begin,
            scratch_size[0],
            global_scratch_ptr + item.get_group(1) * scratch_size[1],
            scratch_size[1], item, item.get_group_linear_id(),
            item.get_group_range(1));
        if constexpr (std::is_void<work_tag>::value)
          functor_wrapper.get_functor()(team_member);
        else
          functor_wrapper.get_functor()(work_tag(), team_member);
      };

      static sycl::kernel kernel = [&] {
        sycl::kernel_id functor_kernel_id =
            sycl::get_kernel_id<decltype(lambda)>();
        auto kernel_bundle =
            sycl::get_kernel_bundle<sycl::bundle_state::executable>(
                q.get_context(), std::vector{functor_kernel_id});
        return kernel_bundle.get_kernel(functor_kernel_id);
      }();
      auto max_sg_size =
          kernel
              .get_info<sycl::info::kernel_device_specific::max_sub_group_size>(
                  q.get_device());
      auto final_vector_size = std::min<int>(m_vector_size, max_sg_size);
      // FIXME_SYCL For some reason, explicitly enforcing the kernel bundle to
      // be used gives a runtime error.
      // cgh.use_kernel_bundle(kernel_bundle);

#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
      cgh.depends_on(memcpy_event);
#else
      (void)memcpy_event;
#endif
      cgh.parallel_for(
          sycl::nd_range<2>(
              sycl::range<2>(m_team_size, m_league_size * final_vector_size),
              sycl::range<2>(m_team_size, final_vector_size)),
          lambda);
    };

#ifdef SYCL_EXT_ONEAPI_GRAPH
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

 public:
  inline void execute() const {
    if (m_league_size == 0) return;

    auto& instance = *m_policy.space().impl_internal_space_instance();

    // Only let one instance at a time resize the instance's scratch memory
    // allocations.
    std::scoped_lock<std::mutex> team_scratch_lock(
        instance.m_team_scratch_mutex);

    // Functor's reduce memory, team scan memory, and team shared memory depend
    // upon team size.
    int scratch_pool_id = instance.acquire_team_scratch_space();
    const sycl_device_ptr<char> global_scratch_ptr =
        static_cast<sycl_device_ptr<char>>(instance.resize_team_scratch_space(
            scratch_pool_id,
            static_cast<ptrdiff_t>(m_scratch_size[1]) * m_league_size));

    Kokkos::Experimental::Impl::SYCLInternal::IndirectKernelMem&
        indirectKernelMem = instance.get_indirect_kernel_mem();

    auto functor_wrapper = Experimental::Impl::make_sycl_function_wrapper(
        m_functor, indirectKernelMem);

    sycl::event event = sycl_direct_launch(global_scratch_ptr, functor_wrapper,
                                           functor_wrapper.get_copy_event());
    functor_wrapper.register_event(event);
    instance.register_team_scratch_event(scratch_pool_id, event);
  }

  ParallelFor(FunctorType const& arg_functor, Policy const& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_league_size(arg_policy.league_size()),
        m_team_size(arg_policy.team_size()),
        m_vector_size(arg_policy.impl_vector_length()) {
    // FIXME_SYCL optimize
    if (m_team_size < 0)
      m_team_size =
          m_policy.team_size_recommended(arg_functor, ParallelForTag{});

    m_shmem_begin = (sizeof(double) * (m_team_size + 2));
    m_shmem_size =
        (m_policy.scratch_size(0, m_team_size) +
         FunctorTeamShmemSize<FunctorType>::value(m_functor, m_team_size));
    m_scratch_size[0] = m_shmem_size;
    m_scratch_size[1] = m_policy.scratch_size(1, m_team_size);

    const auto& instance = *m_policy.space().impl_internal_space_instance();
    if (static_cast<int>(instance.m_maxShmemPerBlock) <
        m_shmem_size - m_shmem_begin) {
      std::stringstream out;
      out << "Kokkos::Impl::ParallelFor<SYCL> insufficient shared memory! "
             "Requested "
          << m_shmem_size - m_shmem_begin << " bytes but maximum is "
          << instance.m_maxShmemPerBlock << '\n';
      Kokkos::Impl::throw_runtime_exception(out.str());
    }

    const auto max_team_size =
        m_policy.team_size_max(arg_functor, ParallelForTag{});
    if (m_team_size > m_policy.team_size_max(arg_functor, ParallelForTag{}))
      Kokkos::Impl::throw_runtime_exception(
          "Kokkos::Impl::ParallelFor<SYCL> requested too large team size. The "
          "maximal team_size is " +
          std::to_string(max_team_size) + '!');
  }
};

#endif
