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

#ifndef KOKKOS_SYCL_PARALLEL_FOR_RANGE_HPP_
#define KOKKOS_SYCL_PARALLEL_FOR_RANGE_HPP_

#ifdef SYCL_EXT_ONEAPI_AUTO_LOCAL_RANGE
#include <Kokkos_BitManipulation.hpp>
#endif
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
#include <vector>
#endif

namespace Kokkos::Impl {
#ifndef SYCL_EXT_ONEAPI_AUTO_LOCAL_RANGE
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
#endif

// Same as above but for a user-provided workgroup size
template <typename FunctorWrapper, typename Policy>
struct FunctorWrapperRangePolicyParallelForCustom {
  using WorkTag = typename Policy::work_tag;

  void operator()(sycl::nd_item<1> item) const {
    const typename Policy::index_type id = item.get_global_linear_id();
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
  sycl::event sycl_direct_launch(const Policy& policy, const Functor& functor,
                                 const sycl::event& memcpy_event) const {
    // Convenience references
    const Kokkos::Experimental::SYCL& space = policy.space();
    sycl::queue& q                          = space.sycl_queue();

    desul::ensure_sycl_lock_arrays_on_device(q);

    auto cgh_lambda = [&](sycl::handler& cgh) {
#ifndef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
      cgh.depends_on(memcpy_event);
#else
      (void)memcpy_event;
#endif

      if (policy.chunk_size() <= 1) {
#ifdef SYCL_EXT_ONEAPI_AUTO_LOCAL_RANGE
        const auto actual_range = policy.end() - policy.begin();
        FunctorWrapperRangePolicyParallelForCustom<Functor, Policy> f{
            policy.begin(), functor, actual_range};
        // Round the actual range up to the closest power of two not exceeding
        // the maximum workgroup size
        const auto max_wgroup_size =
            q.get_device().get_info<sycl::info::device::max_work_group_size>();
        const auto wgroup_size_multiple = Kokkos::bit_floor(
            std::min<std::size_t>(max_wgroup_size, actual_range));

        const auto launch_range = (actual_range + wgroup_size_multiple - 1) /
                                  wgroup_size_multiple * wgroup_size_multiple;
        sycl::nd_range<1> range(
            launch_range, sycl::ext::oneapi::experimental::auto_range<1>());
        cgh.parallel_for<
            FunctorWrapperRangePolicyParallelForCustom<Functor, Policy>>(range,
                                                                         f);
#else
        FunctorWrapperRangePolicyParallelFor<Functor, Policy> f{policy.begin(),
                                                                functor};
        sycl::range<1> range(policy.end() - policy.begin());
        cgh.parallel_for<FunctorWrapperRangePolicyParallelFor<Functor, Policy>>(
            range, f);
#endif
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

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

#endif  // KOKKOS_SYCL_PARALLEL_FOR_RANGE_HPP_
