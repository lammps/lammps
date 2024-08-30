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

#ifndef KOKKOS_THREADS_PARALLEL_SCAN_RANGE_HPP
#define KOKKOS_THREADS_PARALLEL_SCAN_RANGE_HPP

#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Threads> {
 private:
  using Policy    = Kokkos::RangePolicy<Traits...>;
  using WorkRange = typename Policy::WorkRange;
  using WorkTag   = typename Policy::work_tag;
  using Member    = typename Policy::member_type;
  using Analysis  = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                         Policy, FunctorType, void>;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update, const bool final) {
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(i, update, final);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update, const bool final) {
    const TagType t{};
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(t, i, update, final);
    }
  }

  static void exec(ThreadsInternal &instance, const void *arg) {
    const ParallelScan &self = *((const ParallelScan *)arg);

    const WorkRange range(self.m_policy, instance.pool_rank(),
                          instance.pool_size());

    typename Analysis::Reducer final_reducer(self.m_functor);

    reference_type update =
        final_reducer.init(static_cast<pointer_type>(instance.reduce_memory()));

    ParallelScan::template exec_range<WorkTag>(self.m_functor, range.begin(),
                                               range.end(), update, false);

    instance.scan_small(final_reducer);

    ParallelScan::template exec_range<WorkTag>(self.m_functor, range.begin(),
                                               range.end(), update, true);

    instance.fan_in();
  }

 public:
  inline void execute() const {
    ThreadsInternal::resize_scratch(2 * Analysis::value_size(m_functor), 0);
    ThreadsInternal::start(&ParallelScan::exec, this);
    ThreadsInternal::fence();
  }

  ParallelScan(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Threads> {
 private:
  using Policy    = Kokkos::RangePolicy<Traits...>;
  using WorkRange = typename Policy::WorkRange;
  using WorkTag   = typename Policy::work_tag;
  using Member    = typename Policy::member_type;

  using Analysis = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                         Policy, FunctorType, ReturnType>;

  using value_type     = typename Analysis::value_type;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;
  const pointer_type m_result_ptr;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update, const bool final) {
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(i, update, final);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update, const bool final) {
    const TagType t{};
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(t, i, update, final);
    }
  }

  static void exec(ThreadsInternal &instance, const void *arg) {
    const ParallelScanWithTotal &self = *((const ParallelScanWithTotal *)arg);

    const WorkRange range(self.m_policy, instance.pool_rank(),
                          instance.pool_size());

    typename Analysis::Reducer final_reducer(self.m_functor);

    reference_type update =
        final_reducer.init(static_cast<pointer_type>(instance.reduce_memory()));

    ParallelScanWithTotal::template exec_range<WorkTag>(
        self.m_functor, range.begin(), range.end(), update, false);

    instance.scan_small(final_reducer);

    ParallelScanWithTotal::template exec_range<WorkTag>(
        self.m_functor, range.begin(), range.end(), update, true);

    instance.fan_in();

    if (instance.pool_rank() == instance.pool_size() - 1) {
      *self.m_result_ptr = update;
    }
  }

 public:
  inline void execute() const {
    ThreadsInternal::resize_scratch(2 * Analysis::value_size(m_functor), 0);
    ThreadsInternal::start(&ParallelScanWithTotal::exec, this);
    ThreadsInternal::fence();
  }

  template <class ViewType>
  ParallelScanWithTotal(const FunctorType &arg_functor,
                        const Policy &arg_policy,
                        const ViewType &arg_result_view)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_result_ptr(arg_result_view.data()) {
    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::Threads parallel_scan result must be host-accessible!");
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
