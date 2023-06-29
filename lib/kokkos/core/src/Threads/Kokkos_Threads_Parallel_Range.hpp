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

#ifndef KOKKOS_THREADS_PARALLEL_RANGE_HPP
#define KOKKOS_THREADS_PARALLEL_RANGE_HPP

#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>,
                  Kokkos::Threads> {
 private:
  using Policy    = Kokkos::RangePolicy<Traits...>;
  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  const FunctorType m_functor;
  const Policy m_policy;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member ibeg, const Member iend) {
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(i);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member ibeg, const Member iend) {
    const TagType t{};
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(t, i);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    exec_schedule<typename Policy::schedule_type::type>(exec, arg);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Static>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    ParallelFor::template exec_range<WorkTag>(self.m_functor, range.begin(),
                                              range.end());

    exec.fan_in();
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    exec.set_work_range(range.begin() - self.m_policy.begin(),
                        range.end() - self.m_policy.begin(),
                        self.m_policy.chunk_size());
    exec.reset_steal_target();
    exec.barrier();

    long work_index = exec.get_work_index();

    while (work_index != -1) {
      const Member begin =
          static_cast<Member>(work_index) * self.m_policy.chunk_size() +
          self.m_policy.begin();
      const Member end =
          begin + self.m_policy.chunk_size() < self.m_policy.end()
              ? begin + self.m_policy.chunk_size()
              : self.m_policy.end();
      ParallelFor::template exec_range<WorkTag>(self.m_functor, begin, end);
      work_index = exec.get_work_index();
    }

    exec.fan_in();
  }

 public:
  inline void execute() const {
    ThreadsExec::start(&ParallelFor::exec, this);
    ThreadsExec::fence();
  }

  ParallelFor(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

template <class CombinedFunctorReducerType, class... Traits>
class ParallelReduce<CombinedFunctorReducerType, Kokkos::RangePolicy<Traits...>,
                     Kokkos::Threads> {
 private:
  using Policy      = Kokkos::RangePolicy<Traits...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag   = typename Policy::work_tag;
  using WorkRange = typename Policy::WorkRange;
  using Member    = typename Policy::member_type;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update) {
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(i, update);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_range(
      const FunctorType &functor, const Member &ibeg, const Member &iend,
      reference_type update) {
    const TagType t{};
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP)
#pragma ivdep
#endif
    for (Member i = ibeg; i < iend; ++i) {
      functor(t, i, update);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    exec_schedule<typename Policy::schedule_type::type>(exec, arg);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Static>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelReduce &self = *((const ParallelReduce *)arg);
    const WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    const ReducerType &reducer = self.m_functor_reducer.get_reducer();

    ParallelReduce::template exec_range<WorkTag>(
        self.m_functor_reducer.get_functor(), range.begin(), range.end(),
        reducer.init(static_cast<pointer_type>(exec.reduce_memory())));

    exec.fan_in_reduce(reducer);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_schedule(ThreadsExec &exec, const void *arg) {
    const ParallelReduce &self = *((const ParallelReduce *)arg);
    const WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    exec.set_work_range(range.begin() - self.m_policy.begin(),
                        range.end() - self.m_policy.begin(),
                        self.m_policy.chunk_size());
    exec.reset_steal_target();
    exec.barrier();

    long work_index            = exec.get_work_index();
    const ReducerType &reducer = self.m_functor_reducer.get_reducer();

    reference_type update =
        reducer.init(static_cast<pointer_type>(exec.reduce_memory()));
    while (work_index != -1) {
      const Member begin =
          static_cast<Member>(work_index) * self.m_policy.chunk_size() +
          self.m_policy.begin();
      const Member end =
          begin + self.m_policy.chunk_size() < self.m_policy.end()
              ? begin + self.m_policy.chunk_size()
              : self.m_policy.end();
      ParallelReduce::template exec_range<WorkTag>(
          self.m_functor_reducer.get_functor(), begin, end, update);
      work_index = exec.get_work_index();
    }

    exec.fan_in_reduce(reducer);
  }

 public:
  inline void execute() const {
    const ReducerType &reducer = m_functor_reducer.get_reducer();

    if (m_policy.end() <= m_policy.begin()) {
      if (m_result_ptr) {
        reducer.init(m_result_ptr);
        reducer.final(m_result_ptr);
      }
    } else {
      ThreadsExec::resize_scratch(reducer.value_size(), 0);

      ThreadsExec::start(&ParallelReduce::exec, this);

      ThreadsExec::fence();

      if (m_result_ptr) {
        const pointer_type data =
            (pointer_type)ThreadsExec::root_reduce_scratch();

        const unsigned n = reducer.value_count();
        for (unsigned i = 0; i < n; ++i) {
          m_result_ptr[i] = data[i];
        }
      }
    }
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType &arg_functor_reducer,
                 const Policy &arg_policy, const ViewType &arg_result_view)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result_view.data()) {
    static_assert(Kokkos::is_view<ViewType>::value,
                  "Kokkos::Threads reduce result must be a View");

    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::Threads reduce result must be a View accessible from "
        "HostSpace");
  }
};

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

  static void exec(ThreadsExec &exec, const void *arg) {
    const ParallelScan &self = *((const ParallelScan *)arg);

    const WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    typename Analysis::Reducer final_reducer(self.m_functor);

    reference_type update =
        final_reducer.init(static_cast<pointer_type>(exec.reduce_memory()));

    ParallelScan::template exec_range<WorkTag>(self.m_functor, range.begin(),
                                               range.end(), update, false);

    //  exec.template scan_large( final_reducer );
    exec.scan_small(final_reducer);

    ParallelScan::template exec_range<WorkTag>(self.m_functor, range.begin(),
                                               range.end(), update, true);

    exec.fan_in();
  }

 public:
  inline void execute() const {
    ThreadsExec::resize_scratch(2 * Analysis::value_size(m_functor), 0);
    ThreadsExec::start(&ParallelScan::exec, this);
    ThreadsExec::fence();
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

  static void exec(ThreadsExec &exec, const void *arg) {
    const ParallelScanWithTotal &self = *((const ParallelScanWithTotal *)arg);

    const WorkRange range(self.m_policy, exec.pool_rank(), exec.pool_size());

    typename Analysis::Reducer final_reducer(self.m_functor);

    reference_type update =
        final_reducer.init(static_cast<pointer_type>(exec.reduce_memory()));

    ParallelScanWithTotal::template exec_range<WorkTag>(
        self.m_functor, range.begin(), range.end(), update, false);

    //  exec.template scan_large(final_reducer);
    exec.scan_small(final_reducer);

    ParallelScanWithTotal::template exec_range<WorkTag>(
        self.m_functor, range.begin(), range.end(), update, true);

    exec.fan_in();

    if (exec.pool_rank() == exec.pool_size() - 1) {
      *self.m_result_ptr = update;
    }
  }

 public:
  inline void execute() const {
    ThreadsExec::resize_scratch(2 * Analysis::value_size(m_functor), 0);
    ThreadsExec::start(&ParallelScanWithTotal::exec, this);
    ThreadsExec::fence();
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
