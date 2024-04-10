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

#ifndef KOKKOS_THREADS_PARALLEL_REDUCE_RANGE_HPP
#define KOKKOS_THREADS_PARALLEL_REDUCE_RANGE_HPP

#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {

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

  static void exec(ThreadsInternal &instance, const void *arg) {
    exec_schedule<typename Policy::schedule_type::type>(instance, arg);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Static>::value>
  exec_schedule(ThreadsInternal &instance, const void *arg) {
    const ParallelReduce &self = *((const ParallelReduce *)arg);
    const WorkRange range(self.m_policy, instance.pool_rank(),
                          instance.pool_size());

    const ReducerType &reducer = self.m_functor_reducer.get_reducer();

    ParallelReduce::template exec_range<WorkTag>(
        self.m_functor_reducer.get_functor(), range.begin(), range.end(),
        reducer.init(static_cast<pointer_type>(instance.reduce_memory())));

    instance.fan_in_reduce(reducer);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_schedule(ThreadsInternal &instance, const void *arg) {
    const ParallelReduce &self = *((const ParallelReduce *)arg);
    const WorkRange range(self.m_policy, instance.pool_rank(),
                          instance.pool_size());

    instance.set_work_range(range.begin() - self.m_policy.begin(),
                            range.end() - self.m_policy.begin(),
                            self.m_policy.chunk_size());
    instance.reset_steal_target();
    instance.barrier();

    long work_index            = instance.get_work_index();
    const ReducerType &reducer = self.m_functor_reducer.get_reducer();

    reference_type update =
        reducer.init(static_cast<pointer_type>(instance.reduce_memory()));
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
      work_index = instance.get_work_index();
    }

    instance.fan_in_reduce(reducer);
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
      ThreadsInternal::resize_scratch(reducer.value_size(), 0);

      ThreadsInternal::start(&ParallelReduce::exec, this);

      ThreadsInternal::fence();

      if (m_result_ptr) {
        const pointer_type data =
            (pointer_type)ThreadsInternal::root_reduce_scratch();

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

}  // namespace Impl
}  // namespace Kokkos

#endif
