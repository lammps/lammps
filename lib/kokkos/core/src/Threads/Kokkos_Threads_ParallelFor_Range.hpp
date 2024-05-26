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

#ifndef KOKKOS_THREADS_PARALLEL_FOR_RANGE_HPP
#define KOKKOS_THREADS_PARALLEL_FOR_RANGE_HPP

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

  static void exec(ThreadsInternal &instance, const void *arg) {
    exec_schedule<typename Policy::schedule_type::type>(instance, arg);
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Static>::value>
  exec_schedule(ThreadsInternal &instance, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    WorkRange range(self.m_policy, instance.pool_rank(), instance.pool_size());

    ParallelFor::template exec_range<WorkTag>(self.m_functor, range.begin(),
                                              range.end());

    instance.fan_in();
  }

  template <class Schedule>
  static std::enable_if_t<std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_schedule(ThreadsInternal &instance, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    WorkRange range(self.m_policy, instance.pool_rank(), instance.pool_size());

    instance.set_work_range(range.begin() - self.m_policy.begin(),
                            range.end() - self.m_policy.begin(),
                            self.m_policy.chunk_size());
    instance.reset_steal_target();
    instance.barrier();

    long work_index = instance.get_work_index();

    while (work_index != -1) {
      const Member begin =
          static_cast<Member>(work_index) * self.m_policy.chunk_size() +
          self.m_policy.begin();
      const Member end =
          begin + self.m_policy.chunk_size() < self.m_policy.end()
              ? begin + self.m_policy.chunk_size()
              : self.m_policy.end();
      ParallelFor::template exec_range<WorkTag>(self.m_functor, begin, end);
      work_index = instance.get_work_index();
    }

    instance.fan_in();
  }

 public:
  inline void execute() const {
    ThreadsInternal::start(&ParallelFor::exec, this);
    ThreadsInternal::fence();
  }

  ParallelFor(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

}  // namespace Impl
}  // namespace Kokkos

#endif
