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

#ifndef KOKKOS_THREADS_PARALLEL_TEAM_HPP
#define KOKKOS_THREADS_PARALLEL_TEAM_HPP

#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {

template <class FunctorType, class... Properties>
class ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                  Kokkos::Threads> {
 private:
  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::Threads, Properties...>;
  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  const FunctorType m_functor;
  const Policy m_policy;
  const size_t m_shared;

  template <class TagType, class Schedule>
  inline static std::enable_if_t<std::is_void<TagType>::value &&
                                 std::is_same<Schedule, Kokkos::Static>::value>
  exec_team(const FunctorType &functor, Member member) {
    for (; member.valid_static(); member.next_static()) {
      functor(member);
    }
  }

  template <class TagType, class Schedule>
  inline static std::enable_if_t<!std::is_void<TagType>::value &&
                                 std::is_same<Schedule, Kokkos::Static>::value>
  exec_team(const FunctorType &functor, Member member) {
    const TagType t{};
    for (; member.valid_static(); member.next_static()) {
      functor(t, member);
    }
  }

  template <class TagType, class Schedule>
  inline static std::enable_if_t<std::is_void<TagType>::value &&
                                 std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_team(const FunctorType &functor, Member member) {
    for (; member.valid_dynamic(); member.next_dynamic()) {
      functor(member);
    }
  }

  template <class TagType, class Schedule>
  inline static std::enable_if_t<!std::is_void<TagType>::value &&
                                 std::is_same<Schedule, Kokkos::Dynamic>::value>
  exec_team(const FunctorType &functor, Member member) {
    const TagType t{};
    for (; member.valid_dynamic(); member.next_dynamic()) {
      functor(t, member);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    const ParallelFor &self = *((const ParallelFor *)arg);

    ParallelFor::exec_team<WorkTag, typename Policy::schedule_type::type>(
        self.m_functor, Member(&exec, self.m_policy, self.m_shared));

    exec.barrier();
    exec.fan_in();
  }
  template <typename Policy>
  Policy fix_policy(Policy policy) {
    if (policy.impl_vector_length() < 0) {
      policy.impl_set_vector_length(1);
    }
    if (policy.team_size() < 0) {
      policy.impl_set_team_size(
          policy.team_size_recommended(m_functor, ParallelForTag{}));
    }
    return policy;
  }

 public:
  inline void execute() const {
    ThreadsExec::resize_scratch(
        0, Policy::member_type::team_reduce_size() + m_shared);

    ThreadsExec::start(&ParallelFor::exec, this);

    ThreadsExec::fence();
  }

  ParallelFor(const FunctorType &arg_functor, const Policy &arg_policy)
      : m_functor(arg_functor),
        m_policy(fix_policy(arg_policy)),
        m_shared(m_policy.scratch_size(0) + m_policy.scratch_size(1) +
                 FunctorTeamShmemSize<FunctorType>::value(
                     arg_functor, m_policy.team_size())) {}
};

template <class CombinedFunctorReducerType, class... Properties>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::TeamPolicy<Properties...>, Kokkos::Threads> {
 private:
  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::Threads, Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;
  using WorkTag     = typename Policy::work_tag;
  using Member      = typename Policy::member_type;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const size_t m_shared;

  template <class TagType>
  inline static std::enable_if_t<std::is_void<TagType>::value> exec_team(
      const FunctorType &functor, Member member, reference_type update) {
    for (; member.valid_static(); member.next_static()) {
      functor(member, update);
    }
  }

  template <class TagType>
  inline static std::enable_if_t<!std::is_void<TagType>::value> exec_team(
      const FunctorType &functor, Member member, reference_type update) {
    const TagType t{};
    for (; member.valid_static(); member.next_static()) {
      functor(t, member, update);
    }
  }

  static void exec(ThreadsExec &exec, const void *arg) {
    const ParallelReduce &self = *((const ParallelReduce *)arg);

    ParallelReduce::template exec_team<WorkTag>(
        self.m_functor_reducer.get_functor(),
        Member(&exec, self.m_policy, self.m_shared),
        self.m_functor_reducer.get_reducer().init(
            static_cast<pointer_type>(exec.reduce_memory())));

    exec.fan_in_reduce(self.m_functor_reducer.get_reducer());
  }

 public:
  inline void execute() const {
    const ReducerType &reducer = m_functor_reducer.get_reducer();

    if (m_policy.league_size() * m_policy.team_size() == 0) {
      if (m_result_ptr) {
        reducer.init(m_result_ptr);
        reducer.final(m_result_ptr);
      }
    } else {
      ThreadsExec::resize_scratch(
          reducer.value_size(),
          Policy::member_type::team_reduce_size() + m_shared);

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

  template <typename Policy>
  Policy fix_policy(Policy policy) {
    if (policy.impl_vector_length() < 0) {
      policy.impl_set_vector_length(1);
    }
    if (policy.team_size() < 0) {
      policy.impl_set_team_size(policy.team_size_recommended(
          m_functor_reducer.get_functor(), m_functor_reducer.get_reducer(),
          ParallelReduceTag{}));
    }
    return policy;
  }

  template <class ViewType>
  inline ParallelReduce(const CombinedFunctorReducerType &arg_functor_reducer,
                        const Policy &arg_policy, const ViewType &arg_result)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(fix_policy(arg_policy)),
        m_result_ptr(arg_result.data()),
        m_shared(m_policy.scratch_size(0) + m_policy.scratch_size(1) +
                 FunctorTeamShmemSize<FunctorType>::value(
                     arg_functor_reducer.get_functor(), m_policy.team_size())) {
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
