/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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

template <class FunctorType, class ReducerType, class... Properties>
class ParallelReduce<FunctorType, Kokkos::TeamPolicy<Properties...>,
                     ReducerType, Kokkos::Threads> {
 private:
  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::Threads, Properties...>;
  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  using ReducerConditional =
      Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                         FunctorType, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using WorkTagFwd =
      typename Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                                  WorkTag, void>::type;

  using Analysis = Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                                         Policy, ReducerTypeFwd>;
  using pointer_type   = typename Analysis::pointer_type;
  using reference_type = typename Analysis::reference_type;

  const FunctorType m_functor;
  const Policy m_policy;
  const ReducerType m_reducer;
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

    typename Analysis::Reducer reducer(
        &ReducerConditional::select(self.m_functor, self.m_reducer));

    ParallelReduce::template exec_team<WorkTag>(
        self.m_functor, Member(&exec, self.m_policy, self.m_shared),
        reducer.init(static_cast<pointer_type>(exec.reduce_memory())));

    exec.fan_in_reduce(reducer);
  }

 public:
  inline void execute() const {
    if (m_policy.league_size() * m_policy.team_size() == 0) {
      if (m_result_ptr) {
        typename Analysis::Reducer final_reducer(
            &ReducerConditional::select(m_functor, m_reducer));
        final_reducer.init(m_result_ptr);
        final_reducer.final(m_result_ptr);
      }
    } else {
      ThreadsExec::resize_scratch(
          Analysis::value_size(
              ReducerConditional::select(m_functor, m_reducer)),
          Policy::member_type::team_reduce_size() + m_shared);

      ThreadsExec::start(&ParallelReduce::exec, this);

      ThreadsExec::fence();

      if (m_result_ptr) {
        const pointer_type data =
            (pointer_type)ThreadsExec::root_reduce_scratch();

        const unsigned n = Analysis::value_count(
            ReducerConditional::select(m_functor, m_reducer));
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
          m_functor, m_reducer, ParallelReduceTag{}));
    }
    return policy;
  }

  template <class ViewType>
  inline ParallelReduce(
      const FunctorType &arg_functor, const Policy &arg_policy,
      const ViewType &arg_result,
      std::enable_if_t<Kokkos::is_view<ViewType>::value &&
                           !Kokkos::is_reducer<ReducerType>::value,
                       void *> = nullptr)
      : m_functor(arg_functor),
        m_policy(fix_policy(arg_policy)),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result.data()),
        m_shared(m_policy.scratch_size(0) + m_policy.scratch_size(1) +
                 FunctorTeamShmemSize<FunctorType>::value(
                     arg_functor, m_policy.team_size())) {}

  inline ParallelReduce(const FunctorType &arg_functor, Policy arg_policy,
                        const ReducerType &reducer)
      : m_functor(arg_functor),
        m_policy(fix_policy(arg_policy)),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()),
        m_shared(m_policy.scratch_size(0) + m_policy.scratch_size(1) +
                 FunctorTeamShmemSize<FunctorType>::value(
                     arg_functor, m_policy.team_size())) {
    /*static_assert( std::is_same< typename ViewType::memory_space
                            , Kokkos::HostSpace >::value
    , "Reduction result on Kokkos::OpenMP must be a Kokkos::View in HostSpace"
    );*/
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
