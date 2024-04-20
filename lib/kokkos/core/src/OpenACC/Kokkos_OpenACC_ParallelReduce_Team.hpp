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

#ifndef KOKKOS_OPENACC_PARALLEL_REDUCE_TEAM_HPP
#define KOKKOS_OPENACC_PARALLEL_REDUCE_TEAM_HPP

#include <OpenACC/Kokkos_OpenACC_Team.hpp>
#include <OpenACC/Kokkos_OpenACC_FunctorAdapter.hpp>
#include <OpenACC/Kokkos_OpenACC_Macros.hpp>

#ifdef KOKKOS_ENABLE_OPENACC_COLLAPSE_HIERARCHICAL_CONSTRUCTS
#define KOKKOS_IMPL_OPENACC_LOOP_CLAUSE \
  Kokkos::Experimental::Impl::RoutineClause::seq
#else
#define KOKKOS_IMPL_OPENACC_LOOP_CLAUSE \
  Kokkos::Experimental::Impl::RoutineClause::worker
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Hierarchical Parallelism -> Team level implementation
namespace Kokkos::Experimental::Impl {

// primary template: catch-all non-implemented custom reducers
template <class Functor, class Reducer, class Policy,
          bool = std::is_arithmetic_v<typename Reducer::value_type>>
struct OpenACCParallelReduceTeamHelper {
  OpenACCParallelReduceTeamHelper(Functor const&, Reducer const&,
                                  Policy const&) {
    static_assert(Kokkos::Impl::always_false<Functor>::value,
                  "not implemented");
  }
};

}  // namespace Kokkos::Experimental::Impl

template <class CombinedFunctorReducerType, class... Properties>
class Kokkos::Impl::ParallelReduce<CombinedFunctorReducerType,
                                   Kokkos::TeamPolicy<Properties...>,
                                   Kokkos::Experimental::OpenACC> {
 private:
  using Policy =
      TeamPolicyInternal<Kokkos::Experimental::OpenACC, Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using value_type   = typename ReducerType::value_type;
  using pointer_type = typename ReducerType::pointer_type;

  CombinedFunctorReducerType m_functor_reducer;
  Policy m_policy;
  pointer_type m_result_ptr;
  bool m_result_ptr_on_device;

 public:
  void execute() const {
    auto league_size   = m_policy.league_size();
    auto team_size     = m_policy.team_size();
    auto vector_length = m_policy.impl_vector_length();

    int const async_arg = m_policy.space().acc_async_queue();
    value_type val;
    const ReducerType& reducer = m_functor_reducer.get_reducer();
    reducer.init(&val);
    if (league_size <= 0) {
      if (m_result_ptr_on_device == false) {
        *m_result_ptr = val;
      } else {
        acc_memcpy_to_device(m_result_ptr, &val, sizeof(value_type));
      }
      return;
    }

    Kokkos::Experimental::Impl::OpenACCParallelReduceTeamHelper(
        Kokkos::Experimental::Impl::FunctorAdapter<
            FunctorType, Policy, KOKKOS_IMPL_OPENACC_LOOP_CLAUSE>(
            m_functor_reducer.get_functor()),
        std::conditional_t<
            std::is_same_v<FunctorType, typename ReducerType::functor_type>,
            Sum<value_type>, typename ReducerType::functor_type>(val),
        m_policy);

    // OpenACC backend supports only built-in Reducer types; thus
    // reducer.final() below is a no-op.
    reducer.final(&val);
    // acc_wait(async_arg) in the below if-else statements is needed because the
    // above OpenACC compute kernel can be executed asynchronously and val is a
    // local host variable.
    if (m_result_ptr_on_device == false) {
      acc_wait(async_arg);
      *m_result_ptr = val;
    } else {
      acc_memcpy_to_device_async(m_result_ptr, &val, sizeof(value_type),
                                 async_arg);
      acc_wait(async_arg);
    }
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const Policy& arg_policy, const ViewType& arg_result_view)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result_view.data()),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenACCSpace,
                              typename ViewType::memory_space>::accessible) {}
};

namespace Kokkos {

// Hierarchical Parallelism -> Team thread level implementation
// FIXME_OPENACC: custom reduction is not implemented.
template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::OpenACCTeamMember>&
        loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& init_result) {
  static_assert(Kokkos::Impl::always_false<Lambda>::value,
                "custom reduction is not implemented");
}

// Hierarchical Parallelism -> Thread vector level implementation
// FIXME_OPENACC: custom reduction is not implemented.
template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::OpenACCTeamMember>& loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& init_result) {
  static_assert(Kokkos::Impl::always_false<Lambda>::value,
                "custom reduction is not implemented");
}

}  // namespace Kokkos

#ifdef KOKKOS_ENABLE_OPENACC_COLLAPSE_HIERARCHICAL_CONSTRUCTS

#define KOKKOS_IMPL_ACC_REDUCE_TEAM_PRAGMA \
  vector vector_length(team_size* vector_length)
#define KOKKOS_IMPL_ACC_REDUCE_TEAM_ITRS league_size* team_size* vector_length
#define KOKKOS_IMPL_ACC_REDUCE_TEAM_LEAGUE_ID_INIT \
  i / (team_size * vector_length)

namespace Kokkos {

// Hierarchical Parallelism -> Team thread level implementation
#pragma acc routine seq
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer_v<ValueType>>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::OpenACCTeamMember>& loop_boundaries,
                const Lambda& lambda, ValueType& result) {
  ValueType tmp = ValueType();
  iType j_start =
      loop_boundaries.team.team_rank() / loop_boundaries.team.vector_length();
  if (j_start == 0) {
#pragma acc loop seq
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++)
      lambda(i, tmp);
    result = tmp;
  }
}

#pragma acc routine seq
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer_v<ReducerType>>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::OpenACCTeamMember>& loop_boundaries,
                const Lambda& lambda, const ReducerType& reducer) {
  using ValueType = typename ReducerType::value_type;
  ValueType tmp;
  reducer.init(tmp);
  iType j_start =
      loop_boundaries.team.team_rank() / loop_boundaries.team.vector_length();
  if (j_start == 0) {
#pragma acc loop seq
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++)
      lambda(i, tmp);
    reducer.reference() = tmp;
  }
}

// Hierarchical Parallelism -> Thread vector level implementation
#pragma acc routine seq
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer_v<ValueType>>
parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<
                    iType, Impl::OpenACCTeamMember>& loop_boundaries,
                const Lambda& lambda, ValueType& result) {
  ValueType tmp = ValueType();
  iType j_start =
      loop_boundaries.team.team_rank() % loop_boundaries.team.vector_length();
  if (j_start == 0) {
#pragma acc loop seq
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      lambda(i, tmp);
    }
    result = tmp;
  }
}

#pragma acc routine seq
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer_v<ReducerType>>
parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<
                    iType, Impl::OpenACCTeamMember>& loop_boundaries,
                const Lambda& lambda, const ReducerType& reducer) {
  using ValueType = typename ReducerType::value_type;
  ValueType tmp;
  reducer.init(tmp);
  iType j_start =
      loop_boundaries.team.team_rank() % loop_boundaries.team.vector_length();
  if (j_start == 0) {
#pragma acc loop seq
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      lambda(i, tmp);
    }
    reducer.reference() = tmp;
  }
}

// Hierarchical Parallelism -> Team vector level implementation
#pragma acc routine seq
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamVectorRangeBoundariesStruct<iType, Impl::OpenACCTeamMember>&
        loop_boundaries,
    const Lambda& lambda, ValueType& result) {
  ValueType tmp = ValueType();
  iType j_start =
      loop_boundaries.team.team_rank() % loop_boundaries.team.vector_length();
  if (j_start == 0) {
#pragma acc loop seq
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      lambda(i, tmp);
    }
    result = tmp;
  }
}

}  // namespace Kokkos

#else /* #ifdef KOKKOS_ENABLE_OPENACC_COLLAPSE_HIERARCHICAL_CONSTRUCTS */

#define KOKKOS_IMPL_ACC_REDUCE_TEAM_PRAGMA \
  num_workers(team_size) vector_length(vector_length)
#define KOKKOS_IMPL_ACC_REDUCE_TEAM_ITRS league_size
#define KOKKOS_IMPL_ACC_REDUCE_TEAM_LEAGUE_ID_INIT i

// FIXME_OPENACC: below implementation conforms to the OpenACC standard, but
// the NVHPC compiler (V22.11) fails due to the lack of support for lambda
// expressions containing parallel loops.

namespace Kokkos {

// Hierarchical Parallelism -> Team thread level implementation
#pragma acc routine worker
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer_v<ValueType>>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::OpenACCTeamMember>& loop_boundaries,
                const Lambda& lambda, ValueType& result) {
  ValueType tmp = ValueType();
#pragma acc loop worker reduction(+ : tmp)
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++)
    lambda(i, tmp);
  result = tmp;
}

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD(REDUCER, OPERATOR)   \
  KOKKOS_IMPL_ACC_PRAGMA(routine worker)                                     \
  template <typename iType, class Lambda, class Scalar, class Space>         \
  KOKKOS_INLINE_FUNCTION                                                     \
      std::enable_if_t<Kokkos::is_reducer_v<Kokkos::REDUCER<Scalar, Space>>> \
      parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<           \
                          iType, Impl::OpenACCTeamMember>& loop_boundaries,  \
                      const Lambda& lambda,                                  \
                      const Kokkos::REDUCER<Scalar, Space>& reducer) {       \
    using ValueType = typename Kokkos::REDUCER<Scalar, Space>::value_type;   \
    ValueType tmp   = ValueType();                                           \
    reducer.init(tmp);                                                       \
    KOKKOS_IMPL_ACC_PRAGMA(loop worker reduction(OPERATOR : tmp))            \
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++)      \
      lambda(i, tmp);                                                        \
    reducer.reference() = tmp;                                               \
  }

KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD(Sum, +);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD(Prod, *);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD(Min, min);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD(Max, max);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD(LAnd, &&);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD(LOr, ||);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD(BAnd, &);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD(BOr, |);

// Hierarchical Parallelism -> Thread vector level implementation
#pragma acc routine vector
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer_v<ValueType>>
parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<
                    iType, Impl::OpenACCTeamMember>& loop_boundaries,
                const Lambda& lambda, ValueType& result) {
  ValueType tmp = ValueType();
#pragma acc loop vector reduction(+ : tmp)
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, tmp);
  }
  result = tmp;
}

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR(REDUCER, OPERATOR) \
  KOKKOS_IMPL_ACC_PRAGMA(routine vector)                                     \
  template <typename iType, class Lambda, class Scalar, class Space>         \
  KOKKOS_INLINE_FUNCTION                                                     \
      std::enable_if_t<Kokkos::is_reducer_v<Kokkos::REDUCER<Scalar, Space>>> \
      parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<         \
                          iType, Impl::OpenACCTeamMember>& loop_boundaries,  \
                      const Lambda& lambda,                                  \
                      const Kokkos::REDUCER<Scalar, Space>& reducer) {       \
    using ValueType = typename Kokkos::REDUCER<Scalar, Space>::value_type;   \
    ValueType tmp;                                                           \
    reducer.init(tmp);                                                       \
    KOKKOS_IMPL_ACC_PRAGMA(loop vector reduction(OPERATOR : tmp))            \
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {    \
      lambda(i, tmp);                                                        \
    }                                                                        \
    reducer.reference() = tmp;                                               \
  }

KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR(Sum, +);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR(Prod, *);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR(Min, min);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR(Max, max);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR(LAnd, &&);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR(LOr, ||);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR(BAnd, &);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR(BOr, |);

// Hierarchical Parallelism -> Team vector level implementation
#pragma acc routine vector
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamVectorRangeBoundariesStruct<iType, Impl::OpenACCTeamMember>&
        loop_boundaries,
    const Lambda& lambda, ValueType& result) {
  ValueType tmp = ValueType();
#pragma acc loop vector reduction(+ : tmp)
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, tmp);
  }
  result = tmp;
}

}  // namespace Kokkos

#undef KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_THREAD
#undef KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_THREAD_VECTOR

#endif /* #ifdef KOKKOS_ENABLE_OPENACC_COLLAPSE_HIERARCHICAL_CONSTRUCTS */

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_SCHEDULE(REDUCER,     \
                                                              OPERATOR)    \
  namespace Kokkos::Experimental::Impl {                                   \
  template <class Policy, class ValueType, class Functor>                  \
  void OpenACCParallelReduceTeam##REDUCER(Policy const policy,             \
                                          ValueType& aval,                 \
                                          Functor const& afunctor,         \
                                          int async_arg) {                 \
    auto const functor       = afunctor;                                   \
    auto val                 = aval;                                       \
    auto const league_size   = policy.league_size();                       \
    auto const team_size     = policy.team_size();                         \
    auto const vector_length = policy.impl_vector_length();                \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang num_gangs(league_size) KOKKOS_IMPL_ACC_REDUCE_TEAM_PRAGMA reduction(OPERATOR : val) copyin(functor) async(async_arg))                                               \
    /* clang-format on */                                                  \
    for (int i = 0; i < KOKKOS_IMPL_ACC_REDUCE_TEAM_ITRS; i++) {           \
      int league_id = KOKKOS_IMPL_ACC_REDUCE_TEAM_LEAGUE_ID_INIT;          \
      typename Policy::member_type team(league_id, league_size, team_size, \
                                        vector_length);                    \
      functor(team, val);                                                  \
    }                                                                      \
    acc_wait(async_arg);                                                   \
    aval = val;                                                            \
  }                                                                        \
  }  // namespace Kokkos::Experimental::Impl

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER(REDUCER, OPERATOR) \
  KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_SCHEDULE(REDUCER, OPERATOR) \
                                                                           \
  template <class Functor, class Scalar, class Space, class... Traits>     \
  struct Kokkos::Experimental::Impl::OpenACCParallelReduceTeamHelper<      \
      Functor, Kokkos::REDUCER<Scalar, Space>,                             \
      Kokkos::Impl::TeamPolicyInternal<Traits...>, true> {                 \
    using Policy    = Kokkos::Impl::TeamPolicyInternal<Traits...>;         \
    using Reducer   = REDUCER<Scalar, Space>;                              \
    using ValueType = typename Reducer::value_type;                        \
                                                                           \
    OpenACCParallelReduceTeamHelper(Functor const& functor,                \
                                    Reducer const& reducer,                \
                                    Policy const& policy) {                \
      auto league_size   = policy.league_size();                           \
      auto team_size     = policy.team_size();                             \
      auto vector_length = policy.impl_vector_length();                    \
                                                                           \
      if (league_size <= 0) {                                              \
        return;                                                            \
      }                                                                    \
                                                                           \
      ValueType val;                                                       \
      reducer.init(val);                                                   \
                                                                           \
      int const async_arg = policy.space().acc_async_queue();              \
                                                                           \
      OpenACCParallelReduceTeam##REDUCER(policy, val, functor, async_arg); \
                                                                           \
      reducer.reference() = val;                                           \
    }                                                                      \
  }

KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER(Sum, +);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER(Prod, *);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER(Min, min);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER(Max, max);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER(LAnd, &&);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER(LOr, ||);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER(BAnd, &);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER(BOr, |);

#undef KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_TEAM_HELPER
#undef KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_SCHEDULE
#undef KOKKOS_IMPL_ACC_REDUCE_TEAM_PRAGMA
#undef KOKKOS_IMPL_ACC_REDUCE_TEAM_ITRS
#undef KOKKOS_IMPL_ACC_REDUCE_TEAM_LEAGUE_ID_INIT

#endif /* #ifndef KOKKOS_OPENACC_PARALLEL_REDUCE_TEAM_HPP */
