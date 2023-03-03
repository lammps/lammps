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

#ifndef KOKKOS_OPENACC_PARALLEL_REDUCE_RANGE_HPP
#define KOKKOS_OPENACC_PARALLEL_REDUCE_RANGE_HPP

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_Macros.hpp>
#include <OpenACC/Kokkos_OpenACC_FunctorAdapter.hpp>
#include <OpenACC/Kokkos_OpenACC_ScheduleType.hpp>
#include <Kokkos_Parallel.hpp>
#include <type_traits>

namespace Kokkos::Experimental::Impl {

// primary template: catch-all non-implemented custom reducers
template <class Functor, class Reducer, class Policy,
          bool = std::is_arithmetic_v<typename Reducer::value_type>>
struct OpenACCParallelReduceHelper {
  OpenACCParallelReduceHelper(Functor const&, Reducer const&, Policy const&) {
    static_assert(!Kokkos::Impl::always_true<Functor>::value,
                  "not implemented");
  }
};

}  // namespace Kokkos::Experimental::Impl

template <class Functor, class ReducerType, class... Traits>
class Kokkos::Impl::ParallelReduce<Functor, Kokkos::RangePolicy<Traits...>,
                                   ReducerType, Kokkos::Experimental::OpenACC> {
  using Policy = RangePolicy<Traits...>;

  using ReducerConditional =
      if_c<std::is_same_v<InvalidType, ReducerType>, Functor, ReducerType>;
  using ReducerTypeFwd = typename ReducerConditional::type;
  using Analysis =
      FunctorAnalysis<FunctorPatternInterface::REDUCE, Policy, ReducerTypeFwd>;

  using Pointer   = typename Analysis::pointer_type;
  using ValueType = typename Analysis::value_type;

  Functor m_functor;
  Policy m_policy;
  ReducerType m_reducer;
  Pointer m_result_ptr;

 public:
  ParallelReduce(Functor const& functor, Policy const& policy,
                 ReducerType const& reducer)
      : m_functor(functor),
        m_policy(policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()) {}

  template <class ViewType>
  ParallelReduce(
      const Functor& functor, const Policy& policy, const ViewType& result,
      std::enable_if_t<Kokkos::is_view<ViewType>::value, void*> = nullptr)
      : m_functor(functor),
        m_policy(policy),
        m_reducer(InvalidType()),
        m_result_ptr(result.data()) {}

  void execute() const {
    auto const begin = m_policy.begin();
    auto const end   = m_policy.end();

    if (end <= begin) {
      return;
    }

    ValueType val;
    typename Analysis::Reducer final_reducer(
        &ReducerConditional::select(m_functor, m_reducer));
    final_reducer.init(&val);

    Kokkos::Experimental::Impl::OpenACCParallelReduceHelper(
        Kokkos::Experimental::Impl::FunctorAdapter<Functor, Policy>(m_functor),
        std::conditional_t<is_reducer_v<ReducerType>, ReducerType,
                           Sum<ValueType>>(val),
        m_policy);

    *m_result_ptr = val;
  }
};

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_SCHEDULE(REDUCER,    \
                                                              OPERATOR)   \
  namespace Kokkos::Experimental::Impl {                                  \
  template <class IndexType, class ValueType, class Functor>              \
  void OpenACCParallelReduce##REDUCER(Schedule<Static>, int chunk_size,   \
                                      IndexType begin, IndexType end,     \
                                      ValueType& aval,                    \
                                      Functor const& afunctor,            \
                                      int async_arg) {                    \
    /* FIXME_OPENACC FIXME_NVHPC workaround compiler bug (incorrect scope \
       analysis)                                                          \
       NVC++-S-1067-Cannot determine bounds for array - functor */        \
    auto const functor(afunctor);                                         \
    auto val = aval;                                                      \
    if (chunk_size >= 1) {                                                \
      /* clang-format off */ \
      KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang(static:chunk_size) vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                            \
      /* clang-format on */                                               \
      for (auto i = begin; i < end; i++) {                                \
        functor(i, val);                                                  \
      }                                                                   \
    } else {                                                              \
      /* clang-format off */ \
      KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang(static:*) vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                            \
      /* clang-format on */                                               \
      for (auto i = begin; i < end; i++) {                                \
        functor(i, val);                                                  \
      }                                                                   \
    }                                                                     \
    aval = val;                                                           \
  }                                                                       \
                                                                          \
  template <class IndexType, class ValueType, class Functor>              \
  void OpenACCParallelReduce##REDUCER(Schedule<Dynamic>, int chunk_size,  \
                                      IndexType begin, IndexType end,     \
                                      ValueType& aval,                    \
                                      Functor const& afunctor,            \
                                      int async_arg) {                    \
    /* FIXME_OPENACC FIXME_NVHPC workaround compiler bug (incorrect scope \
       analysis)                                                          \
       NVC++-S-1067-Cannot determine bounds for array - functor */        \
    auto const functor(afunctor);                                         \
    auto val = aval;                                                      \
    if (chunk_size >= 1) {                                                \
      /* clang-format off */ \
      KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang(static:chunk_size) vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                            \
      /* clang-format on */                                               \
      for (auto i = begin; i < end; i++) {                                \
        functor(i, val);                                                  \
      }                                                                   \
    } else {                                                              \
      /* clang-format off */ \
      KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                            \
      /* clang-format on */                                               \
      for (auto i = begin; i < end; i++) {                                \
        functor(i, val);                                                  \
      }                                                                   \
    }                                                                     \
    aval = val;                                                           \
  }                                                                       \
  }  // namespace Kokkos::Experimental::Impl

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(REDUCER, OPERATOR)          \
  KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_SCHEDULE(REDUCER, OPERATOR)     \
  template <class Functor, class Scalar, class Space, class... Traits>         \
  struct Kokkos::Experimental::Impl::OpenACCParallelReduceHelper<              \
      Functor, Kokkos::REDUCER<Scalar, Space>, Kokkos::RangePolicy<Traits...>, \
      true> {                                                                  \
    using Policy = RangePolicy<Traits...>;                                     \
    using ScheduleType =                                                       \
        Kokkos::Experimental::Impl::OpenACCScheduleType<Policy>;               \
    using Reducer   = REDUCER<Scalar, Space>;                                  \
    using ValueType = typename Reducer::value_type;                            \
                                                                               \
    OpenACCParallelReduceHelper(Functor const& functor,                        \
                                Reducer const& reducer,                        \
                                Policy const& policy) {                        \
      auto const begin = policy.begin();                                       \
      auto const end   = policy.end();                                         \
                                                                               \
      if (end <= begin) {                                                      \
        return;                                                                \
      }                                                                        \
                                                                               \
      ValueType val;                                                           \
      reducer.init(val);                                                       \
                                                                               \
      int const async_arg  = policy.space().acc_async_queue();                 \
      int const chunk_size = policy.chunk_size();                              \
                                                                               \
      OpenACCParallelReduce##REDUCER(ScheduleType(), chunk_size, begin, end,   \
                                     val, functor, async_arg);                 \
                                                                               \
      reducer.reference() = val;                                               \
    }                                                                          \
  }

KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(Sum, +);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(Prod, *);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(Min, min);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(Max, max);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(LAnd, &&);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(LOr, ||);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(BAnd, &);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(BOr, |);

#undef KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER
#undef KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_SCHEDULE

#endif
