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

template <class CombinedFunctorReducerType, class... Traits>
class Kokkos::Impl::ParallelReduce<CombinedFunctorReducerType,
                                   Kokkos::RangePolicy<Traits...>,
                                   Kokkos::Experimental::OpenACC> {
  using Policy      = RangePolicy<Traits...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using Pointer   = typename ReducerType::pointer_type;
  using ValueType = typename ReducerType::value_type;

  CombinedFunctorReducerType m_functor_reducer;
  Policy m_policy;
  Pointer m_result_ptr;
  bool m_result_ptr_on_device;

 public:
  template <class ViewType>
  ParallelReduce(CombinedFunctorReducerType const& functor_reducer,
                 Policy const& policy, ViewType const& result)
      : m_functor_reducer(functor_reducer),
        m_policy(policy),
        m_result_ptr(result.data()),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenACCSpace,
                              typename ViewType::memory_space>::accessible) {}

  void execute() const {
    auto const begin = m_policy.begin();
    auto const end   = m_policy.end();

    ValueType val;
    ReducerType const& reducer = m_functor_reducer.get_reducer();
    reducer.init(&val);

    if (end <= begin) {
      if (m_result_ptr_on_device == false) {
        *m_result_ptr = val;
      } else {
        acc_memcpy_to_device(m_result_ptr, &val, sizeof(ValueType));
      }
      return;
    }

    int const async_arg = m_policy.space().acc_async_queue();

    Kokkos::Experimental::Impl::OpenACCParallelReduceHelper(
        Kokkos::Experimental::Impl::FunctorAdapter<
            FunctorType, Policy,
            Kokkos::Experimental::Impl::RoutineClause::seq>(
            m_functor_reducer.get_functor()),
        std::conditional_t<
            std::is_same_v<FunctorType, typename ReducerType::functor_type>,
            Sum<ValueType>, typename ReducerType::functor_type>(val),
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
      acc_memcpy_to_device_async(m_result_ptr, &val, sizeof(ValueType),
                                 async_arg);
      acc_wait(async_arg);
    }
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
