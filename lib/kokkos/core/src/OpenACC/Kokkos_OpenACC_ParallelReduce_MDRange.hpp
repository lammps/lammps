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

#ifndef KOKKOS_OPENACC_PARALLEL_REDUCE_MDRANGE_HPP
#define KOKKOS_OPENACC_PARALLEL_REDUCE_MDRANGE_HPP

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_Macros.hpp>
#include <OpenACC/Kokkos_OpenACC_FunctorAdapter.hpp>
#include <OpenACC/Kokkos_OpenACC_MDRangePolicy.hpp>
#include <Kokkos_Parallel.hpp>

namespace Kokkos::Experimental::Impl {

// primary template: catch-all non-implemented custom reducers
template <class Functor, class Reducer, class Policy,
          bool = std::is_arithmetic_v<typename Reducer::value_type>>
struct OpenACCParallelReduceMDRangeHelper {
  OpenACCParallelReduceMDRangeHelper(Functor const&, Reducer const&,
                                     Policy const&) {
    static_assert(Kokkos::Impl::always_false<Functor>::value,
                  "not implemented");
  }
};
}  // namespace Kokkos::Experimental::Impl

template <class CombinedFunctorReducerType, class... Traits>
class Kokkos::Impl::ParallelReduce<CombinedFunctorReducerType,
                                   Kokkos::MDRangePolicy<Traits...>,
                                   Kokkos::Experimental::OpenACC> {
  using Policy      = MDRangePolicy<Traits...>;
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
  ParallelReduce(const CombinedFunctorReducerType& functor_reducer,
                 const Policy& policy, const ViewType& result)
      : m_functor_reducer(functor_reducer),
        m_policy(policy),
        m_result_ptr(result.data()),
        m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenACCSpace,
                              typename ViewType::memory_space>::accessible) {}

  void execute() const {
    static_assert(1 < Policy::rank && Policy::rank < 7);
    static_assert(Policy::inner_direction == Iterate::Left ||
                  Policy::inner_direction == Iterate::Right);
    constexpr int rank = Policy::rank;
    ValueType val;
    const ReducerType& reducer = m_functor_reducer.get_reducer();
    reducer.init(&val);

    for (int i = 0; i < rank; ++i) {
      if (m_policy.m_lower[i] >= m_policy.m_upper[i]) {
        if (m_result_ptr_on_device) {
          acc_memcpy_to_device(m_result_ptr, &val, sizeof(ValueType));
        } else {
          *m_result_ptr = val;
        }
        return;
      }
    }

    int const async_arg = m_policy.space().acc_async_queue();

    Kokkos::Experimental::Impl::OpenACCParallelReduceMDRangeHelper(
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
    if (m_result_ptr_on_device) {
      acc_memcpy_to_device_async(m_result_ptr, &val, sizeof(ValueType),
                                 async_arg);
      acc_wait(async_arg);
    } else {
      acc_wait(async_arg);
      *m_result_ptr = val;
    }
  }
};

#if defined(KOKKOS_ENABLE_OPENACC_COLLAPSE_MDRANGE_LOOPS)

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_ITERATE(REDUCER,       \
                                                             OPERATOR)      \
  namespace Kokkos::Experimental::Impl {                                    \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,  \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<2> const& begin,  \
                                      OpenACCMDRangeEnd<2> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim1 * dim0;                                              \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto i1 = m / dim0 + begin1;                                          \
      auto i0 = m % dim0 + begin0;                                          \
      functor(i0, i1, val);                                                 \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
                                                                            \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval, \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<2> const& begin,  \
                                      OpenACCMDRangeEnd<2> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim1 * dim0;                                              \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto i0 = m / dim1 + begin0;                                          \
      auto i1 = m % dim1 + begin1;                                          \
      functor(i0, i1, val);                                                 \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
                                                                            \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,  \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<3> const& begin,  \
                                      OpenACCMDRangeEnd<3> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin2 = begin[2];                                                 \
    auto end2   = end[2];                                                   \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto dim2   = end2 - begin2;                                            \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim2 * dim1 * dim0;                                       \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto tmp1 = dim1 * dim0;                                              \
      auto i2   = m / tmp1 + begin2;                                        \
      auto tmp2 = m % tmp1;                                                 \
      auto i1   = tmp2 / dim0 + begin1;                                     \
      auto i0   = tmp2 % dim0 + begin0;                                     \
      functor(i0, i1, i2, val);                                             \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
                                                                            \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval, \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<3> const& begin,  \
                                      OpenACCMDRangeEnd<3> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto begin2 = begin[2];                                                 \
    auto end2   = end[2];                                                   \
    auto dim2   = end2 - begin2;                                            \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim2 * dim1 * dim0;                                       \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto tmp1 = dim2 * dim1;                                              \
      auto i0   = m / tmp1 + begin0;                                        \
      auto tmp2 = m % tmp1;                                                 \
      auto i1   = tmp2 / dim2 + begin1;                                     \
      auto i2   = tmp2 % dim2 + begin2;                                     \
      functor(i0, i1, i2, val);                                             \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
                                                                            \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,  \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<4> const& begin,  \
                                      OpenACCMDRangeEnd<4> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin3 = begin[3];                                                 \
    auto end3   = end[3];                                                   \
    auto begin2 = begin[2];                                                 \
    auto end2   = end[2];                                                   \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto dim3   = end3 - begin3;                                            \
    auto dim2   = end2 - begin2;                                            \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim3 * dim2 * dim1 * dim0;                                \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto tmp1 = dim2 * dim1 * dim0;                                       \
      auto i3   = m / tmp1 + begin3;                                        \
      auto tmp2 = m % tmp1;                                                 \
      tmp1      = dim1 * dim0;                                              \
      auto i2   = tmp2 / tmp1 + begin2;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      auto i1   = tmp2 / dim0 + begin1;                                     \
      auto i0   = tmp2 % dim0 + begin0;                                     \
      functor(i0, i1, i2, i3, val);                                         \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
                                                                            \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval, \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<4> const& begin,  \
                                      OpenACCMDRangeEnd<4> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto begin2 = begin[2];                                                 \
    auto end2   = end[2];                                                   \
    auto begin3 = begin[3];                                                 \
    auto end3   = end[3];                                                   \
    auto dim3   = end3 - begin3;                                            \
    auto dim2   = end2 - begin2;                                            \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim3 * dim2 * dim1 * dim0;                                \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto tmp1 = dim3 * dim2 * dim1;                                       \
      auto i0   = m / tmp1 + begin0;                                        \
      auto tmp2 = m % tmp1;                                                 \
      tmp1      = dim3 * dim2;                                              \
      auto i1   = tmp2 / tmp1 + begin1;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      auto i2   = tmp2 / dim3 + begin2;                                     \
      auto i3   = tmp2 % dim3 + begin3;                                     \
      functor(i0, i1, i2, i3, val);                                         \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
                                                                            \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,  \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<5> const& begin,  \
                                      OpenACCMDRangeEnd<5> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin4 = begin[4];                                                 \
    auto end4   = end[4];                                                   \
    auto begin3 = begin[3];                                                 \
    auto end3   = end[3];                                                   \
    auto begin2 = begin[2];                                                 \
    auto end2   = end[2];                                                   \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto dim4   = end4 - begin4;                                            \
    auto dim3   = end3 - begin3;                                            \
    auto dim2   = end2 - begin2;                                            \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim4 * dim3 * dim2 * dim1 * dim0;                         \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto tmp1 = dim3 * dim2 * dim1 * dim0;                                \
      auto i4   = m / tmp1 + begin4;                                        \
      auto tmp2 = m % tmp1;                                                 \
      tmp1      = dim2 * dim1 * dim0;                                       \
      auto i3   = tmp2 / tmp1 + begin3;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      tmp1      = dim1 * dim0;                                              \
      auto i2   = tmp2 / tmp1 + begin2;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      auto i1   = tmp2 / dim0 + begin1;                                     \
      auto i0   = tmp2 % dim0 + begin0;                                     \
      functor(i0, i1, i2, i3, i4, val);                                     \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
                                                                            \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval, \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<5> const& begin,  \
                                      OpenACCMDRangeEnd<5> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto begin2 = begin[2];                                                 \
    auto end2   = end[2];                                                   \
    auto begin3 = begin[3];                                                 \
    auto end3   = end[3];                                                   \
    auto begin4 = begin[4];                                                 \
    auto end4   = end[4];                                                   \
    auto dim4   = end4 - begin4;                                            \
    auto dim3   = end3 - begin3;                                            \
    auto dim2   = end2 - begin2;                                            \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim4 * dim3 * dim2 * dim1 * dim0;                         \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto tmp1 = dim4 * dim3 * dim2 * dim1;                                \
      auto i0   = m / tmp1 + begin0;                                        \
      auto tmp2 = m % tmp1;                                                 \
      tmp1      = dim4 * dim3 * dim2;                                       \
      auto i1   = tmp2 / tmp1 + begin1;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      tmp1      = dim4 * dim3;                                              \
      auto i2   = tmp2 / tmp1 + begin2;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      auto i3   = tmp2 / dim4 + begin3;                                     \
      auto i4   = tmp2 % dim4 + begin4;                                     \
      functor(i0, i1, i2, i3, i4, val);                                     \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
                                                                            \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,  \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<6> const& begin,  \
                                      OpenACCMDRangeEnd<6> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin5 = begin[5];                                                 \
    auto end5   = end[5];                                                   \
    auto begin4 = begin[4];                                                 \
    auto end4   = end[4];                                                   \
    auto begin3 = begin[3];                                                 \
    auto end3   = end[3];                                                   \
    auto begin2 = begin[2];                                                 \
    auto end2   = end[2];                                                   \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto dim5   = end5 - begin5;                                            \
    auto dim4   = end4 - begin4;                                            \
    auto dim3   = end3 - begin3;                                            \
    auto dim2   = end2 - begin2;                                            \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim5 * dim4 * dim3 * dim2 * dim1 * dim0;                  \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto tmp1 = dim4 * dim3 * dim2 * dim1 * dim0;                         \
      auto i5   = m / tmp1 + begin5;                                        \
      auto tmp2 = m % tmp1;                                                 \
      tmp1      = dim3 * dim2 * dim1 * dim0;                                \
      auto i4   = tmp2 / tmp1 + begin4;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      tmp1      = dim2 * dim1 * dim0;                                       \
      auto i3   = tmp2 / tmp1 + begin3;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      tmp1      = dim1 * dim0;                                              \
      auto i2   = tmp2 / tmp1 + begin2;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      auto i1   = tmp2 / dim0 + begin1;                                     \
      auto i0   = tmp2 % dim0 + begin0;                                     \
      functor(i0, i1, i2, i3, i4, i5, val);                                 \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
                                                                            \
  template <class ValueType, class Functor>                                 \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval, \
                                      Functor const& afunctor,              \
                                      OpenACCMDRangeBegin<6> const& begin,  \
                                      OpenACCMDRangeEnd<6> const& end,      \
                                      int async_arg) {                      \
    auto val = aval;                                                        \
    auto const functor(afunctor);                                           \
    auto begin0 = begin[0];                                                 \
    auto end0   = end[0];                                                   \
    auto begin1 = begin[1];                                                 \
    auto end1   = end[1];                                                   \
    auto begin2 = begin[2];                                                 \
    auto end2   = end[2];                                                   \
    auto begin3 = begin[3];                                                 \
    auto end3   = end[3];                                                   \
    auto begin4 = begin[4];                                                 \
    auto end4   = end[4];                                                   \
    auto begin5 = begin[5];                                                 \
    auto end5   = end[5];                                                   \
    auto dim5   = end5 - begin5;                                            \
    auto dim4   = end4 - begin4;                                            \
    auto dim3   = end3 - begin3;                                            \
    auto dim2   = end2 - begin2;                                            \
    auto dim1   = end1 - begin1;                                            \
    auto dim0   = end0 - begin0;                                            \
    auto nIter  = dim5 * dim4 * dim3 * dim2 * dim1 * dim0;                  \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                \
    /* clang-format on */                                                   \
    for (decltype(nIter) m = 0; m < nIter; ++m) {                           \
      auto tmp1 = dim5 * dim4 * dim3 * dim2 * dim1;                         \
      auto i0   = m / tmp1 + begin0;                                        \
      auto tmp2 = m % tmp1;                                                 \
      tmp1      = dim5 * dim4 * dim3 * dim2;                                \
      auto i1   = tmp2 / tmp1 + begin1;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      tmp1      = dim5 * dim4 * dim3;                                       \
      auto i2   = tmp2 / tmp1 + begin2;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      tmp1      = dim5 * dim4;                                              \
      auto i3   = tmp2 / tmp1 + begin3;                                     \
      tmp2      = tmp2 % tmp1;                                              \
      auto i4   = tmp2 / dim5 + begin4;                                     \
      auto i5   = tmp2 % dim5 + begin5;                                     \
      functor(i0, i1, i2, i3, i4, i5, val);                                 \
    }                                                                       \
    acc_wait(async_arg);                                                    \
    aval = val;                                                             \
  }                                                                         \
  }  // namespace Kokkos::Experimental::Impl

#else

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_ITERATE(REDUCER,         \
                                                             OPERATOR)        \
  namespace Kokkos::Experimental::Impl {                                      \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,    \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<2> const& begin,    \
                                      OpenACCMDRangeEnd<2> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(2) reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                  \
    /* clang-format on */                                                     \
    for (auto i1 = begin1; i1 < end1; ++i1) {                                 \
      for (auto i0 = begin0; i0 < end0; ++i0) {                               \
        functor(i0, i1, val);                                                 \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
                                                                              \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval,   \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<2> const& begin,    \
                                      OpenACCMDRangeEnd<2> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(2) reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                  \
    /* clang-format on */                                                     \
    for (auto i0 = begin0; i0 < end0; ++i0) {                                 \
      for (auto i1 = begin1; i1 < end1; ++i1) {                               \
        functor(i0, i1, val);                                                 \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
                                                                              \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,    \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<3> const& begin,    \
                                      OpenACCMDRangeEnd<3> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin2 = begin[2];                                                   \
    auto end2   = end[2];                                                     \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    /* clang-format off */                                                  \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(3) reduction( \
        OPERATOR                                                            \
        : val) copyin(functor) async(async_arg)) \
    /* clang-format on */                                                     \
    for (auto i2 = begin2; i2 < end2; ++i2) {                                 \
      for (auto i1 = begin1; i1 < end1; ++i1) {                               \
        for (auto i0 = begin0; i0 < end0; ++i0) {                             \
          functor(i0, i1, i2, val);                                           \
        }                                                                     \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
                                                                              \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval,   \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<3> const& begin,    \
                                      OpenACCMDRangeEnd<3> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    auto begin2 = begin[2];                                                   \
    auto end2   = end[2];                                                     \
    /* clang-format off */                                                  \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(3) reduction( \
        OPERATOR                                                            \
        : val) copyin(functor) async(async_arg)) \
    /* clang-format on */                                                     \
    for (auto i0 = begin0; i0 < end0; ++i0) {                                 \
      for (auto i1 = begin1; i1 < end1; ++i1) {                               \
        for (auto i2 = begin2; i2 < end2; ++i2) {                             \
          functor(i0, i1, i2, val);                                           \
        }                                                                     \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
                                                                              \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,    \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<4> const& begin,    \
                                      OpenACCMDRangeEnd<4> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin3 = begin[3];                                                   \
    auto end3   = end[3];                                                     \
    auto begin2 = begin[2];                                                   \
    auto end2   = end[2];                                                     \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(4) reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                  \
    /* clang-format on */                                                     \
    for (auto i3 = begin3; i3 < end3; ++i3) {                                 \
      for (auto i2 = begin2; i2 < end2; ++i2) {                               \
        for (auto i1 = begin1; i1 < end1; ++i1) {                             \
          for (auto i0 = begin0; i0 < end0; ++i0) {                           \
            functor(i0, i1, i2, i3, val);                                     \
          }                                                                   \
        }                                                                     \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
                                                                              \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval,   \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<4> const& begin,    \
                                      OpenACCMDRangeEnd<4> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    auto begin2 = begin[2];                                                   \
    auto end2   = end[2];                                                     \
    auto begin3 = begin[3];                                                   \
    auto end3   = end[3];                                                     \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(4) reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                  \
    /* clang-format on */                                                     \
    for (auto i0 = begin0; i0 < end0; ++i0) {                                 \
      for (auto i1 = begin1; i1 < end1; ++i1) {                               \
        for (auto i2 = begin2; i2 < end2; ++i2) {                             \
          for (auto i3 = begin3; i3 < end3; ++i3) {                           \
            functor(i0, i1, i2, i3, val);                                     \
          }                                                                   \
        }                                                                     \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
                                                                              \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,    \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<5> const& begin,    \
                                      OpenACCMDRangeEnd<5> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin4 = begin[4];                                                   \
    auto end4   = end[4];                                                     \
    auto begin3 = begin[3];                                                   \
    auto end3   = end[3];                                                     \
    auto begin2 = begin[2];                                                   \
    auto end2   = end[2];                                                     \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(5) reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                  \
    /* clang-format on */                                                     \
    for (auto i4 = begin4; i4 < end4; ++i4) {                                 \
      for (auto i3 = begin3; i3 < end3; ++i3) {                               \
        for (auto i2 = begin2; i2 < end2; ++i2) {                             \
          for (auto i1 = begin1; i1 < end1; ++i1) {                           \
            for (auto i0 = begin0; i0 < end0; ++i0) {                         \
              functor(i0, i1, i2, i3, i4, val);                               \
            }                                                                 \
          }                                                                   \
        }                                                                     \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
                                                                              \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval,   \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<5> const& begin,    \
                                      OpenACCMDRangeEnd<5> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    auto begin2 = begin[2];                                                   \
    auto end2   = end[2];                                                     \
    auto begin3 = begin[3];                                                   \
    auto end3   = end[3];                                                     \
    auto begin4 = begin[4];                                                   \
    auto end4   = end[4];                                                     \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(5) reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                  \
    /* clang-format on */                                                     \
    for (auto i0 = begin0; i0 < end0; ++i0) {                                 \
      for (auto i1 = begin1; i1 < end1; ++i1) {                               \
        for (auto i2 = begin2; i2 < end2; ++i2) {                             \
          for (auto i3 = begin3; i3 < end3; ++i3) {                           \
            for (auto i4 = begin4; i4 < end4; ++i4) {                         \
              functor(i0, i1, i2, i3, i4, val);                               \
            }                                                                 \
          }                                                                   \
        }                                                                     \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
                                                                              \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateLeft, ValueType& aval,    \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<6> const& begin,    \
                                      OpenACCMDRangeEnd<6> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin5 = begin[5];                                                   \
    auto end5   = end[5];                                                     \
    auto begin4 = begin[4];                                                   \
    auto end4   = end[4];                                                     \
    auto begin3 = begin[3];                                                   \
    auto end3   = end[3];                                                     \
    auto begin2 = begin[2];                                                   \
    auto end2   = end[2];                                                     \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(6) reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                  \
    /* clang-format on */                                                     \
    for (auto i5 = begin5; i5 < end5; ++i5) {                                 \
      for (auto i4 = begin4; i4 < end4; ++i4) {                               \
        for (auto i3 = begin3; i3 < end3; ++i3) {                             \
          for (auto i2 = begin2; i2 < end2; ++i2) {                           \
            for (auto i1 = begin1; i1 < end1; ++i1) {                         \
              for (auto i0 = begin0; i0 < end0; ++i0) {                       \
                functor(i0, i1, i2, i3, i4, i5, val);                         \
              }                                                               \
            }                                                                 \
          }                                                                   \
        }                                                                     \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
                                                                              \
  template <class ValueType, class Functor>                                   \
  void OpenACCParallelReduce##REDUCER(OpenACCIterateRight, ValueType& aval,   \
                                      Functor const& afunctor,                \
                                      OpenACCMDRangeBegin<6> const& begin,    \
                                      OpenACCMDRangeEnd<6> const& end,        \
                                      int async_arg) {                        \
    auto val = aval;                                                          \
    auto const functor(afunctor);                                             \
    auto begin0 = begin[0];                                                   \
    auto end0   = end[0];                                                     \
    auto begin1 = begin[1];                                                   \
    auto end1   = end[1];                                                     \
    auto begin2 = begin[2];                                                   \
    auto end2   = end[2];                                                     \
    auto begin3 = begin[3];                                                   \
    auto end3   = end[3];                                                     \
    auto begin4 = begin[4];                                                   \
    auto end4   = end[4];                                                     \
    auto begin5 = begin[5];                                                   \
    auto end5   = end[5];                                                     \
    /* clang-format off */ \
    KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector collapse(6) reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                  \
    /* clang-format on */                                                     \
    for (auto i0 = begin0; i0 < end0; ++i0) {                                 \
      for (auto i1 = begin1; i1 < end1; ++i1) {                               \
        for (auto i2 = begin2; i2 < end2; ++i2) {                             \
          for (auto i3 = begin3; i3 < end3; ++i3) {                           \
            for (auto i4 = begin4; i4 < end4; ++i4) {                         \
              for (auto i5 = begin5; i5 < end5; ++i5) {                       \
                functor(i0, i1, i2, i3, i4, i5, val);                         \
              }                                                               \
            }                                                                 \
          }                                                                   \
        }                                                                     \
      }                                                                       \
    }                                                                         \
    acc_wait(async_arg);                                                      \
    aval = val;                                                               \
  }                                                                           \
  }  // namespace Kokkos::Experimental::Impl

#endif

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER(REDUCER, OPERATOR) \
  KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_ITERATE(REDUCER, OPERATOR)     \
  template <class Functor, class Scalar, class Space, class... Traits>        \
  struct Kokkos::Experimental::Impl::OpenACCParallelReduceMDRangeHelper<      \
      Functor, Kokkos::REDUCER<Scalar, Space>,                                \
      Kokkos::MDRangePolicy<Traits...>, true> {                               \
    using Policy    = MDRangePolicy<Traits...>;                               \
    using Reducer   = REDUCER<Scalar, Space>;                                 \
    using ValueType = typename Reducer::value_type;                           \
                                                                              \
    OpenACCParallelReduceMDRangeHelper(Functor const& functor,                \
                                       Reducer const& reducer,                \
                                       Policy const& policy) {                \
      ValueType val;                                                          \
      reducer.init(val);                                                      \
                                                                              \
      int const async_arg = policy.space().acc_async_queue();                 \
                                                                              \
      OpenACCParallelReduce##REDUCER(                                         \
          std::integral_constant<Iterate, Policy::inner_direction>(), val,    \
          functor, policy.m_lower, policy.m_upper, async_arg);                \
                                                                              \
      reducer.reference() = val;                                              \
    }                                                                         \
  }

KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER(Sum, +);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER(Prod, *);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER(Min, min);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER(Max, max);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER(LAnd, &&);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER(LOr, ||);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER(BAnd, &);
KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER(BOr, |);

#undef KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_MDRANGE_HELPER
#undef KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_ITERATE

#endif
