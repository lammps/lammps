// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENACC_PARALLEL_REDUCE_RANGE_HPP
#define KOKKOS_OPENACC_PARALLEL_REDUCE_RANGE_HPP

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACC_Macros.hpp>
#include <OpenACC/Kokkos_OpenACC_FunctorAdapter.hpp>
#include <OpenACC/Kokkos_OpenACC_ScheduleType.hpp>
#include <Kokkos_Parallel.hpp>
#include <type_traits>

// Clacc uses an alternative implementation to work around not-yet-implemented
// OpenACC features: Clacc does not fully support private clauses for
// gang-private variables, and the alternative implementation allocates
// the gang-private arrays on GPU global memory using array expansion,
// instead of using the private clause.
#ifdef KOKKOS_COMPILER_CLANG
#define KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS1(THREADID) \
  vector_red_temp[team_id * chunk_size + THREADID]
#define KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(THREADID) \
  vector_red_temp[THREADID]
#define KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_CLAUSE1 \
  create(vector_red_temp [0:n_chunks * chunk_size])
#define KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_CLAUSE2 \
  create(vector_red_temp [0:chunk_size])
#else
#define KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS1(THREADID) \
  vector_red_temp[THREADID]
#define KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(THREADID) \
  vector_red_temp[THREADID]
#define KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_CLAUSE1 \
  private(vector_red_temp [0:chunk_size])
#define KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_CLAUSE2 \
  private(vector_red_temp [0:chunk_size])
#endif

#define KOKKOS_IMPL_OPENACC_CHUNK_SIZE 64

namespace Kokkos::Experimental::Impl {

// primary template: catch-all non-implemented custom reducers
template <
    class Functor, class Reducer, class Policy,
    bool = Kokkos::Impl::FunctorAnalysis<
               Kokkos::Impl::FunctorPatternInterface::REDUCE, Policy,
               typename Functor::functor_type, typename Reducer::value_type>::
               Reducer::has_join_member_function() ||
           !std::is_arithmetic_v<typename Reducer::value_type>>
struct OpenACCParallelReduceHelper {
  OpenACCParallelReduceHelper(Functor const&, Reducer const&, Policy const&,
                              Kokkos::Impl::FunctorAnalysis<
                                  Kokkos::Impl::FunctorPatternInterface::REDUCE,
                                  Policy, typename Functor::functor_type,
                                  typename Reducer::value_type>::pointer_type,
                              bool) {
    static_assert(Kokkos::Impl::always_false<Functor>::value,
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

  static constexpr bool FunctorHasJoin = Impl::FunctorAnalysis<
      Impl::FunctorPatternInterface::REDUCE, Policy, FunctorType,
      typename ReducerType::value_type>::Reducer::has_join_member_function();
  static constexpr bool UseReducer =
      !std::is_same_v<FunctorType, typename ReducerType::functor_type>;

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

    Kokkos::Experimental::Impl::FunctorAdapter<
        FunctorType, Policy, Kokkos::Experimental::Impl::RoutineClause::seq>
        functor(m_functor_reducer.get_functor());
    if constexpr (FunctorHasJoin || !std::is_arithmetic_v<ValueType>) {
      Kokkos::Experimental::Impl::OpenACCParallelReduceHelper(
          functor, m_functor_reducer.get_reducer(), m_policy, m_result_ptr,
          m_result_ptr_on_device);
    } else {
      Kokkos::Experimental::Impl::OpenACCParallelReduceHelper(
          functor,
          std::conditional_t<UseReducer, typename ReducerType::functor_type,
                             Sum<ValueType>>(val),
          m_policy, m_result_ptr, m_result_ptr_on_device);
    }
  }
};

namespace Kokkos::Experimental::Impl {
template <class Policy, class Functor, class Reducer, class Pointer>
void OpenACCParallelReduceCustom(Schedule<Static>, Policy const& apolicy,
                                 Functor const& afunctor,
                                 Reducer const& areducer, Pointer am_result_ptr,
                                 bool m_result_on_device) {
  using ValueType = typename Reducer::value_type;
  using IndexType = typename Policy::index_type;
  auto const policy(apolicy);
  auto const functor(afunctor);
  auto const reducer(areducer);
  auto m_result_ptr(am_result_ptr);
  int const async_arg  = policy.space().acc_async_queue();
  IndexType chunk_size = policy.chunk_size();
  if (chunk_size <= 1) {
    chunk_size = KOKKOS_IMPL_OPENACC_CHUNK_SIZE;
  }
  const IndexType begin    = policy.begin();
  const IndexType end      = policy.end();
  const IndexType N        = end - begin;
  const IndexType n_chunks = (N + chunk_size - 1) / chunk_size;
#ifdef KOKKOS_COMPILER_CLANG
  const IndexType num_elements = n_chunks * chunk_size;
#else
  const IndexType num_elements = chunk_size;
#endif
  Kokkos::View<ValueType, Kokkos::Experimental::OpenACCSpace> m_result_view(
      "Kokkos::OpenACCParallelScan::m_result_view");
  Kokkos::View<ValueType*, Kokkos::Experimental::OpenACCSpace> gang_red_temp(
      "Kokkos::OpenACCParallelReduceCustom::gang_red_temp", n_chunks);
  std::unique_ptr<ValueType[]> vector_red_temp_owner(
      new ValueType[num_elements]);
  ValueType* vector_red_temp = vector_red_temp_owner.get();
  ValueType initVal;
  reducer.init(&initVal);

#pragma acc enter data copyin(functor, reducer, gang_red_temp) async(async_arg)

  /* clang-format off */
KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang num_gangs(n_chunks) num_workers(1) vector_length(chunk_size) KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_CLAUSE1 present(functor, reducer, gang_red_temp) async(async_arg))
  /* clang-format on */
  for (IndexType team_id = 0; team_id < n_chunks; ++team_id) {
    IndexType tSize = chunk_size;
    IndexType tStep;
#pragma acc loop vector
    for (IndexType thread_id = 0; thread_id < chunk_size; ++thread_id) {
      const IndexType local_offset = team_id * chunk_size + begin;
      const IndexType idx          = local_offset + thread_id;
      ValueType temp;
      temp = initVal;
      if (idx < end) {
        functor(idx, temp);
      }
      KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS1(thread_id) = temp;
    }
#pragma acc loop seq
    for (tStep = (chunk_size >> 1); tStep > 0; tStep >>= 1) {
#pragma acc loop vector
      for (IndexType thread_id = 0; thread_id < chunk_size; ++thread_id) {
        if (thread_id < tStep) {
          reducer.join(
              &KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS1(thread_id),
              &KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS1(thread_id + tStep));
        }
        IndexType checkOdd = tSize & 1;
        if (checkOdd == 1) {
          if (thread_id == 0) {
            reducer.join(
                &KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS1(thread_id),
                &KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS1(tSize - 1));
          }
        }
      }
      tSize = tStep;
    }
    gang_red_temp(team_id) = KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS1(0);
  }
  /* clang-format off */
KOKKOS_IMPL_ACC_PRAGMA(parallel num_gangs(1) num_workers(1) vector_length(chunk_size) KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_CLAUSE2 present(reducer, gang_red_temp) copyin(m_result_view) async(async_arg))
  /* clang-format on */
  {
    IndexType tSize = chunk_size;
    IndexType tStep;
#pragma acc loop vector
    for (IndexType thread_id = 0; thread_id < chunk_size; ++thread_id) {
      IndexType idx;
      ValueType temp;
      temp = initVal;
#pragma acc loop seq
      for (idx = thread_id; idx < n_chunks; idx += chunk_size) {
        reducer.join(&temp, &gang_red_temp(idx));
      }
      KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(thread_id) = temp;
    }
#pragma acc loop seq
    for (tStep = (chunk_size >> 1); tStep > 0; tStep >>= 1) {
#pragma acc loop vector
      for (IndexType thread_id = 0; thread_id < chunk_size; ++thread_id) {
        if (thread_id < tStep) {
          reducer.join(
              &KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(thread_id),
              &KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(thread_id + tStep));
        }
        IndexType checkOdd = tSize & 1;
        if (checkOdd == 1) {
          if (thread_id == 0) {
            reducer.join(
                &KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(thread_id),
                &KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(tSize - 1));
          }
        }
      }
      tSize = tStep;
    }
    reducer.final(&KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(0));
    if (m_result_on_device) {
      *m_result_ptr = KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(0);
    } else {
      m_result_view() = KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2(0);
    }
  }

  if (!m_result_on_device) {
    Kokkos::Impl::DeepCopy<HostSpace, Kokkos::Experimental::OpenACCSpace,
                           Kokkos::Experimental::OpenACC>(
        policy.space(), m_result_ptr, m_result_view.data(), sizeof(ValueType));
  }

  acc_wait(async_arg);
#pragma acc exit data delete (functor, reducer, gang_red_temp)
}
}  // namespace Kokkos::Experimental::Impl

template <class Functor, class Reducer, class... Traits>
struct Kokkos::Experimental::Impl::OpenACCParallelReduceHelper<
    Functor, Reducer, Kokkos::RangePolicy<Traits...>, true> {
  using Policy       = RangePolicy<Traits...>;
  using ScheduleType = Kokkos::Experimental::Impl::OpenACCScheduleType<Policy>;
  using Pointer      = Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::REDUCE, Policy,
      typename Functor::functor_type,
      typename Reducer::value_type>::pointer_type;

  OpenACCParallelReduceHelper(Functor const& functor, Reducer const& reducer,
                              Policy const& policy, Pointer m_result_ptr,
                              bool m_result_on_device) {
    auto const begin = policy.begin();
    auto const end   = policy.end();

    if (end <= begin) {
      return;
    }

    OpenACCParallelReduceCustom(ScheduleType(), policy, functor, reducer,
                                m_result_ptr, m_result_on_device);
  }
};

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_SCHEDULE(REDUCER,         \
                                                              OPERATOR)        \
  namespace Kokkos::Experimental::Impl {                                       \
  template <class Policy, class Functor, class Reducer, class Pointer>         \
  void OpenACCParallelReduce##REDUCER(Schedule<Static>, Policy const& apolicy, \
                                      Functor const& afunctor,                 \
                                      Reducer const& areducer,                 \
                                      Pointer am_result_ptr,                   \
                                      bool m_result_on_device) {               \
    using IndexType = typename Policy::index_type;                             \
    using ValueType = typename Reducer::value_type;                            \
    /* FIXME_OPENACC FIXME_NVHPC workaround compiler bug (incorrect scope      \
       analysis)                                                               \
       NVC++-S-1067-Cannot determine bounds for array - functor */             \
    auto const policy(apolicy);                                                \
    auto const functor(afunctor);                                              \
    auto const reducer(areducer);                                              \
    auto const m_result_ptr(am_result_ptr);                                    \
    ValueType val;                                                             \
    int const async_arg        = policy.space().acc_async_queue();             \
    const IndexType chunk_size = policy.chunk_size();                          \
    const IndexType begin      = policy.begin();                               \
    const IndexType end        = policy.end();                                 \
    reducer.init(val);                                                         \
    if (chunk_size > 1) {                                                      \
      /* clang-format off */                                            \
        KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang(static:chunk_size) vector reduction(OPERATOR:val) copyin(functor) async(async_arg))      \
      /* clang-format on */                                                    \
      for (auto i = begin; i < end; i++) {                                     \
        functor(i, val);                                                       \
      }                                                                        \
    } else {                                                                   \
      /* clang-format off */ \
        KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang(static:*) vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                 \
      /* clang-format on */                                                    \
      for (auto i = begin; i < end; i++) {                                     \
        functor(i, val);                                                       \
      }                                                                        \
    }                                                                          \
    acc_wait(async_arg);                                                       \
    if (m_result_on_device) {                                                  \
      acc_memcpy_to_device(m_result_ptr, &val, sizeof(ValueType));             \
    } else {                                                                   \
      *m_result_ptr = val;                                                     \
    }                                                                          \
  }                                                                            \
                                                                               \
  template <class Policy, class Functor, class Reducer, class Pointer>         \
  void OpenACCParallelReduce##REDUCER(Schedule<Dynamic>,                       \
                                      Policy const& apolicy,                   \
                                      Functor const& afunctor,                 \
                                      Reducer const& areducer,                 \
                                      Pointer am_result_ptr,                   \
                                      bool m_result_on_device) {               \
    using IndexType = typename Policy::index_type;                             \
    using ValueType = typename Reducer::value_type;                            \
    auto const policy(apolicy);                                                \
    auto const functor(afunctor);                                              \
    auto const reducer(areducer);                                              \
    auto const m_result_ptr(am_result_ptr);                                    \
    ValueType val;                                                             \
    int const async_arg        = policy.space().acc_async_queue();             \
    const IndexType chunk_size = policy.chunk_size();                          \
    const IndexType begin      = policy.begin();                               \
    const IndexType end        = policy.end();                                 \
    reducer.init(val);                                                         \
    if (chunk_size > 1) {                                                      \
      /* clang-format off */                                            \
        KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang(static:chunk_size) vector reduction(OPERATOR:val) copyin(functor) async(async_arg))      \
      /* clang-format on */                                                    \
      for (auto i = begin; i < end; i++) {                                     \
        functor(i, val);                                                       \
      }                                                                        \
    } else {                                                                   \
      /* clang-format off */ \
        KOKKOS_IMPL_ACC_PRAGMA(parallel loop gang vector reduction(OPERATOR:val) copyin(functor) async(async_arg))                                                 \
      /* clang-format on */                                                    \
      for (auto i = begin; i < end; i++) {                                     \
        functor(i, val);                                                       \
      }                                                                        \
    }                                                                          \
    acc_wait(async_arg);                                                       \
    if (m_result_on_device) {                                                  \
      acc_memcpy_to_device(m_result_ptr, &val, sizeof(ValueType));             \
    } else {                                                                   \
      *m_result_ptr = val;                                                     \
    }                                                                          \
  }                                                                            \
  }  // namespace Kokkos::Experimental::Impl

#define KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_HELPER(REDUCER, OPERATOR)          \
  KOKKOS_IMPL_OPENACC_PARALLEL_REDUCE_DISPATCH_SCHEDULE(REDUCER, OPERATOR)     \
  template <class Functor, class Scalar, class Space, class... Traits>         \
  struct Kokkos::Experimental::Impl::OpenACCParallelReduceHelper<              \
      Functor, Kokkos::REDUCER<Scalar, Space>, Kokkos::RangePolicy<Traits...>, \
      false> {                                                                 \
    using Policy = RangePolicy<Traits...>;                                     \
    using ScheduleType =                                                       \
        Kokkos::Experimental::Impl::OpenACCScheduleType<Policy>;               \
    using Reducer = REDUCER<Scalar, Space>;                                    \
    using Pointer = Kokkos::Impl::FunctorAnalysis<                             \
        Kokkos::Impl::FunctorPatternInterface::REDUCE, Policy,                 \
        typename Functor::functor_type,                                        \
        typename Reducer::value_type>::pointer_type;                           \
                                                                               \
    OpenACCParallelReduceHelper(Functor const& functor,                        \
                                Reducer const& reducer, Policy const& policy,  \
                                Pointer m_result_ptr,                          \
                                bool m_result_on_device) {                     \
      auto const begin = policy.begin();                                       \
      auto const end   = policy.end();                                         \
                                                                               \
      if (end <= begin) {                                                      \
        return;                                                                \
      }                                                                        \
                                                                               \
      OpenACCParallelReduce##REDUCER(ScheduleType(), policy, functor, reducer, \
                                     m_result_ptr, m_result_on_device);        \
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

#undef KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS1
#undef KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_ACCESS2
#undef KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_CLAUSE1
#undef KOKKOS_IMPL_OPENACC_VECTOR_RED_TEMP_CLAUSE2
#undef KOKKOS_IMPL_OPENACC_CHUNK_SIZE

#endif
