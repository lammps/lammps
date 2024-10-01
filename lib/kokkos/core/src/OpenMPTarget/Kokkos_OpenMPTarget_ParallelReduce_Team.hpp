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

#ifndef KOKKOS_OPENMPTARGET_PARALLELREDUCE_TEAM_HPP
#define KOKKOS_OPENMPTARGET_PARALLELREDUCE_TEAM_HPP

#include <omp.h>
#include <sstream>
#include <Kokkos_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Parallel.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTarget_Parallel_Common.hpp>

// FIXME_OPENMPTARGET - Using this macro to implement a workaround for
// hierarchical reducers. It avoids hitting the code path which we wanted to
// write but doesn't work. undef'ed at the end.
// Intel compilers prefer the non-workaround version.
#ifndef KOKKOS_ARCH_INTEL_GPU
#define KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND
#endif

namespace Kokkos {

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team
 * and a summation of val is performed and put into result.
 */

template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer<ValueType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ValueType& result) {
  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamThread_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp barrier
  TeamThread_scratch[0] = ValueType();
#pragma omp barrier

  if constexpr (std::is_arithmetic<ValueType>::value) {
#pragma omp for reduction(+ : TeamThread_scratch[:1])
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      TeamThread_scratch[0] += tmp;
    }
  } else {
#pragma omp declare reduction(custom:ValueType : omp_out += omp_in)

#pragma omp for reduction(custom : TeamThread_scratch[:1])
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      TeamThread_scratch[0] += tmp;
    }
  }

  result = TeamThread_scratch[0];
}

#if !defined(KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND)
// For some reason the actual version we wanted to write doesn't work
// and crashes. We should try this with every new compiler
// This is the variant we actually wanted to write
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType result) {
  using ValueType = typename ReducerType::value_type;

#pragma omp declare reduction(                                               \
    custominner:ValueType                                                    \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamThread_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp barrier
  Impl::OpenMPTargetReducerWrapper<ReducerType>::init(TeamThread_scratch[0]);
#pragma omp barrier

#pragma omp for reduction(custominner : TeamThread_scratch[:1])
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, TeamThread_scratch[0]);
  }
  result.reference() = TeamThread_scratch[0];
}
#else
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType result) {
  using ValueType = typename ReducerType::value_type;

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamThread_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp declare reduction(                                               \
    omp_red_teamthread_reducer:ValueType                                     \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#pragma omp barrier
  ValueType tmp;
  result.init(tmp);
  TeamThread_scratch[0] = tmp;
#pragma omp barrier

  iType team_size = iType(omp_get_num_threads());
#pragma omp for reduction(omp_red_teamthread_reducer \
                          : TeamThread_scratch[:1]) schedule(static, 1)
  for (iType t = 0; t < team_size; t++) {
    ValueType tmp2;
    result.init(tmp2);

    for (iType i = loop_boundaries.start + t; i < loop_boundaries.end;
         i += team_size) {
      lambda(i, tmp2);
    }

    // FIXME_OPENMPTARGET: Join should work but doesn't. Every threads gets a
    // private TeamThread_scratch[0] and at the end of the for-loop the `join`
    // operation is performed by OpenMP itself and hence the simple assignment
    // works.
    //    result.join(TeamThread_scratch[0], tmp2);
    TeamThread_scratch[0] = tmp2;
  }

  result.reference() = TeamThread_scratch[0];
}
#endif  // KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a reduction of val is performed using JoinType(ValueType& val, const
 * ValueType& update) and put into init_result. The input value of init_result
 * is used as initializer for temporary variables of ValueType. Therefore the
 * input value should be the neutral element with respect to the join operation
 * (e.g. '0 for +-' or '1 for *').
 */
template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& init_result) {
  ValueType* TeamThread_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  // FIXME_OPENMPTARGET: Still need to figure out how to get value_count here.
  const int value_count = 1;

#pragma omp barrier
  TeamThread_scratch[0] = init_result;
#pragma omp barrier

#pragma omp for
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, TeamThread_scratch[omp_get_num_threads() * value_count]);
  }

  // Reduce all partial results within a team.
  const int team_size      = omp_get_num_threads();
  int tree_neighbor_offset = 1;
  do {
#pragma omp for
    for (int i = 0; i < team_size - tree_neighbor_offset;
         i += 2 * tree_neighbor_offset) {
      const int neighbor = i + tree_neighbor_offset;
      join(lambda, &TeamThread_scratch[i * value_count],
           &TeamThread_scratch[neighbor * value_count]);
    }
    tree_neighbor_offset *= 2;
  } while (tree_neighbor_offset < team_size);
  init_result = TeamThread_scratch[0];
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a summation of val is performed and put into result.
 */
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda, ValueType& result) {
  ValueType vector_reduce = ValueType();

  if constexpr (std::is_arithmetic<ValueType>::value) {
#pragma omp simd reduction(+ : vector_reduce)
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      vector_reduce += tmp;
    }
  } else {
#pragma omp declare reduction(custom:ValueType : omp_out += omp_in)

#pragma omp simd reduction(custom : vector_reduce)
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      lambda(i, vector_reduce);
    }
  }

  result = vector_reduce;
}

template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType const& result) {
  using ValueType = typename ReducerType::value_type;

#pragma omp declare reduction(                                               \
    custom:ValueType                                                         \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

  ValueType vector_reduce;
  Impl::OpenMPTargetReducerWrapper<ReducerType>::init(vector_reduce);

#pragma omp simd reduction(custom : vector_reduce)
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, vector_reduce);
  }

  result.reference() = vector_reduce;
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a reduction of val is performed using JoinType(ValueType& val, const
 * ValueType& update) and put into init_result. The input value of init_result
 * is used as initializer for temporary variables of ValueType. Therefore the
 * input value should be the neutral element with respect to the join operation
 * (e.g. '0 for +-' or '1 for *').
 */
template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& init_result) {
  ValueType result = init_result;

  // FIXME_OPENMPTARGET think about omp simd
  // join does not work with omp reduction clause
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    ValueType tmp = ValueType();
    lambda(i, tmp);
    join(result, tmp);
  }

  init_result = result;
}

/** \brief  Intra-team vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling team
 * and a summation of val is performed and put into result.
 */
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamVectorRangeBoundariesStruct<
        iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
    const Lambda& lambda, ValueType& result) {
  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamVector_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp barrier
  TeamVector_scratch[0] = ValueType();
#pragma omp barrier

  if constexpr (std::is_arithmetic<ValueType>::value) {
#pragma omp for simd reduction(+ : TeamVector_scratch[:1])
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      TeamVector_scratch[0] += tmp;
    }
  } else {
#pragma omp declare reduction(custom:ValueType : omp_out += omp_in)

#pragma omp for simd reduction(custom : TeamVector_scratch[:1])
    for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
      ValueType tmp = ValueType();
      lambda(i, tmp);
      TeamVector_scratch[0] += tmp;
    }
  }

  result = TeamVector_scratch[0];
}

#if !defined(KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND)
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType const& result) {
  using ValueType = typename ReducerType::value_type;

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

#pragma omp declare reduction(                                               \
    custom:ValueType                                                         \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

  ValueType* TeamVector_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp barrier
  Impl::OpenMPTargetReducerWrapper<ReducerType>::init(TeamVector_scratch[0]);
#pragma omp barrier

#pragma omp for simd reduction(custom : TeamVector_scratch[:1])
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i, TeamVector_scratch[0]);
  }

  result.reference() = TeamVector_scratch[0];
}
#else
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                    iType, Impl::OpenMPTargetExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ReducerType const& result) {
  using ValueType = typename ReducerType::value_type;

  // FIXME_OPENMPTARGET - Make sure that if its an array reduction, number of
  // elements in the array <= 32. For reduction we allocate, 16 bytes per
  // element in the scratch space, hence, 16*32 = 512.
  static_assert(sizeof(ValueType) <=
                Impl::OpenMPTargetExecTeamMember::TEAM_REDUCE_SIZE);

  ValueType* TeamVector_scratch =
      static_cast<ValueType*>(loop_boundaries.team.impl_reduce_scratch());

#pragma omp declare reduction(                                               \
    omp_red_teamthread_reducer:ValueType                                     \
    : Impl::OpenMPTargetReducerWrapper <ReducerType>::join(omp_out, omp_in)) \
    initializer(                                                             \
        Impl::OpenMPTargetReducerWrapper <ReducerType>::init(omp_priv))

#pragma omp barrier
  ValueType tmp;
  result.init(tmp);
  TeamVector_scratch[0] = tmp;
#pragma omp barrier

  iType team_size = iType(omp_get_num_threads());
#pragma omp for simd reduction(omp_red_teamthread_reducer \
                               : TeamVector_scratch[:1]) schedule(static, 1)
  for (iType t = 0; t < team_size; t++) {
    ValueType tmp2;
    result.init(tmp2);

    for (iType i = loop_boundaries.start + t; i < loop_boundaries.end;
         i += team_size) {
      lambda(i, tmp2);
    }
    TeamVector_scratch[0] = tmp2;
  }

  result.reference() = TeamVector_scratch[0];
}
#endif  // KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND

namespace Impl {

template <class CombinedFunctorReducerType, class... Properties>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::TeamPolicy<Properties...>,
                     Kokkos::Experimental::OpenMPTarget> {
 private:
  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::Experimental::OpenMPTarget,
                                       Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using WorkTag = typename Policy::work_tag;
  using Member  = typename Policy::member_type;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;
  using value_type     = typename ReducerType::value_type;

  bool m_result_ptr_on_device;
  const int m_result_ptr_num_elems;

  static constexpr bool FunctorHasJoin = Impl::FunctorAnalysis<
      Impl::FunctorPatternInterface::REDUCE, Policy, FunctorType,
      typename ReducerType::value_type>::Reducer::has_join_member_function();
  static constexpr bool UseReducer =
      !std::is_same_v<FunctorType, typename ReducerType::functor_type>;
  static constexpr bool IsArray = std::is_pointer_v<reference_type>;

  using ParReduceSpecialize =
      ParallelReduceSpecialize<FunctorType, Policy,
                               typename ReducerType::functor_type, pointer_type,
                               typename ReducerType::value_type>;

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const pointer_type m_result_ptr;
  const size_t m_shmem_size;

 public:
  void execute() const {
    // Only let one ParallelReduce instance at a time use the scratch memory.
    std::scoped_lock<std::mutex> scratch_memory_lock(
        OpenMPTargetExec::m_mutex_scratch_ptr);
    const FunctorType& functor = m_functor_reducer.get_functor();
    if constexpr (FunctorHasJoin) {
      ParReduceSpecialize::execute_init_join(functor, m_policy, m_result_ptr,
                                             m_result_ptr_on_device);
    } else if constexpr (UseReducer) {
      ParReduceSpecialize::execute_reducer(functor, m_policy, m_result_ptr,
                                           m_result_ptr_on_device);
    } else if constexpr (IsArray) {
      if (m_result_ptr_num_elems <= 2) {
        ParReduceSpecialize::template execute_array<2>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 4) {
        ParReduceSpecialize::template execute_array<4>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 8) {
        ParReduceSpecialize::template execute_array<8>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 16) {
        ParReduceSpecialize::template execute_array<16>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else if (m_result_ptr_num_elems <= 32) {
        ParReduceSpecialize::template execute_array<32>(
            functor, m_policy, m_result_ptr, m_result_ptr_on_device);
      } else {
        Kokkos::abort("array reduction length must be <= 32");
      }
    } else {
      ParReduceSpecialize::template execute_array<1>(
          functor, m_policy, m_result_ptr, m_result_ptr_on_device);
    }
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const Policy& arg_policy, const ViewType& arg_result)
      : m_result_ptr_on_device(
            MemorySpaceAccess<Kokkos::Experimental::OpenMPTargetSpace,
                              typename ViewType::memory_space>::accessible),
        m_result_ptr_num_elems(arg_result.size()),
        m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_result_ptr(arg_result.data()),
        m_shmem_size(
            arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
            FunctorTeamShmemSize<FunctorType>::value(
                arg_functor_reducer.get_functor(), arg_policy.team_size())) {}
};

}  // namespace Impl

#ifdef KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND
#undef KOKKOS_IMPL_HIERARCHICAL_REDUCERS_WORKAROUND
#endif
}  // namespace Kokkos

#endif
