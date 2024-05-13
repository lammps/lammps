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

#ifndef KOKKOS_IMPL_OPENMP_TASK_HPP
#define KOKKOS_IMPL_OPENMP_TASK_HPP

#if defined(KOKKOS_ENABLE_TASKPOLICY)

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <>
class TaskQueueSpecialization<Kokkos::Experimental::OpenMPTarget> {
 public:
  using execution_space = Kokkos::Experimental::OpenMPTarget;
  using queue_type      = Kokkos::Impl::TaskQueue<execution_space>;
  using task_base_type  = Kokkos::Impl::TaskBase<execution_space, void, void>;

  // Must specify memory space
  using memory_space = Kokkos::HostSpace;

  static void iff_single_thread_recursive_execute(queue_type* const);

  // Must provide task queue execution function
  static void execute(queue_type* const);

  // Must provide mechanism to set function pointer in
  // execution space from the host process.
  template <typename FunctorType>
  static void proc_set_apply(task_base_type::function_type* ptr) {
    using TaskType = TaskBase<Kokkos::Experimental::OpenMPTarget,
                              typename FunctorType::value_type, FunctorType>;
    *ptr           = TaskType::apply;
  }
};

extern template class TaskQueue<Kokkos::Experimental::OpenMPTarget>;

//----------------------------------------------------------------------------

template <>
class TaskExec<Kokkos::Experimental::OpenMPTarget> {
 private:
  TaskExec(TaskExec&&)      = delete;
  TaskExec(TaskExec const&) = delete;
  TaskExec& operator=(TaskExec&&) = delete;
  TaskExec& operator=(TaskExec const&) = delete;

  using PoolExec = Kokkos::Impl::OpenMPTargetExec;

  friend class Kokkos::Impl::TaskQueue<Kokkos::Experimental::OpenMPTarget>;
  friend class Kokkos::Impl::TaskQueueSpecialization<
      Kokkos::Experimental::OpenMPTarget>;

  PoolExec* const m_self_exec;  ///< This thread's thread pool data structure
  PoolExec* const m_team_exec;  ///< Team thread's thread pool data structure
  int64_t m_sync_mask;
  int64_t mutable m_sync_value;
  int mutable m_sync_step;
  int m_group_rank;  ///< Which "team" subset of thread pool
  int m_team_rank;   ///< Which thread within a team
  int m_team_size;

  TaskExec();
  TaskExec(PoolExec& arg_exec, int arg_team_size);

  void team_barrier_impl() const;

 public:
  KOKKOS_FUNCTION void* team_shared() const {
    KOKKOS_IF_ON_HOST(
        (return m_team_exec ? m_team_exec->scratch_thread() : nullptr;))

    KOKKOS_IF_ON_DEVICE((return nullptr;))
  }

  KOKKOS_FUNCTION int team_shared_size() const {
    KOKKOS_IF_ON_HOST(
        (return m_team_exec ? m_team_exec->scratch_thread_size() : 0;))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  /**\brief  Whole team enters this function call
   *         before any teeam member returns from
   *         this function call.
   */
  KOKKOS_FUNCTION void team_barrier() const {
    KOKKOS_IF_ON_HOST((if (1 < m_team_size) { team_barrier_impl(); }))
  }

  KOKKOS_INLINE_FUNCTION
  int team_rank() const { return m_team_rank; }

  KOKKOS_INLINE_FUNCTION
  int team_size() const { return m_team_size; }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template <typename iType>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >
TeamThreadRange(Impl::TaskExec<Kokkos::Experimental::OpenMPTarget>& thread,
                const iType& count) {
  return Impl::TeamThreadRangeBoundariesStruct<
      iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >(thread,
                                                                  count);
}

template <typename iType>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >
TeamThreadRange(Impl::TaskExec<Kokkos::Experimental::OpenMPTarget>& thread,
                const iType& start, const iType& end) {
  return Impl::TeamThreadRangeBoundariesStruct<
      iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >(thread, start,
                                                                  end);
}

/** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team.
 */
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >&
        loop_boundaries,
    const Lambda& lambda) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i);
  }
}

template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >&
        loop_boundaries,
    const Lambda& lambda, ValueType& initialized_result) {
  int team_rank =
      loop_boundaries.thread.team_rank();  // member num within the team
  ValueType result = initialized_result;

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }

  if (1 < loop_boundaries.thread.team_size()) {
    ValueType* shared = (ValueType*)loop_boundaries.thread.team_shared();

    loop_boundaries.thread.team_barrier();
    shared[team_rank] = result;

    loop_boundaries.thread.team_barrier();

    // reduce across threads to thread 0
    if (team_rank == 0) {
      for (int i = 1; i < loop_boundaries.thread.team_size(); i++) {
        shared[0] += shared[i];
      }
    }

    loop_boundaries.thread.team_barrier();

    // broadcast result
    initialized_result = shared[0];
  } else {
    initialized_result = result;
  }
}

template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >&
        loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& initialized_result) {
  int team_rank =
      loop_boundaries.thread.team_rank();  // member num within the team
  ValueType result = initialized_result;

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }

  if (1 < loop_boundaries.thread.team_size()) {
    ValueType* shared = (ValueType*)loop_boundaries.thread.team_shared();

    loop_boundaries.thread.team_barrier();
    shared[team_rank] = result;

    loop_boundaries.thread.team_barrier();

    // reduce across threads to thread 0
    if (team_rank == 0) {
      for (int i = 1; i < loop_boundaries.thread.team_size(); i++) {
        join(shared[0], shared[i]);
      }
    }

    loop_boundaries.thread.team_barrier();

    // broadcast result
    initialized_result = shared[0];
  } else {
    initialized_result = result;
  }
}

// placeholder for future function
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >&
        loop_boundaries,
    const Lambda& lambda, ValueType& initialized_result) {}

// placeholder for future function
template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >&
        loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& initialized_result) {
}

template <typename ValueType, typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >&
        loop_boundaries,
    const Lambda& lambda) {
  ValueType accum = 0;
  ValueType val, local_total;
  ValueType* shared = (ValueType*)loop_boundaries.thread.team_shared();
  int team_size     = loop_boundaries.thread.team_size();
  int team_rank =
      loop_boundaries.thread.team_rank();  // member num within the team

  // Intra-member scan
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    local_total = 0;
    lambda(i, local_total, false);
    val = accum;
    lambda(i, val, true);
    accum += local_total;
  }

  shared[team_rank] = accum;
  loop_boundaries.thread.team_barrier();

  // Member 0 do scan on accumulated totals
  if (team_rank == 0) {
    for (iType i = 1; i < team_size; i += 1) {
      shared[i] += shared[i - 1];
    }
    accum = 0;  // Member 0 set accum to 0 in preparation for inter-member scan
  }

  loop_boundaries.thread.team_barrier();

  // Inter-member scan adding in accumulated totals
  if (team_rank != 0) {
    accum = shared[team_rank - 1];
  }
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    local_total = 0;
    lambda(i, local_total, false);
    val = accum;
    lambda(i, val, true);
    accum += local_total;
  }
}

// placeholder for future function
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Experimental::OpenMPTarget> >&
        loop_boundaries,
    const Lambda& lambda) {}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKPOLICY ) */
#endif /* #ifndef KOKKOS_IMPL_OPENMP_TASK_HPP */
