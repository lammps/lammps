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

#ifndef KOKKOS_IMPL_TASKPOLICYDATA_HPP
#define KOKKOS_IMPL_TASKPOLICYDATA_HPP

//----------------------------------------------------------------------------

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_TaskScheduler_fwd.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template <int TaskEnum, typename DepFutureType>
struct TaskPolicyWithPredecessor {
 private:
  DepFutureType m_predecessor;
  Kokkos::TaskPriority m_priority;

 public:
  KOKKOS_INLINE_FUNCTION
  TaskPolicyWithPredecessor(DepFutureType arg_predecessor,
                            Kokkos::TaskPriority arg_priority)
      : m_predecessor(std::move(arg_predecessor)), m_priority(arg_priority) {}

  TaskPolicyWithPredecessor() = delete;

  KOKKOS_DEFAULTED_FUNCTION
  TaskPolicyWithPredecessor(TaskPolicyWithPredecessor const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskPolicyWithPredecessor(TaskPolicyWithPredecessor&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskPolicyWithPredecessor& operator=(TaskPolicyWithPredecessor const&) =
      default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskPolicyWithPredecessor& operator=(TaskPolicyWithPredecessor&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~TaskPolicyWithPredecessor() = default;

  KOKKOS_INLINE_FUNCTION
  DepFutureType&& predecessor() && { return std::move(m_predecessor); }

  KOKKOS_INLINE_FUNCTION
  constexpr TaskPriority priority() const { return m_priority; }

  KOKKOS_INLINE_FUNCTION
  static constexpr int task_type() noexcept { return TaskEnum; }
};

// TODO @tasking @cleanup DSH clean this up. Using nullptr_t here is too clever
template <int TaskEnum, typename Scheduler,
          typename PredecessorFuture = std::nullptr_t>
struct TaskPolicyWithScheduler {
 public:
  using predecessor_future_type = PredecessorFuture;

 private:
  Scheduler m_scheduler;
  Kokkos::TaskPriority m_priority;
  predecessor_future_type m_predecessor;

 public:
  KOKKOS_INLINE_FUNCTION
  TaskPolicyWithScheduler(Scheduler arg_scheduler,
                          Kokkos::TaskPriority arg_priority)
      : m_scheduler(std::move(arg_scheduler)), m_priority(arg_priority) {}

  KOKKOS_INLINE_FUNCTION
  TaskPolicyWithScheduler(Scheduler arg_scheduler,
                          predecessor_future_type arg_predecessor,
                          Kokkos::TaskPriority arg_priority)
      : m_scheduler(std::move(arg_scheduler)),
        m_priority(arg_priority),
        m_predecessor(std::move(arg_predecessor)) {}

  TaskPolicyWithScheduler() = delete;

  KOKKOS_DEFAULTED_FUNCTION
  TaskPolicyWithScheduler(TaskPolicyWithScheduler const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskPolicyWithScheduler(TaskPolicyWithScheduler&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskPolicyWithScheduler& operator=(TaskPolicyWithScheduler const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskPolicyWithScheduler& operator=(TaskPolicyWithScheduler&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  ~TaskPolicyWithScheduler() = default;

  KOKKOS_INLINE_FUNCTION
  Scheduler& scheduler() & { return m_scheduler; }

  KOKKOS_INLINE_FUNCTION
  constexpr TaskPriority priority() const { return m_priority; }

  KOKKOS_INLINE_FUNCTION
  predecessor_future_type& predecessor() & { return m_predecessor; }

  KOKKOS_INLINE_FUNCTION
  static constexpr bool has_predecessor() noexcept {
    return !std::is_same<PredecessorFuture, std::nullptr_t>::value;
  }

  KOKKOS_INLINE_FUNCTION
  static constexpr int task_type() noexcept { return TaskEnum; }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_TASKPOLICYDATA_HPP */
