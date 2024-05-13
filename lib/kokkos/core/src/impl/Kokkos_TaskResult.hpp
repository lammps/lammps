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

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_IMPL_TASKRESULT_HPP
#define KOKKOS_IMPL_TASKRESULT_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_TaskScheduler_fwd.hpp>
#include <Kokkos_Core_fwd.hpp>

#include <impl/Kokkos_TaskBase.hpp>
#include <impl/Kokkos_TaskNode.hpp>

#include <string>
#include <typeinfo>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <typename ResultType>
struct TaskResult {
  enum : int32_t { size = sizeof(ResultType) };

  using reference_type = ResultType&;

  template <class CountType>
  KOKKOS_INLINE_FUNCTION static ResultType* ptr(
      PoolAllocatedObjectBase<CountType>* task) {
    return reinterpret_cast<ResultType*>(reinterpret_cast<char*>(task) +
                                         task->get_allocation_size() -
                                         sizeof(ResultType));
  }

  KOKKOS_INLINE_FUNCTION static ResultType* ptr(TaskBase* task) {
    return reinterpret_cast<ResultType*>(reinterpret_cast<char*>(task) +
                                         task->m_alloc_size -
                                         sizeof(ResultType));
  }

  KOKKOS_INLINE_FUNCTION static reference_type get(TaskBase* task) {
    return *ptr(task);
  }

  template <class TaskQueueTraits>
  KOKKOS_INLINE_FUNCTION static reference_type get(
      TaskNode<TaskQueueTraits>* task) {
    return *ptr(task);
  }

  KOKKOS_INLINE_FUNCTION static void destroy(TaskBase* task) {
    get(task).~ResultType();
  }

  // template <class TaskQueueTraits>
  // KOKKOS_INLINE_FUNCTION static
  // void destroy( TaskNode<TaskQueueTraits>* task )
  //{ get(task).~ResultType(); }
};

template <>
struct TaskResult<void> {
  enum : int32_t { size = 0 };

  using reference_type = void;

  template <class TaskQueueTraits>
  KOKKOS_INLINE_FUNCTION static void* ptr(TaskNode<TaskQueueTraits>* /*task*/) {
    return nullptr;
  }

  KOKKOS_INLINE_FUNCTION static void* ptr(TaskBase*) { return nullptr; }

  template <class TaskQueueTraits>
  KOKKOS_INLINE_FUNCTION static reference_type get(
      TaskNode<TaskQueueTraits>* /*task*/) { /* Should never be called */
  }

  KOKKOS_INLINE_FUNCTION static reference_type get(TaskBase*) {}

  KOKKOS_INLINE_FUNCTION static void destroy(TaskBase* /*task*/) {}

  // template <class TaskQueueTraits>
  // KOKKOS_INLINE_FUNCTION static
  // void destroy( TaskNode<TaskQueueTraits>* task )
  //{ }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_TASKRESULT_HPP */
