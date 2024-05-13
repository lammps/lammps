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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_TASKSCHEDULER_FWD_HPP
#define KOKKOS_TASKSCHEDULER_FWD_HPP

//----------------------------------------------------------------------------

#include <cstddef>
#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_Core_fwd.hpp>
//----------------------------------------------------------------------------

namespace Kokkos {

// Forward declarations used in Impl::TaskQueue

template <typename ValueType, typename Scheduler>
class BasicFuture;

template <class Space, class Queue>
class SimpleTaskScheduler;

template <class Space, class Queue>
class BasicTaskScheduler;

template <typename Space>
struct is_scheduler : public std::false_type {};

template <class Space, class Queue>
struct is_scheduler<BasicTaskScheduler<Space, Queue>> : public std::true_type {
};

template <class Space, class Queue>
struct is_scheduler<SimpleTaskScheduler<Space, Queue>> : public std::true_type {
};

enum class TaskPriority : int { High = 0, Regular = 1, Low = 2 };

}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

template <class Device>
class MemoryPool;

namespace Impl {

template <class TaskQueueTraits>
class TaskNode;

class TaskBase;

/*\brief  Implementation data for task data management, access, and execution.
 *  (Deprecated)
 *  CRTP Inheritance structure to allow static_cast from the
 *  task root type and a task's FunctorType.
 *
 *    TaskBase< Space , ResultType , FunctorType >
 *      : TaskBase< Space , ResultType , void >
 *      , FunctorType
 *      { ... };
 *
 *    TaskBase< Space , ResultType , void >
 *      : TaskBase< Space , void , void >
 *      { ... };
 */
template <typename Space, typename ResultType, typename FunctorType>
class Task;

class TaskQueueBase;

template <typename Space, typename MemorySpace>
class TaskQueue;

template <typename ExecSpace, typename MemorySpace>
class TaskQueueMultiple;

template <typename ExecSpace, typename MemSpace, typename TaskQueueTraits,
          class MemoryPool =
              Kokkos::MemoryPool<Kokkos::Device<ExecSpace, MemSpace>>>
class SingleTaskQueue;

template <typename ExecSpace, typename MemSpace, typename TaskQueueTraits,
          class MemoryPool>
class MultipleTaskQueue;

struct TaskQueueTraitsLockBased;

template <size_t CircularBufferSize = 64>
struct TaskQueueTraitsChaseLev;

template <typename ResultType>
struct TaskResult;

struct TaskSchedulerBase;

template <class ExecSpace>
struct default_tasking_memory_space_for_execution_space {
  using type = typename ExecSpace::memory_space;
};

#if defined(KOKKOS_ENABLE_CUDA)
template <>
struct default_tasking_memory_space_for_execution_space<Kokkos::Cuda> {
  using type = Kokkos::CudaUVMSpace;
};
#endif

template <class ExecSpace>
using default_tasking_memory_space_for_execution_space_t =
    typename default_tasking_memory_space_for_execution_space<ExecSpace>::type;

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

template <typename Space>
using DeprecatedTaskScheduler = BasicTaskScheduler<
    Space,
    Impl::TaskQueue<
        Space,
        Impl::default_tasking_memory_space_for_execution_space_t<Space>>>;

template <typename Space>
using DeprecatedTaskSchedulerMultiple = BasicTaskScheduler<
    Space,
    Impl::TaskQueueMultiple<
        Space,
        Impl::default_tasking_memory_space_for_execution_space_t<Space>>>;

template <typename Space>
using TaskScheduler = SimpleTaskScheduler<
    Space,
    Impl::SingleTaskQueue<
        Space, Impl::default_tasking_memory_space_for_execution_space_t<Space>,
        Impl::TaskQueueTraitsLockBased>>;

template <typename Space>
using TaskSchedulerMultiple = SimpleTaskScheduler<
    Space,
    Impl::MultipleTaskQueue<
        Space, Impl::default_tasking_memory_space_for_execution_space_t<Space>,
        Impl::TaskQueueTraitsLockBased,
        Kokkos::MemoryPool<Kokkos::Device<
            Space,
            Impl::default_tasking_memory_space_for_execution_space_t<Space>>>>>;

template <typename Space>
using ChaseLevTaskScheduler = SimpleTaskScheduler<
    Space,
    Impl::MultipleTaskQueue<
        Space, Impl::default_tasking_memory_space_for_execution_space_t<Space>,
        Impl::TaskQueueTraitsChaseLev<>,
        Kokkos::MemoryPool<Kokkos::Device<
            Space,
            Impl::default_tasking_memory_space_for_execution_space_t<Space>>>>>;

template <class Space, class QueueType>
void wait(BasicTaskScheduler<Space, QueueType> const&);

namespace Impl {

struct TaskSchedulerBase {};

class TaskQueueBase {};

template <typename Scheduler, typename EnableIfConstraint = void>
class TaskQueueSpecializationConstrained {};

template <typename Scheduler>
struct TaskQueueSpecialization : TaskQueueSpecializationConstrained<Scheduler> {
};

template <int, typename>
struct TaskPolicyData;

}  // end namespace Impl

}  // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_TASKSCHEDULER_FWD_HPP */
