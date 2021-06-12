/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_IMPL_TASKBASE_HPP
#define KOKKOS_IMPL_TASKBASE_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_TaskScheduler_fwd.hpp>
#include <Kokkos_Core_fwd.hpp>

#include <impl/Kokkos_LIFO.hpp>

#include <string>
#include <typeinfo>
#include <stdexcept>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Base class for task management, access, and execution.
 *
 *  Inheritance structure to allow static_cast from the task root type
 *  and a task's FunctorType.
 *
 *    // Enable a functor to access the base class
 *    // and provide memory for result value.
 *    TaskBase< Space , ResultType , FunctorType >
 *      : TaskBase< void , void , void >
 *      , FunctorType
 *      { ... };
 *    Followed by memory allocated for result value.
 *
 *
 *  States of a task:
 *
 *    Constructing State, NOT IN a linked list
 *      m_wait == 0
 *      m_next == 0
 *
 *    Scheduling transition : Constructing -> Waiting
 *      before:
 *        m_wait == 0
 *        m_next == this task's initial dependence, 0 if none
 *      after:
 *        m_wait == EndTag
 *        m_next == EndTag
 *
 *    Waiting State, IN a linked list
 *      m_apply != 0
 *      m_queue != 0
 *      m_ref_count > 0
 *      m_wait == head of linked list of tasks waiting on this task
 *      m_next == next of linked list of tasks
 *
 *    transition : Waiting -> Executing
 *      before:
 *        m_next == EndTag
 *      after::
 *        m_next == LockTag
 *
 *    Executing State, NOT IN a linked list
 *      m_apply != 0
 *      m_queue != 0
 *      m_ref_count > 0
 *      m_wait == head of linked list of tasks waiting on this task
 *      m_next == LockTag
 *
 *    Respawn transition : Executing -> Executing-Respawn
 *      before:
 *        m_next == LockTag
 *      after:
 *        m_next == this task's updated dependence, 0 if none
 *
 *    Executing-Respawn State, NOT IN a linked list
 *      m_apply != 0
 *      m_queue != 0
 *      m_ref_count > 0
 *      m_wait == head of linked list of tasks waiting on this task
 *      m_next == this task's updated dependence, 0 if none
 *
 *    transition : Executing -> Complete
 *      before:
 *        m_wait == head of linked list
 *      after:
 *        m_wait == LockTag
 *
 *    Complete State, NOT IN a linked list
 *      m_wait == LockTag: cannot add dependence (<=> complete)
 *      m_next == LockTag: not a member of a wait queue
 *
 */
class TaskBase {
 public:
  enum : int16_t { TaskTeam = 0, TaskSingle = 1, Aggregate = 2 };
  enum : uintptr_t { LockTag = ~uintptr_t(0), EndTag = ~uintptr_t(1) };

  template <typename, typename>
  friend class Kokkos::BasicTaskScheduler;

  using queue_type = TaskQueueBase;

  using function_type = void (*)(TaskBase*, void*);
  using destroy_type  = void (*)(TaskBase*);

  // sizeof(TaskBase) == 48

  function_type m_apply = nullptr;  ///< Apply function pointer
  queue_type* m_queue   = nullptr;  ///< Pointer to the scheduler
  TaskBase* m_next      = nullptr;  ///< next in linked list of ready tasks
  TaskBase* m_wait      = nullptr;  ///< Queue of tasks waiting on this
  int32_t m_ref_count   = 0;
  int32_t m_alloc_size  = 0;
  int32_t m_dep_count;  ///< Aggregate's number of dependences
  int16_t m_task_type;  ///< Type of task
  int16_t m_priority;   ///< Priority of runnable task

  TaskBase(TaskBase&&)      = delete;
  TaskBase(const TaskBase&) = delete;
  TaskBase& operator=(TaskBase&&) = delete;
  TaskBase& operator=(const TaskBase&) = delete;

  KOKKOS_DEFAULTED_FUNCTION ~TaskBase() = default;

  KOKKOS_INLINE_FUNCTION constexpr TaskBase()
      : m_apply(nullptr),
        m_queue(nullptr),
        m_next(nullptr),
        m_wait(nullptr),
        m_ref_count(0),
        m_alloc_size(0),
        m_dep_count(0),
        m_task_type(0),
        m_priority(0) {}

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  TaskBase* volatile* aggregate_dependences() volatile {
    return reinterpret_cast<TaskBase* volatile*>(this + 1);
  }

  KOKKOS_INLINE_FUNCTION
  bool requested_respawn() {
    // This should only be called when a task has finished executing and is
    // in the transition to either the complete or executing-respawn state.
    TaskBase* const lock = reinterpret_cast<TaskBase*>(LockTag);
    return lock != m_next;
  }

  KOKKOS_INLINE_FUNCTION
  void add_dependence(TaskBase* dep) {
    // Precondition: lock == m_next

    TaskBase* const lock = (TaskBase*)LockTag;

    // Assign dependence to m_next.  It will be processed in the subsequent
    // call to schedule.  Error if the dependence is reset.
    if (lock != Kokkos::atomic_exchange(&m_next, dep)) {
      Kokkos::abort("TaskScheduler ERROR: resetting task dependence");
    }

    if (nullptr != dep) {
      // The future may be destroyed upon returning from this call
      // so increment reference count to track this assignment.
      Kokkos::atomic_increment(&(dep->m_ref_count));
    }
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  int32_t reference_count() const {
    return *((int32_t volatile*)(&m_ref_count));
  }
};

//------------------------------------------------------------------------------
// <editor-fold desc="Verify the size of TaskBase is as expected"> {{{2

// Workaround: some compilers implement int16_t as 4 bytes, so the size might
// not actually be 48 bytes.
// There's not a lot of reason to keep checking this here; the program will
// work fine if this isn't true. I think this check was originally here to
// emphasize the fact that adding to the size of TaskBase could have a
// significant performance penalty, since doing so could substantially decrease
// the number of full task types that fit into a cache line.  We'll leave it
// here for now, though, since we're probably going to be ripping all of the
// old TaskBase stuff out eventually anyway.
constexpr size_t unpadded_task_base_size = 44 + 2 * sizeof(int16_t);
// don't forget padding:
constexpr size_t task_base_misalignment =
    unpadded_task_base_size % alignof(void*);
constexpr size_t task_base_padding_size =
    (alignof(void*) - task_base_misalignment) % alignof(void*);
constexpr size_t expected_task_base_size =
    unpadded_task_base_size + task_base_padding_size;

// Produce a more readable compiler error message than the plain static assert
template <size_t Size>
struct verify_task_base_size_is_48_note_actual_size_is_ {};
template <>
struct verify_task_base_size_is_48_note_actual_size_is_<
    expected_task_base_size> {
  using type = int;
};
static constexpr
    typename verify_task_base_size_is_48_note_actual_size_is_<sizeof(
        TaskBase)>::type verify = {};

static_assert(sizeof(TaskBase) == expected_task_base_size,
              "Verifying expected sizeof(TaskBase)");

// </editor-fold> end Verify the size of TaskBase is as expected }}}2
//------------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class Scheduler, typename ResultType, class FunctorType>
class Task : public TaskBase, public FunctorType {
 public:
  Task()            = delete;
  Task(Task&&)      = delete;
  Task(const Task&) = delete;
  Task& operator=(Task&&) = delete;
  Task& operator=(const Task&) = delete;

  using root_type    = TaskBase;
  using functor_type = FunctorType;
  using result_type  = ResultType;

  using specialization = TaskQueueSpecialization<Scheduler>;
  using member_type    = typename specialization::member_type;

  KOKKOS_INLINE_FUNCTION
  void apply_functor(member_type* const member, void*) {
    this->functor_type::operator()(*member);
  }

  template <typename T>
  KOKKOS_INLINE_FUNCTION void apply_functor(member_type* const member,
                                            T* const result) {
    this->functor_type::operator()(*member, *result);
  }

  KOKKOS_FUNCTION static void destroy(root_type* root) {
    TaskResult<result_type>::destroy(root);
  }

  KOKKOS_FUNCTION static void apply(root_type* root, void* exec) {
    Task* const task          = static_cast<Task*>(root);
    member_type* const member = reinterpret_cast<member_type*>(exec);
    result_type* const result = TaskResult<result_type>::ptr(task);

    // Task may be serial or team.
    // If team then must synchronize before querying if respawn was requested.
    // If team then only one thread calls destructor.

    const bool only_one_thread =
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
        0 == threadIdx.x && 0 == threadIdx.y;
#else
        0 == member->team_rank();
#endif

    task->apply_functor(member, result);

    member->team_barrier();

    if (only_one_thread && !(task->requested_respawn())) {
      // Did not respawn, destroy the functor to free memory.
      task->functor_type::~functor_type();
      // Cannot destroy and deallocate the task until its dependences
      // have been processed.
    }
  }

  // Constructor for runnable task
  KOKKOS_INLINE_FUNCTION constexpr Task(FunctorType&& arg_functor)
      : root_type(), functor_type(std::move(arg_functor)) {}

  KOKKOS_INLINE_FUNCTION
  ~Task() = delete;
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_TASKBASE_HPP */
