/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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

#ifndef KOKKOS_IMPL_SERIAL_TASK_HPP
#define KOKKOS_IMPL_SERIAL_TASK_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_TASKDAG )

#include <Kokkos_TaskScheduler_fwd.hpp>

#include <impl/Kokkos_TaskQueue.hpp>
#include <Kokkos_Serial.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<class QueueType>
class TaskQueueSpecialization<
  SimpleTaskScheduler<Kokkos::Serial, QueueType>
>
{
public:

  // Note: Scheduler may be an incomplete type at class scope (but not inside
  // of the methods, obviously)

  using execution_space = Kokkos::Serial;
  using memory_space = Kokkos::HostSpace;
  using scheduler_type = SimpleTaskScheduler<Kokkos::Serial, QueueType>;
  using member_type = TaskTeamMemberAdapter<
    HostThreadTeamMember<Kokkos::Serial>, scheduler_type
  >;

  static
  void execute(scheduler_type const& scheduler)
  {
    using task_base_type = typename scheduler_type::task_base_type;

    // Set default buffers
    serial_resize_thread_team_data(
      0,   /* global reduce buffer */
      512, /* team reduce buffer */
      0,   /* team shared buffer */
      0    /* thread local buffer */
    );

    Impl::HostThreadTeamData& self = *Impl::serial_get_thread_team_data();

    auto& queue = scheduler.queue();
    auto team_scheduler = scheduler.get_team_scheduler(0);

    member_type member(scheduler, self);

    auto current_task = OptionalRef<task_base_type>(nullptr);

    while(not queue.is_done()) {

      // Each team lead attempts to acquire either a thread team task
      // or a single thread task for the team.

      // pop a task off
      current_task = queue.pop_ready_task(team_scheduler.team_scheduler_info());

      // run the task
      if(current_task) {
        current_task->as_runnable_task().run(member);
        // Respawns are handled in the complete function
        queue.complete(
          (*std::move(current_task)).as_runnable_task(),
          team_scheduler.team_scheduler_info()
        );
      }

    }

  }

  static constexpr uint32_t
  get_max_team_count(execution_space const&) noexcept
  {
    return 1;
  }

  template <typename TaskType>
  static void
  get_function_pointer(
    typename TaskType::function_type& ptr,
    typename TaskType::destroy_type& dtor
  )
  {
    ptr = TaskType::apply;
    dtor = TaskType::destroy;
  }
};

//----------------------------------------------------------------------------

template<class Scheduler>
class TaskQueueSpecializationConstrained<
  Scheduler,
  typename std::enable_if<
    std::is_same<typename Scheduler::execution_space, Kokkos::Serial>::value
  >::type
>
{
public:

  // Note: Scheduler may be an incomplete type at class scope (but not inside
  // of the methods, obviously)

  using execution_space = Kokkos::Serial;
  using memory_space = Kokkos::HostSpace;
  using scheduler_type = Scheduler;
  using member_type = TaskTeamMemberAdapter<
    HostThreadTeamMember<Kokkos::Serial>, scheduler_type
  >;

  static
  void iff_single_thread_recursive_execute(scheduler_type const& scheduler) {
    using task_base_type = TaskBase;
    using queue_type = typename scheduler_type::queue_type;

    task_base_type * const end = (task_base_type *) task_base_type::EndTag ;

    Impl::HostThreadTeamData * const data = Impl::serial_get_thread_team_data();

    member_type exec( scheduler, *data );

    // Loop until no runnable task

    task_base_type * task = end ;

    auto* const queue = scheduler.m_queue;

    do {

      task = end ;

      for ( int i = 0 ; i < queue_type::NumQueue && end == task ; ++i ) {
        for ( int j = 0 ; j < 2 && end == task ; ++j ) {
          task = queue_type::pop_ready_task( & queue->m_ready[i][j] );
        }
      }

      if ( end == task ) break ;

      (*task->m_apply)( task , & exec );

      queue->complete( task );

    } while(1);

  }

  static
  void execute(scheduler_type const& scheduler)
  {
    using task_base_type = TaskBase;
    using queue_type = typename scheduler_type::queue_type;

    task_base_type * const end = (task_base_type *) task_base_type::EndTag ;

    // Set default buffers
    serial_resize_thread_team_data(
      0,   /* global reduce buffer */
      512, /* team reduce buffer */
      0,   /* team shared buffer */
      0    /* thread local buffer */
    );

    auto* const queue = scheduler.m_queue;

    Impl::HostThreadTeamData * const data = Impl::serial_get_thread_team_data();

    member_type exec( scheduler, *data );

    // Loop until all queues are empty
    while ( 0 < queue->m_ready_count ) {

      task_base_type * task = end ;

      for ( int i = 0 ; i < queue_type::NumQueue && end == task ; ++i ) {
        for ( int j = 0 ; j < 2 && end == task ; ++j ) {
          task = queue_type::pop_ready_task( & queue->m_ready[i][j] );
        }
      }

      if ( end != task ) {

        // pop_ready_task resulted in lock == task->m_next
        // In the executing state

        (*task->m_apply)( task , & exec );

#if 0
        printf( "TaskQueue<Serial>::executed: 0x%lx { 0x%lx 0x%lx %d %d %d }\n"
        , uintptr_t(task)
        , uintptr_t(task->m_wait)
        , uintptr_t(task->m_next)
        , task->m_task_type
        , task->m_priority
        , task->m_ref_count );
#endif

        // If a respawn then re-enqueue otherwise the task is complete
        // and all tasks waiting on this task are updated.
        queue->complete( task );
      }
      else if ( 0 != queue->m_ready_count ) {
        Kokkos::abort("TaskQueue<Serial>::execute ERROR: ready_count");
      }
    }
  }

  template <typename TaskType>
  static void
  get_function_pointer(
    typename TaskType::function_type& ptr,
    typename TaskType::destroy_type& dtor
  )
  {
    ptr = TaskType::apply;
    dtor = TaskType::destroy;
  }
};

extern template class TaskQueue< Kokkos::Serial, typename Kokkos::Serial::memory_space > ;

}} /* namespace Kokkos::Impl */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_SERIAL_TASK_HPP */

