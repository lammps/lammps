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

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_SERIAL ) && defined( KOKKOS_ENABLE_TASKDAG )

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_Serial_Task.hpp>
#include <impl/Kokkos_TaskQueue_impl.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template class TaskQueue< Kokkos::Serial > ;

void TaskQueueSpecialization< Kokkos::Serial >::execute
  ( TaskQueue< Kokkos::Serial > * const queue )
{
  using exec_space = Kokkos::Serial ;
  using tqs_queue_type      = TaskQueue< exec_space > ;
  using task_root_type  = TaskBase< void , void , void > ;
  using Member          = Impl::HostThreadTeamMember< exec_space > ;

  task_root_type * const end = (task_root_type *) task_root_type::EndTag ;

  // Set default buffers
  serial_resize_thread_team_data( 0   /* global reduce buffer */
                                , 512 /* team reduce buffer */
                                , 0   /* team shared buffer */
                                , 0   /* thread local buffer */
                                );

  Impl::HostThreadTeamData * const data = Impl::serial_get_thread_team_data();

  Member exec( *data );

  // Loop until all queues are empty
  while ( 0 < queue->m_ready_count ) {

    task_root_type * task = end ;

    for ( int i = 0 ; i < tqs_queue_type::NumQueue && end == task ; ++i ) {
      for ( int j = 0 ; j < 2 && end == task ; ++j ) {
        task = tqs_queue_type::pop_ready_task( & queue->m_ready[i][j] );
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

void TaskQueueSpecialization< Kokkos::Serial > ::
  iff_single_thread_recursive_execute(
    TaskQueue< Kokkos::Serial > * const queue )
{
  using exec_space = Kokkos::Serial ;
  using tqs_queue_type      = TaskQueue< exec_space > ;
  using task_root_type  = TaskBase< void , void , void > ;
  using Member          = Impl::HostThreadTeamMember< exec_space > ;

  task_root_type * const end = (task_root_type *) task_root_type::EndTag ;

  Impl::HostThreadTeamData * const data = Impl::serial_get_thread_team_data();

  Member exec( *data );

  // Loop until no runnable task

  task_root_type * task = end ;

  do {

    task = end ;

    for ( int i = 0 ; i < tqs_queue_type::NumQueue && end == task ; ++i ) {
      for ( int j = 0 ; j < 2 && end == task ; ++j ) {
        task = tqs_queue_type::pop_ready_task( & queue->m_ready[i][j] );
      }
    }

    if ( end == task ) break ;

    (*task->m_apply)( task , & exec );

    queue->complete( task );

  } while(1);
}

}} /* namespace Kokkos::Impl */

#else
void KOKKOS_CORE_SRC_IMPL_SERIAL_TASK_PREVENT_LINK_ERROR() {}
#endif /* #if defined( KOKKOS_ENABLE_SERIAL ) && defined( KOKKOS_ENABLE_TASKDAG ) */

