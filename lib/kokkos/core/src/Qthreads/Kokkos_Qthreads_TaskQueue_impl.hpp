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

#ifndef KOKKOS_QTHREADS_TASKQUEUE_IMPL_HPP
#define KOKKOS_QTHREADS_TASKQUEUE_IMPL_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_QTHREADS ) && defined( KOKKOS_ENABLE_TASKPOLICY )

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ExecSpace >
void TaskQueue< ExecSpace >::Destroy::destroy_shared_allocation()
{
  m_queue->~TaskQueue();
}

//----------------------------------------------------------------------------

template< typename ExecSpace >
TaskQueue< ExecSpace >::TaskQueue
  ( const TaskQueue< ExecSpace >::memory_space & arg_space,
    unsigned const arg_memory_pool_capacity,
    unsigned const arg_memory_pool_superblock_capacity_log2 )
  : m_memory( arg_space,
              arg_memory_pool_capacity,
              arg_memory_pool_superblock_capacity_log2 )
    m_team_size( unsigned( qthread_num_workers_local(NO_SHEPHERD) ) ),
    m_accum_alloc(0),
    m_count_alloc(0),
    m_max_alloc(0),
    m_ready_count(0)
{}

//----------------------------------------------------------------------------

template< typename ExecSpace >
TaskQueue< ExecSpace >::~TaskQueue()
{
  // Verify that ready count is zero.
  if ( 0 != m_ready_count ) {
    Kokkos::abort("TaskQueue::~TaskQueue ERROR: has ready or executing tasks");
  }
}

//----------------------------------------------------------------------------

template< typename ExecSpace >
KOKKOS_FUNCTION
size_t TaskQueue< ExecSpace >::allocate_block_size( size_t n )
{
  return m_memory.allocate_block_size( n );
}

//----------------------------------------------------------------------------

template< typename ExecSpace >
KOKKOS_FUNCTION
void * TaskQueue< ExecSpace >::allocate( size_t n )
{
  void * const p = m_memory.allocate(n);

  if ( p ) {
    Kokkos::atomic_increment( & m_accum_alloc );
    Kokkos::atomic_increment( & m_count_alloc );

    if ( m_max_alloc < m_count_alloc ) m_max_alloc = m_count_alloc ;
  }

  return p ;
}

//----------------------------------------------------------------------------

template< typename ExecSpace >
KOKKOS_FUNCTION
void TaskQueue< ExecSpace >::deallocate( void * p, size_t n )
{
  m_memory.deallocate( p, n );
  Kokkos::atomic_decrement( & m_count_alloc );
}

//----------------------------------------------------------------------------

template< typename ExecSpace >
KOKKOS_FUNCTION
void TaskQueue< ExecSpace >::schedule
  ( TaskQueue< ExecSpace >::task_root_type * const task )
{
#if 0
  printf( "schedule( 0x%lx { %d %d %d }\n",
          uintptr_t(task),
          task->m_task_type,
          task->m_priority,
          task->m_ref_count );
#endif

  // The task has been constructed and is waiting to be executed.
  task->m_state = TASK_STATE_WAITING ;

  if ( task->m_task_type != task_root_type::Aggregate ) {
    // Scheduling a single or team task.

    // Increment active task count before spawning.
    Kokkos::atomic_increment( m_ready_count );

    if ( task->m_dep == 0 ) {
      // Schedule a task with no dependences.

      if ( task_root_type::TaskTeam == task->m_task_type && m_team_size > 1 ) {
        // If more than one shepherd spawn on a shepherd other than this shepherd
        const int num_shepherd  = qthread_num_shepherds();
        const int this_shepherd = qthread_shep();
        int spawn_shepherd      = ( this_shepherd + 1 ) % num_shepherd ;

#if 0
        fprintf( stdout,
                 "worker(%d.%d) task 0x%.12lx spawning on shepherd(%d) clone(%d)\n",
                 qthread_shep(),
                 qthread_worker_local(NULL),
                 reinterpret_cast<unsigned long>(this),
                 spawn_shepherd,
                 m_team_size - 1
               );
        fflush(stdout);
#endif

        qthread_spawn_cloneable(
          & task_root_type::qthread_func,
          task,
          0,
          NULL,
          0, // no depenedences
          0, // dependences array
          spawn_shepherd,
          unsigned( QTHREAD_SPAWN_SIMPLE | QTHREAD_SPAWN_LOCAL_PRIORITY ),
          m_team_size - 1
        );
      }
      else {
        qthread_spawn(
          & task_root_type::qthread_func,
          task,
          0,
          NULL,
          0, // no depenedences
          0, // dependences array
          NO_SHEPHERD,
          QTHREAD_SPAWN_SIMPLE /* allows optimization for non-blocking task */
        );
      }
    }
    else if ( task->m_dep->m_task_type != task_root_type::Aggregate )
    // Malloc the precondition array to pass to qthread_spawn().  For
    // non-aggregate tasks, it is a single pointer since there are no
    // dependences.  Qthreads will eventually free this allocation so memory will
    // not be leaked. Is malloc thread-safe?  Should this call be guarded?  The
    // memory can't be allocated from the pool allocator because Qthreads frees
    // it using free().
    aligned_t ** qprecon = (aligned_t **) malloc( sizeof(aligned_t *) );

    *qprecon = reinterpret_cast<aligned_t *>( uintptr_t(m_dep_size) );

    if ( task->m_task_type == task_root_type::TaskTeam && m_team_size > 1) {
      // If more than one shepherd spawn on a shepherd other than this shepherd
      const int num_shepherd  = qthread_num_shepherds();
      const int this_shepherd = qthread_shep();
      int spawn_shepherd      = ( this_shepherd + 1 ) % num_shepherd ;

#if 0
  fprintf( stdout,
           "worker(%d.%d) task 0x%.12lx spawning on shepherd(%d) clone(%d)\n",
           qthread_shep(),
           qthread_worker_local(NULL),
           reinterpret_cast<unsigned long>(this),
           spawn_shepherd,
           m_team_size - 1
         );
  fflush(stdout);
#endif

      qthread_spawn_cloneable(
        & Task::qthread_func,
        this,
        0,
        NULL,
        m_dep_size,
        qprecon, /* dependences */
        spawn_shepherd,
        unsigned( QTHREAD_SPAWN_SIMPLE | QTHREAD_SPAWN_LOCAL_PRIORITY ),
        m_team_size - 1
      );
    }
    else {
      qthread_spawn(
        & Task::qthread_func, /* function */
        this,                 /* function argument */
        0,
        NULL,
        m_dep_size,
        qprecon, /* dependences */
        NO_SHEPHERD,
        QTHREAD_SPAWN_SIMPLE /* allows optimization for non-blocking task */
      );
    }
  }
  else {
    // GEM: How do I handle an aggregate (when_all) task?
  }
}

//----------------------------------------------------------------------------

template< typename ExecSpace >
KOKKOS_FUNCTION
void TaskQueue< ExecSpace >::reschedule( task_root_type * task )
{
  // Precondition:
  //   task is in Executing state
  //   task->m_next == LockTag
  //
  // Postcondition:
  //   task is in Executing-Respawn state
  //   task->m_next == 0 (no dependence)

  task_root_type * const zero = (task_root_type *) 0 ;
  task_root_type * const lock = (task_root_type *) task_root_type::LockTag ;

  if ( lock != Kokkos::atomic_exchange( & task->m_next, zero ) ) {
    Kokkos::abort("TaskScheduler::respawn ERROR: already respawned");
  }
}

//----------------------------------------------------------------------------

template< typename ExecSpace >
KOKKOS_FUNCTION
void TaskQueue< ExecSpace >::complete
  ( TaskQueue< ExecSpace >::task_root_type * task )
{
  // Complete a runnable task that has finished executing
  // or a when_all task when all of its dependeneces are complete.

  task_root_type * const zero = (task_root_type *) 0 ;
  task_root_type * const lock = (task_root_type *) task_root_type::LockTag ;
  task_root_type * const end  = (task_root_type *) task_root_type::EndTag ;

#if 0
  printf( "complete( 0x%lx { 0x%lx 0x%lx %d %d %d }\n",
          uintptr_t(task),
          uintptr_t(task->m_wait),
          uintptr_t(task->m_next),
          task->m_task_type,
          task->m_priority,
          task->m_ref_count
        );
  fflush( stdout );
#endif

  const bool runnable = task_root_type::Aggregate != task->m_task_type ;

  //----------------------------------------

  if ( runnable && lock != task->m_next ) {
    // Is a runnable task has finished executing and requested respawn.
    // Schedule the task for subsequent execution.

    schedule( task );
  }
  //----------------------------------------
  else {
    // Is either an aggregate or a runnable task that executed
    // and did not respawn.  Transition this task to complete.

    // If 'task' is an aggregate then any of the runnable tasks that
    // it depends upon may be attempting to complete this 'task'.
    // Must only transition a task once to complete status.
    // This is controled by atomically locking the wait queue.

    // Stop other tasks from adding themselves to this task's wait queue
    // by locking the head of this task's wait queue.

    task_root_type * x = Kokkos::atomic_exchange( & task->m_wait, lock );

    if ( x != (task_root_type *) lock ) {

      // This thread has transitioned this 'task' to complete.
      // 'task' is no longer in a queue and is not executing
      // so decrement the reference count from 'task's creation.
      // If no other references to this 'task' then it will be deleted.

      TaskQueue::assign( & task, zero );

      // This thread has exclusive access to the wait list so
      // the concurrency-safe pop_task function is not needed.
      // Schedule the tasks that have been waiting on the input 'task',
      // which may have been deleted.

      while ( x != end ) {

        // Set x->m_next = zero  <=  no dependence

        task_root_type * const next =
          (task_root_type *) Kokkos::atomic_exchange( & x->m_next, zero );

        schedule( x );

        x = next ;
      }
    }
  }

  if ( runnable ) {
    // A runnable task was popped from a ready queue and executed.
    // If respawned into a ready queue then the ready count was incremented
    // so decrement whether respawned or not.
    Kokkos::atomic_decrement( & m_ready_count );
  }
}

//----------------------------------------------------------------------------

template<>
aligned_t
TaskBase< Kokkos::Qthreads, void, void >::qthread_func( void * arg )
{
  using execution_space = Kokkos::Qthreads ;
  using task_root_type  = TaskBase< execution_space , void , void > ;
  using Member          = Kokkos::Impl::QthreadsTeamPolicyMember;

  task_root_type * const task = reinterpret_cast< task_root_type * >( arg );

  // First member of the team change state to executing.
  // Use compare-exchange to avoid race condition with a respawn.
  Kokkos::atomic_compare_exchange_strong( & task->m_state,
                                          queue_type::TASK_STATE_WAITING,
                                          queue_type::TASK_STATE_EXECUTING
                                        );

  if ( task_root_type::TaskTeam == task->m_task_type )
  {
    if ( 1 < task->m_queue->m_team_size ) {
      // Team task with team size of more than 1.
      Member::TaskTeam task_team_tag ;

      // Initialize team size and rank with shephered info
      Member member( task_team_tag );

      (*task->m_apply)( task , & member );

#if 0
      fprintf( stdout,
              "worker(%d.%d) task 0x%.12lx executed by member(%d:%d)\n",
              qthread_shep(),
              qthread_worker_local(NULL),
              reinterpret_cast<unsigned long>(task),
              member.team_rank(),
              member.team_size()
            );
      fflush(stdout);
#endif

      member.team_barrier();
      if ( member.team_rank() == 0 ) task->closeout();
      member.team_barrier();
    }
    else {
      // Team task with team size of 1.
      Member member ;
      (*task->m_apply)( task , & member );
      task->closeout();
    }
  }
  else {
    (*task->m_apply)( task );
    task->closeout();
  }

#if 0
fprintf( stdout
       , "worker(%d.%d) task 0x%.12lx return\n"
       , qthread_shep()
       , qthread_worker_local(NULL)
       , reinterpret_cast<unsigned long>(task)
       );
fflush(stdout);
#endif

  return 0 ;
}

} /* namespace Impl */
} /* namespace Kokkos */


#endif /* #if defined( KOKKOS_ENABLE_TASKPOLICY ) */
#endif // KOKKOS_QTHREADS_TASKQUEUE_IMPL_HPP

