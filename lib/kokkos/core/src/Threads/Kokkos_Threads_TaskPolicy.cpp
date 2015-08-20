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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

// Experimental unified task-data parallel manycore LDRD

#include <stdio.h>
#include <iostream>
#include <sstream>
#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

#if defined( KOKKOS_HAVE_PTHREAD )

namespace Kokkos {
namespace Experimental {
namespace Impl {

typedef TaskMember< Kokkos::Threads , void , void > Task ;

namespace {

int    volatile s_count_serial = 0 ;
int    volatile s_count_team   = 0 ;
Task * volatile s_ready_team   = 0 ;
Task * volatile s_ready_serial = 0 ;
Task * const    s_lock   = reinterpret_cast<Task*>( ~((unsigned long)0) );
Task * const    s_denied = reinterpret_cast<Task*>( ~((unsigned long)0) - 1 );

} /* namespace */
} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

namespace Kokkos {
namespace Experimental {

TaskPolicy< Kokkos::Threads >::TaskPolicy
  ( const unsigned arg_default_dependence_capacity
  , const unsigned arg_team_size
  )
  : m_default_dependence_capacity( arg_default_dependence_capacity )
  , m_team_size( arg_team_size )
{
  const int threads_total    = Threads::thread_pool_size(0);
  const int threads_per_numa = Threads::thread_pool_size(1);
  const int threads_per_core = Threads::thread_pool_size(2);

  if ( 0 == arg_team_size ) {
    // If a team task then claim for execution until count is zero
    // Issue: team collectives cannot assume which pool members are in the team.
    // Issue: team must only span a single NUMA region.

    // If more than one thread per core then map cores to work team,
    // else  map numa to work team.

    if      ( 1 < threads_per_core ) m_team_size = threads_per_core ;
    else if ( 1 < threads_per_numa ) m_team_size = threads_per_numa ;
    else                             m_team_size = 1 ;
  }

  // Verify a valid team size
  const bool valid_team_size =
    ( 0 < m_team_size && m_team_size <= threads_total ) &&
    (
      ( 1                == m_team_size ) ||
      ( threads_per_core == m_team_size ) ||
      ( threads_per_numa == m_team_size )
    );

  if ( ! valid_team_size ) {
    std::ostringstream msg ;

    msg << "Kokkos::Experimental::TaskPolicy< Kokkos::Threads > ERROR"
        << " invalid team_size(" << m_team_size << ")"
        << " threads_per_core(" << threads_per_core << ")"
        << " threads_per_numa(" << threads_per_numa << ")"
        << " threads_total(" << threads_total << ")"
        ;

    Kokkos::Impl::throw_runtime_exception( msg.str() );

  }
}

TaskPolicy< Kokkos::Threads >::member_type &
TaskPolicy< Kokkos::Threads >::member_single()
{
  static member_type s ;
  return s ;
}

void wait( Kokkos::Experimental::TaskPolicy< Kokkos::Threads > & policy )
{
  typedef Kokkos::Impl::ThreadsExecTeamMember member_type ;

  enum { BASE_SHMEM = 1024 };

  void * const arg = reinterpret_cast<void*>( long( policy.m_team_size ) );

  Kokkos::Impl::ThreadsExec::resize_scratch( 0 , member_type::team_reduce_size() + BASE_SHMEM );
  Kokkos::Impl::ThreadsExec::start( & Impl::Task::execute_ready_tasks_driver , arg );
  Kokkos::Impl::ThreadsExec::fence();
}

} /* namespace Experimental */
} /* namespace Kokkos */

namespace Kokkos {
namespace Experimental {
namespace Impl {

//----------------------------------------------------------------------------

void Task::throw_error_verify_type()
{
  Kokkos::Impl::throw_runtime_exception("TaskMember< Threads >::verify_type ERROR");
}

void Task::deallocate( void * ptr )
{
  free( ptr );
}

void * Task::allocate( const unsigned n )
{
  void * const ptr = malloc(n);

  return ptr ;
}

Task::~TaskMember()
{
}

//----------------------------------------------------------------------------

void Task::reschedule()
{
  // Reschedule transitions from executing back to waiting.
  const int old_state = atomic_compare_exchange( & m_state , int(TASK_STATE_EXECUTING) , int(TASK_STATE_WAITING) );

  if ( old_state != int(TASK_STATE_EXECUTING) ) {

fprintf( stderr
       , "reschedule ERROR task[%lx] state(%d)\n"
       , (unsigned long) this
       , old_state
       );
fflush(stderr);

  }
}

void Task::schedule()
{
  //----------------------------------------
  // State is either constructing or already waiting.
  // If constructing then transition to waiting.

  {
    const int old_state = atomic_compare_exchange( & m_state , int(TASK_STATE_CONSTRUCTING) , int(TASK_STATE_WAITING) );
    Task * const waitTask = *((Task * volatile const *) & m_wait );
    Task * const next = *((Task * volatile const *) & m_next );

    if ( s_denied == waitTask || 0 != next ||
         ( old_state != int(TASK_STATE_CONSTRUCTING) &&
           old_state != int(TASK_STATE_WAITING) ) ) {
      fprintf(stderr,"Task::schedule task(0x%lx) STATE ERROR: state(%d) wait(0x%lx) next(0x%lx)\n"
                    , (unsigned long) this
                    , old_state
                    , (unsigned long) waitTask
                    , (unsigned long) next );
      fflush(stderr);
      Kokkos::Impl::throw_runtime_exception("Kokkos::Impl::Task spawn or respawn state error");
    }
  }

  //----------------------------------------
  // Insert this task into another dependence that is not complete
  // Push on to the wait queue, fails if ( s_denied == m_dep[i]->m_wait )

  bool insert_in_ready_queue = true ;

  for ( int i = 0 ; i < m_dep_size && insert_in_ready_queue ; ) {

    Task * const task_dep = m_dep[i] ;
    Task * const head_value_old = *((Task * volatile *) & task_dep->m_wait );

    if ( s_denied == head_value_old ) {
      // Wait queue is closed, try again with the next queue
      ++i ;
    }
    else {

      // Wait queue is open and not locked.
      // If CAS succeeds then have acquired the lock.

      // Have exclusive access to this task.
      // Assign m_next assuming a successfull insertion into the queue.
      // Fence the memory assignment before attempting the CAS.

      *((Task * volatile *) & m_next ) = head_value_old ;

      memory_fence();

      // Attempt to insert this task into the queue

      Task * const wait_queue_head = atomic_compare_exchange( & task_dep->m_wait , head_value_old , this );

      if ( head_value_old == wait_queue_head ) {
        insert_in_ready_queue = false ;
      }
    }
  }

  //----------------------------------------
  // All dependences are complete, insert into the ready list

  if ( insert_in_ready_queue ) {

    // Increment the count of ready tasks.
    // Count is decremented when task is complete.

    Task * volatile * queue = 0 ;

    if ( m_serial ) {
      atomic_increment( & s_count_serial );
      queue = & s_ready_serial ;
    }
    else {
      atomic_increment( & s_count_team );
      queue = & s_ready_team ;
    }

    while ( insert_in_ready_queue ) {

      Task * const head_value_old = *queue ;

      if ( s_lock != head_value_old ) {
        // Read the head of ready queue, if same as previous value then CAS locks the ready queue
        // Only access via CAS

        // Have exclusive access to this task, assign to head of queue, assuming successful insert
        // Fence assignment before attempting insert.
        *((Task * volatile *) & m_next ) = head_value_old ;

        memory_fence();

        Task * const ready_queue_head = atomic_compare_exchange( queue , head_value_old , this );

        if ( head_value_old == ready_queue_head ) {
          // Successful insert
          insert_in_ready_queue = false ; // done
        }
      }
    }
  }
}

//----------------------------------------------------------------------------

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )

void Task::assign( Task ** const lhs_ptr , Task * rhs )
{
  // Increment rhs reference count.
  if ( rhs ) { atomic_increment( & rhs->m_ref_count ); }

  // Assign the pointer and retrieve the previous value.

  Task * const old_lhs = atomic_exchange( lhs_ptr , rhs );

  if ( old_lhs ) {

    // Decrement former lhs reference count.
    // If reference count is zero task must be complete, then delete task.
    // Task is ready for deletion when  wait == s_denied

    int const count = atomic_fetch_add( & (old_lhs->m_ref_count) , -1 ) - 1 ;

    // if 'count != 0' then 'old_lhs' may be deallocated before dereferencing
    Task * const wait = count == 0 ? *((Task * const volatile *) & old_lhs->m_wait ) : (Task*) 0 ;

    if ( count < 0 || ( count == 0 && wait != s_denied ) ) {

      static const char msg_error_header[]  = "Kokkos::Impl::TaskManager<Kokkos::Threads>::assign ERROR deleting" ;

      fprintf( stderr , "%s task(0x%lx) m_ref_count(%d) , m_wait(0x%ld)\n"
                      , msg_error_header
                      , (unsigned long) old_lhs
                      , count
                      , (unsigned long) wait );
      fflush(stderr);

      Kokkos::Impl::throw_runtime_exception( msg_error_header );
    }

    if ( count == 0 ) {
      // When 'count == 0' this thread has exclusive access to 'old_lhs'
      const Task::function_dealloc_type d = old_lhs->m_dealloc ;
      (*d)( old_lhs );
    }
  }
}

#endif

//----------------------------------------------------------------------------

Task * Task::get_dependence( int i ) const
{
  Task * const t = m_dep[i] ;

  if ( Kokkos::Experimental::TASK_STATE_EXECUTING != m_state || i < 0 || m_dep_size <= i || 0 == t ) {

fprintf( stderr
       , "TaskMember< Threads >::get_dependence ERROR : task[%lx]{ state(%d) dep_size(%d) dep[%d] = %lx }\n"
       , (unsigned long) this
       , m_state
       , m_dep_size
       , i
       , (unsigned long) t
       );
fflush( stderr );

    Kokkos::Impl::throw_runtime_exception("TaskMember< Threads >::get_dependence ERROR");
  }

  return t ;
}

//----------------------------------------------------------------------------

void Task::add_dependence( Task * before )
{
  if ( before != 0 ) {

    int const state = *((volatile const int *) & m_state );

    // Can add dependence during construction or during execution

    if ( ( Kokkos::Experimental::TASK_STATE_CONSTRUCTING == state ||
           Kokkos::Experimental::TASK_STATE_EXECUTING    == state ) &&
         m_dep_size < m_dep_capacity ) {

      ++m_dep_size ;

      assign( m_dep + (m_dep_size-1) , before );

      memory_fence();
    }
    else {

fprintf( stderr
       , "TaskMember< Threads >::add_dependence ERROR : task[%lx]{ state(%d) dep_size(%d) m_dep_capacity(%d) }\n"
       , (unsigned long) this
       , m_state
       , m_dep_size
       , m_dep_capacity
       );
fflush( stderr );

      Kokkos::Impl::throw_runtime_exception("TaskMember< Threads >::add_dependence ERROR");
    }
  }
}

//----------------------------------------------------------------------------

void Task::clear_dependence()
{
  for ( int i = m_dep_size - 1 ; 0 <= i ; --i ) {
    assign( m_dep + i , 0 );
  }

  *((volatile int *) & m_dep_size ) = 0 ;

  memory_fence();
}

//----------------------------------------------------------------------------

Task * Task::pop_ready_task( Task * volatile * const queue )
{
  Task * const task_old = *queue ;

  if ( s_lock != task_old && 0 != task_old ) {

    Task * const task = atomic_compare_exchange( queue , task_old , s_lock );

    if ( task_old == task ) {

      // May have acquired the lock and task.
      // One or more other threads may have acquired this same task and lock
      // due to respawning ABA race condition.
      // Can only be sure of acquire with a successful state transition from waiting to executing

      const int old_state = atomic_compare_exchange( & task->m_state, int(TASK_STATE_WAITING), int(TASK_STATE_EXECUTING) );

      if ( old_state == int(TASK_STATE_WAITING) ) {

        // Transitioned this task from waiting to executing
        // Update the queue to the next entry and release the lock

        Task * const next_old = *((Task * volatile *) & task->m_next );

        Task * const s = atomic_compare_exchange( queue , s_lock , next_old );

        if ( s != s_lock ) {
          fprintf(stderr,"Task::pop_ready_task( 0x%lx ) UNLOCK ERROR\n", (unsigned long) queue );
          fflush(stderr);
        }

        *((Task * volatile *) & task->m_next ) = 0 ;

        return task ;
      }
      else {
        fprintf(stderr,"Task::pop_ready_task( 0x%lx ) task(0x%lx) state(%d) ERROR\n"
                      , (unsigned long) queue
                      , (unsigned long) task
                      , old_state );
        fflush(stderr);
      }
    }
  }

  return (Task *) 0 ;
}


void Task::complete_executed_task( Task * task , volatile int * const queue_count )
{
  // State is either executing or if respawned then waiting,
  // try to transition from executing to complete.
  // Reads the current value.

  const int state_old =
    atomic_compare_exchange( & task->m_state
                           , int(Kokkos::Experimental::TASK_STATE_EXECUTING)
                           , int(Kokkos::Experimental::TASK_STATE_COMPLETE) );

  if ( Kokkos::Experimental::TASK_STATE_WAITING == state_old ) {
    task->schedule(); /* Task requested a respawn so reschedule it */
  }
  else if ( Kokkos::Experimental::TASK_STATE_EXECUTING != state_old ) {
    fprintf( stderr
           , "TaskMember< Threads >::execute_serial completion ERROR : task[%lx]{ state_old(%d) dep_size(%d) }\n"
           , (unsigned long) & task
           , state_old
           , task->m_dep_size
           );
    fflush( stderr );
  }
  else {

    // Clear dependences of this task before locking wait queue

    task->clear_dependence();

    // Stop other tasks from adding themselves to this task's wait queue.
    // The wait queue is updated concurrently so guard with an atomic.
    // Setting the wait queue to denied denotes delete-ability of the task by any thread.
    // Therefore, once 'denied' the task pointer must be treated as invalid.

    Task * wait_queue     = *((Task * volatile *) & task->m_wait );
    Task * wait_queue_old = 0 ;

    do {
      wait_queue_old = wait_queue ;
      wait_queue     = atomic_compare_exchange( & task->m_wait , wait_queue_old , s_denied );
    } while ( wait_queue_old != wait_queue );

    task = 0 ;

    // Pop waiting tasks and schedule them
    while ( wait_queue ) {
      Task * const x = wait_queue ; wait_queue = x->m_next ; x->m_next = 0 ;
      x->schedule();
    }
  }

  atomic_decrement( queue_count );
}

//----------------------------------------------------------------------------

void Task::execute_ready_tasks_driver( Kokkos::Impl::ThreadsExec & exec , const void * arg )
{
  typedef Kokkos::Impl::ThreadsExecTeamMember member_type ;

  // Whole pool is calling this function

  // Create the thread team member with shared memory for the given task.
  const int team_size = reinterpret_cast<long>( arg );

  member_type member( & exec , TeamPolicy< Kokkos::Threads >( 1 , team_size ) , 0 );

  Kokkos::Impl::ThreadsExec & exec_team_base = member.threads_exec_team_base();

  Task * volatile * const task_team_ptr = reinterpret_cast<Task**>( exec_team_base.reduce_memory() );

  if ( member.team_fan_in() ) {
    *task_team_ptr = 0 ;
    Kokkos::memory_fence();
  }
  member.team_fan_out();

  long int iteration_count = 0 ;

  // Each team must iterate this loop synchronously to insure team-execution of team-task

  while ( 0 < s_count_serial || 0 < s_count_team ) {

    if ( member.team_rank() == 0 ) {
      // Only one team member attempts to pop a team task
      *task_team_ptr = pop_ready_task( & s_ready_team );
    }

    // Query if team acquired a team task
    Task * const task_team = *task_team_ptr ;

    if ( task_team ) {
      // Set shared memory
      member.set_league_shmem( 0 , 1 , task_team->m_shmem_size );

      (*task_team->m_team)( task_team , member );

      // Do not proceed until all members have completed the task,
      // the task has been completed or rescheduled, and
      // the team task pointer has been cleared.
      if ( member.team_fan_in() ) {
        complete_executed_task( task_team , & s_count_team );
        *task_team_ptr = 0 ;
        Kokkos::memory_fence();
      }
      member.team_fan_out();
    }
    else {
      Task * const task_serial = pop_ready_task( & s_ready_serial );

      if ( task_serial ) {
        if ( task_serial->m_serial ) (*task_serial->m_serial)( task_serial );

        complete_executed_task( task_serial , & s_count_serial );
      }
    }

    ++iteration_count ;
  }

  exec.fan_in();
}

} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_HAVE_PTHREAD ) */

