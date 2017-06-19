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

#ifndef KOKKOS_QTHREADS_TASKQUEUE_HPP
#define KOKKOS_QTHREADS_TASKQUEUE_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_QTHREADS ) && defined( KOKKOS_ENABLE_TASKPOLICY )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Manage task allocation, deallocation, and scheduling.
 *
 *  Task execution is handled here directly for the Qthread implementation.
 */
template<>
class TaskQueue< Kokkos::Qthread > {
private:

  using execution_space = Kokkos::Qthread ;
  using memory_space    = Kokkos::HostSpace
  using device_type     = Kokkos::Device< execution_space, memory_space > ;
  using memory_pool     = Kokkos::MemoryPool< device_type > ;
  using task_root_type  = Kokkos::Impl::TaskBase< execution_space, void, void > ;

  friend class Kokkos::TaskScheduler< execution_space > ;

  struct Destroy {
    TaskQueue * m_queue ;
    void destroy_shared_allocation();
  };

  //----------------------------------------

  enum : int { TASK_STATE_NULL         =  0,  ///<  Does not exist
               TASK_STATE_CONSTRUCTING =  1,  ///<  Is under construction
               TASK_STATE_WAITING      =  2,  ///<  Is waiting for execution
               TASK_STATE_EXECUTING    =  4,  ///<  Is executing
               TASK_STATE_RESPAWN      =  8,  ///<  Requested respawn
               TASK_STATE_COMPLETE     = 16   ///<  Execution is complete
             };

  // Queue is organized as [ priority ][ type ]

  memory_pool  m_memory ;
  unsigned     m_team_size ;   // Number of threads in a team
  long         m_accum_alloc ; // Accumulated number of allocations
  int          m_count_alloc ; // Current number of allocations
  int          m_max_alloc ;   // Maximum number of allocations
  int          m_ready_count ; // Number of ready or executing

  //----------------------------------------

  ~TaskQueue();
  TaskQueue() = delete ;
  TaskQueue( TaskQueue && ) = delete ;
  TaskQueue( TaskQueue const & ) = delete ;
  TaskQueue & operator = ( TaskQueue && ) = delete ;
  TaskQueue & operator = ( TaskQueue const & ) = delete ;

  TaskQueue
    ( const memory_space & arg_space,
      unsigned const arg_memory_pool_capacity,
      unsigned const arg_memory_pool_superblock_capacity_log2
    );

  // Schedule a task
  //   Precondition:
  //     task is not executing
  //     task->m_next is the dependence or zero
  //   Postcondition:
  //     task->m_next is linked list membership
  KOKKOS_FUNCTION
  void schedule( task_root_type * const );

  // Reschedule a task
  //   Precondition:
  //     task is in Executing state
  //     task->m_next == LockTag
  //   Postcondition:
  //     task is in Executing-Respawn state
  //     task->m_next == 0 (no dependence)
  KOKKOS_FUNCTION
  void reschedule( task_root_type * );

  // Complete a task
  //   Precondition:
  //     task is not executing
  //     task->m_next == LockTag  =>  task is complete
  //     task->m_next != LockTag  =>  task is respawn
  //   Postcondition:
  //     task->m_wait == LockTag  =>  task is complete
  //     task->m_wait != LockTag  =>  task is waiting
  KOKKOS_FUNCTION
  void complete( task_root_type * );

public:

  // If and only if the execution space is a single thread
  // then execute ready tasks.
  KOKKOS_INLINE_FUNCTION
  void iff_single_thread_recursive_execute()
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      specialization::iff_single_thread_recursive_execute( this );
#endif
    }

  void execute() { specialization::execute( this ); }

  template< typename FunctorType >
  void proc_set_apply( typename task_root_type::function_type * ptr )
    {
      specialization::template proc_set_apply< FunctorType >( ptr );
    }

  // Assign task pointer with reference counting of assigned tasks
  template< typename LV, typename RV >
  KOKKOS_FUNCTION static
  void assign( TaskBase< execution_space, LV, void > ** const lhs,
               TaskBase< execution_space, RV, void > *  const rhs )
    {
      using task_lhs = TaskBase< execution_space, LV, void > ;
#if 0
  {
    printf( "assign( 0x%lx { 0x%lx %d %d }, 0x%lx { 0x%lx %d %d } )\n",
            uintptr_t( lhs ? *lhs : 0 ),
            uintptr_t( lhs && *lhs ? (*lhs)->m_next : 0 ),
            int( lhs && *lhs ? (*lhs)->m_task_type : 0 ),
            int( lhs && *lhs ? (*lhs)->m_ref_count : 0 ),
            uintptr_t(rhs),
            uintptr_t( rhs ? rhs->m_next : 0 ),
            int( rhs ? rhs->m_task_type : 0 ),
            int( rhs ? rhs->m_ref_count : 0 )
          );
    fflush( stdout );
  }
#endif

      if ( *lhs )
      {
        const int count = Kokkos::atomic_fetch_add( &((*lhs)->m_ref_count), -1 );

        if ( ( 1 == count ) && ( (*lhs)->m_state == TASK_STATE_COMPLETE ) ) {
          // Reference count is zero and task is complete, deallocate.
          (*lhs)->m_queue->deallocate( *lhs, (*lhs)->m_alloc_size );
        }
        else if ( count <= 1 ) {
          Kokkos::abort("TaskScheduler task has negative reference count or is incomplete" );
        }

        // GEM: Should I check that there are no dependences here?  Can the state
        //      be set to complete while there are still dependences?
      }

      if ( rhs ) { Kokkos::atomic_fetch_add( &(rhs->m_ref_count), 1 ); }

      // Force write of *lhs

      *static_cast< task_lhs * volatile * >(lhs) = rhs ;

      Kokkos::memory_fence();
    }

  KOKKOS_FUNCTION
  size_t allocate_block_size( size_t n ); ///< Actual block size allocated

  KOKKOS_FUNCTION
  void * allocate( size_t n ); ///< Allocate from the memory pool

  KOKKOS_FUNCTION
  void deallocate( void * p, size_t n ); ///< Deallocate to the memory pool
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
class TaskBase< Kokkos::Qthread, void, void >
{
public:

  enum : int16_t   { TaskTeam   = TaskBase< void, void, void >::TaskTeam,
                     TaskSingle = TaskBase< void, void, void >::TaskSingle,
                     Aggregate  = TaskBase< void, void, void >::Aggregate };

  enum : uintptr_t { LockTag = TaskBase< void, void, void >::LockTag,
                     EndTag  = TaskBase< void, void, void >::EndTag };

  using execution_space = Kokkos::Qthread ;
  using queue_type      = TaskQueue< execution_space > ;

  template< typename > friend class Kokkos::TaskScheduler ;

  typedef void (* function_type) ( TaskBase *, void * );

  // sizeof(TaskBase) == 48

  function_type  m_apply ;       ///< Apply function pointer
  queue_type   * m_queue ;       ///< Queue in which this task resides
  TaskBase     * m_dep ;         ///< Dependence
  int32_t        m_ref_count ;   ///< Reference count
  int32_t        m_alloc_size ;  ///< Allocation size
  int32_t        m_dep_count ;   ///< Aggregate's number of dependences
  int16_t        m_task_type ;   ///< Type of task
  int16_t        m_priority ;    ///< Priority of runnable task
  aligned_t      m_qfeb ;        ///< Qthread full/empty bit
  int            m_state ;       ///< State of the task

  TaskBase( TaskBase && ) = delete ;
  TaskBase( const TaskBase & ) = delete ;
  TaskBase & operator = ( TaskBase && ) = delete ;
  TaskBase & operator = ( const TaskBase & ) = delete ;

  KOKKOS_INLINE_FUNCTION ~TaskBase() = default ;

  KOKKOS_INLINE_FUNCTION
  constexpr TaskBase() noexcept
    : m_apply(0),
      m_queue(0),
      m_dep(0),
      m_ref_count(0),
      m_alloc_size(0),
      m_dep_count(0),
      m_task_type( TaskSingle ),
      m_priority( 1 /* TaskRegularPriority */ ),
      m_qfeb(0),
      m_state( queue_type::TASK_STATE_CONSTRUCTING )
    {
      qthread_empty( & m_qfeb ); // Set to full when complete
    }

  //----------------------------------------

  static aligned_t qthread_func( void * arg );

  KOKKOS_INLINE_FUNCTION
  TaskBase ** aggregate_dependences()
    { return reinterpret_cast<TaskBase**>( this + 1 ); }

  KOKKOS_INLINE_FUNCTION
  void requested_respawn()
    { return m_state == queue_type::TASK_STATE_RESPAWN; }

  KOKKOS_INLINE_FUNCTION
  void add_dependence( TaskBase* dep )
    {
      // Assign dependence to m_dep.  It will be processed in the subsequent
      // call to schedule.  Error if the dependence is reset.
      if ( 0 != Kokkos::atomic_exchange( & m_dep, dep ) ) {
        Kokkos::abort("TaskScheduler ERROR: resetting task dependence");
      }

      if ( 0 != dep ) {
        // The future may be destroyed upon returning from this call
        // so increment reference count to track this assignment.
        Kokkos::atomic_fetch_add( &(dep->m_ref_count), 1 );
      }
    }

  using get_return_type = void ;

  KOKKOS_INLINE_FUNCTION
  get_return_type get() const {}
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKPOLICY ) */
#endif // KOKKOS_QTHREADS_TASKQUEUE_HPP

