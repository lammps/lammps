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

#ifndef KOKKOS_IMPL_TASKQUEUE_HPP
#define KOKKOS_IMPL_TASKQUEUE_HPP

#if defined( KOKKOS_ENABLE_TASKDAG )

#include <string>
#include <typeinfo>
#include <stdexcept>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/*\brief  Implementation data for task data management, access, and execution.
 *
 *  Curiously recurring template pattern (CRTP)
 *  to allow static_cast from the
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
template< typename Space , typename ResultType , typename FunctorType >
class TaskBase ;

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename Space >
class TaskQueueSpecialization ;

/** \brief  Manage task allocation, deallocation, and scheduling.
 *
 *  Task execution is deferred to the TaskQueueSpecialization.
 *  All other aspects of task management have shared implementation.
 */
template< typename ExecSpace >
class TaskQueue {
private:

  friend class TaskQueueSpecialization< ExecSpace > ;
  friend class Kokkos::TaskScheduler< ExecSpace > ;

  using execution_space = ExecSpace ;
  using specialization  = TaskQueueSpecialization< execution_space > ;
  using memory_space    = typename specialization::memory_space ;
  using device_type     = Kokkos::Device< execution_space , memory_space > ;
  using memory_pool     = Kokkos::Experimental::MemoryPool< device_type > ;
  using task_root_type  = Kokkos::Impl::TaskBase<execution_space,void,void> ;

  struct Destroy {
    TaskQueue * m_queue ;
    void destroy_shared_allocation();
  };

  //----------------------------------------

  enum : int { NumQueue = 3 };

  // Queue is organized as [ priority ][ type ]

  memory_pool               m_memory ;
  task_root_type * volatile m_ready[ NumQueue ][ 2 ];
  long                      m_accum_alloc ; // Accumulated number of allocations
  int                       m_count_alloc ; // Current number of allocations
  int                       m_max_alloc ;   // Maximum number of allocations
  int                       m_ready_count ; // Number of ready or executing

  //----------------------------------------

  ~TaskQueue();
  TaskQueue() = delete ;
  TaskQueue( TaskQueue && ) = delete ;
  TaskQueue( TaskQueue const & ) = delete ;
  TaskQueue & operator = ( TaskQueue && ) = delete ;
  TaskQueue & operator = ( TaskQueue const & ) = delete ;

  TaskQueue
    ( const memory_space & arg_space
    , unsigned const arg_memory_pool_capacity
    , unsigned const arg_memory_pool_superblock_capacity_log2
    );

  // Schedule a task
  //   Precondition:
  //     task is not executing
  //     task->m_next is the dependence or zero
  //   Postcondition:
  //     task->m_next is linked list membership
  KOKKOS_FUNCTION void schedule_runnable(  task_root_type * const );
  KOKKOS_FUNCTION void schedule_aggregate( task_root_type * const );

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

  KOKKOS_FUNCTION
  static bool push_task( task_root_type * volatile * const
                       , task_root_type * const );

  KOKKOS_FUNCTION
  static task_root_type * pop_ready_task( task_root_type * volatile * const );

  KOKKOS_FUNCTION static
  void decrement( task_root_type * task );

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
  template< typename LV , typename RV >
  KOKKOS_FUNCTION static
  void assign( TaskBase< execution_space,LV,void> ** const lhs
             , TaskBase< execution_space,RV,void> *  const rhs )
    {
      using task_lhs = TaskBase< execution_space,LV,void> ;
#if 0
  {
    printf( "assign( 0x%lx { 0x%lx %d %d } , 0x%lx { 0x%lx %d %d } )\n"
          , uintptr_t( lhs ? *lhs : 0 )
          , uintptr_t( lhs && *lhs ? (*lhs)->m_next : 0 )
          , int( lhs && *lhs ? (*lhs)->m_task_type : 0 )
          , int( lhs && *lhs ? (*lhs)->m_ref_count : 0 )
          , uintptr_t(rhs)
          , uintptr_t( rhs ? rhs->m_next : 0 )
          , int( rhs ? rhs->m_task_type : 0 )
          , int( rhs ? rhs->m_ref_count : 0 )
          );
    fflush( stdout );
  }
#endif

      if ( *lhs ) decrement( *lhs );
      if ( rhs ) { Kokkos::atomic_increment( &(rhs->m_ref_count) ); }

      // Force write of *lhs

      *static_cast< task_lhs * volatile * >(lhs) = rhs ;

      Kokkos::memory_fence();
    }

  KOKKOS_FUNCTION
  size_t allocate_block_size( size_t n ); ///< Actual block size allocated

  KOKKOS_FUNCTION
  void * allocate( size_t n ); ///< Allocate from the memory pool

  KOKKOS_FUNCTION
  void deallocate( void * p , size_t n ); ///< Deallocate to the memory pool
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
class TaskBase< void , void , void > {
public:
  enum : int16_t   { TaskTeam = 0 , TaskSingle = 1 , Aggregate = 2 };
  enum : uintptr_t { LockTag = ~uintptr_t(0) , EndTag = ~uintptr_t(1) };
};

/** \brief  Base class for task management, access, and execution.
 *
 *  Inheritance structure to allow static_cast from the task root type
 *  and a task's FunctorType.
 *
 *    // Enable a Future to access result data
 *    TaskBase< Space , ResultType , void >
 *      : TaskBase< void , void , void >
 *      { ... };
 *
 *    // Enable a functor to access the base class
 *    TaskBase< Space , ResultType , FunctorType >
 *      : TaskBase< Space , ResultType , void >
 *      , FunctorType
 *      { ... };
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
 *      m_wait == LockTag: cannot add dependence
 *      m_next == LockTag: not a member of a wait queue
 *
 */
template< typename ExecSpace >
class TaskBase< ExecSpace , void , void >
{
public:

  enum : int16_t   { TaskTeam   = TaskBase<void,void,void>::TaskTeam
                   , TaskSingle = TaskBase<void,void,void>::TaskSingle
                   , Aggregate  = TaskBase<void,void,void>::Aggregate };

  enum : uintptr_t { LockTag = TaskBase<void,void,void>::LockTag
                   , EndTag  = TaskBase<void,void,void>::EndTag };

  using execution_space = ExecSpace ;
  using queue_type      = TaskQueue< execution_space > ;

  template< typename > friend class Kokkos::TaskScheduler ;

  typedef void (* function_type) ( TaskBase * , void * );

  // sizeof(TaskBase) == 48

  function_type  m_apply ;       ///< Apply function pointer
  queue_type   * m_queue ;       ///< Queue in which this task resides
  TaskBase     * m_wait ;        ///< Linked list of tasks waiting on this
  TaskBase     * m_next ;        ///< Waiting linked-list next
  int32_t        m_ref_count ;   ///< Reference count
  int32_t        m_alloc_size ;  ///< Allocation size
  int32_t        m_dep_count ;   ///< Aggregate's number of dependences
  int16_t        m_task_type ;   ///< Type of task
  int16_t        m_priority ;    ///< Priority of runnable task

  TaskBase() = delete ;
  TaskBase( TaskBase && ) = delete ;
  TaskBase( const TaskBase & ) = delete ;
  TaskBase & operator = ( TaskBase && ) = delete ;
  TaskBase & operator = ( const TaskBase & ) = delete ;

  KOKKOS_INLINE_FUNCTION ~TaskBase() = default ;

  // Constructor for a runnable task
  KOKKOS_INLINE_FUNCTION
  constexpr TaskBase( function_type arg_apply
                    , queue_type  * arg_queue
                    , TaskBase    * arg_dependence
                    , int           arg_ref_count
                    , int           arg_alloc_size
                    , int           arg_task_type
                    , int           arg_priority
                    ) noexcept
    : m_apply(      arg_apply )
    , m_queue(      arg_queue )
    , m_wait( 0 )
    , m_next(       arg_dependence )
    , m_ref_count(  arg_ref_count )
    , m_alloc_size( arg_alloc_size )
    , m_dep_count( 0 )
    , m_task_type(  arg_task_type )
    , m_priority(   arg_priority )
    {}

  // Constructor for an aggregate task
  KOKKOS_INLINE_FUNCTION
  constexpr TaskBase( queue_type  * arg_queue
                    , int           arg_ref_count
                    , int           arg_alloc_size
                    , int           arg_dep_count
                    ) noexcept
    : m_apply( 0 )
    , m_queue( arg_queue )
    , m_wait( 0 )
    , m_next( 0 )
    , m_ref_count(  arg_ref_count )
    , m_alloc_size( arg_alloc_size )
    , m_dep_count(  arg_dep_count )
    , m_task_type(  Aggregate )
    , m_priority( 0 )
    {}

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  TaskBase ** aggregate_dependences()
    { return reinterpret_cast<TaskBase**>( this + 1 ); }

  KOKKOS_INLINE_FUNCTION
  bool requested_respawn()
    {
      // This should only be called when a task has finished executing and is
      // in the transition to either the complete or executing-respawn state.
      TaskBase * const lock = reinterpret_cast< TaskBase * >( LockTag );
      return lock != m_next;
    }

  KOKKOS_INLINE_FUNCTION
  void add_dependence( TaskBase* dep )
    {
      // Precondition: lock == m_next

      TaskBase * const lock = (TaskBase *) LockTag ;

      // Assign dependence to m_next.  It will be processed in the subsequent
      // call to schedule.  Error if the dependence is reset.
      if ( lock != Kokkos::atomic_exchange( & m_next, dep ) ) {
        Kokkos::abort("TaskScheduler ERROR: resetting task dependence");
      }

      if ( 0 != dep ) {
        // The future may be destroyed upon returning from this call
        // so increment reference count to track this assignment.
        Kokkos::atomic_increment( &(dep->m_ref_count) );
      }
    }

  using get_return_type = void ;

  KOKKOS_INLINE_FUNCTION
  get_return_type get() const {}
};

template < typename ExecSpace , typename ResultType >
class TaskBase< ExecSpace , ResultType , void >
  : public TaskBase< ExecSpace , void , void >
{
private:

  using root_type     = TaskBase<ExecSpace,void,void> ;
  using function_type = typename root_type::function_type ;
  using queue_type    = typename root_type::queue_type ;

  static_assert( sizeof(root_type) == 48 , "" );

  TaskBase() = delete ;
  TaskBase( TaskBase && ) = delete ;
  TaskBase( const TaskBase & ) = delete ;
  TaskBase & operator = ( TaskBase && ) = delete ;
  TaskBase & operator = ( const TaskBase & ) = delete ;

public:

  ResultType   m_result ;

  KOKKOS_INLINE_FUNCTION ~TaskBase() = default ;

  // Constructor for runnable task
  KOKKOS_INLINE_FUNCTION
  constexpr TaskBase( function_type arg_apply
                    , queue_type  * arg_queue
                    , root_type   * arg_dependence
                    , int           arg_ref_count
                    , int           arg_alloc_size
                    , int           arg_task_type
                    , int           arg_priority
                    )
    : root_type( arg_apply 
               , arg_queue
               , arg_dependence
               , arg_ref_count
               , arg_alloc_size
               , arg_task_type
               , arg_priority
               )
    , m_result()
    {}

  using get_return_type = ResultType const & ;

  KOKKOS_INLINE_FUNCTION
  get_return_type get() const { return m_result ; }
};

template< typename ExecSpace , typename ResultType , typename FunctorType >
class TaskBase
  : public TaskBase< ExecSpace , ResultType , void >
  , public FunctorType
{
private:

  TaskBase() = delete ;
  TaskBase( TaskBase && ) = delete ;
  TaskBase( const TaskBase & ) = delete ;
  TaskBase & operator = ( TaskBase && ) = delete ;
  TaskBase & operator = ( const TaskBase & ) = delete ;

public:

  using root_type       = TaskBase< ExecSpace , void , void > ;
  using base_type       = TaskBase< ExecSpace , ResultType , void > ;
  using specialization  = TaskQueueSpecialization< ExecSpace > ;
  using function_type   = typename root_type::function_type ;
  using queue_type      = typename root_type::queue_type ;
  using member_type     = typename specialization::member_type ;
  using functor_type    = FunctorType ;
  using result_type     = ResultType ;

  template< typename Type >
  KOKKOS_INLINE_FUNCTION static
  void apply_functor
    ( Type * const task
    , typename std::enable_if
        < std::is_same< typename Type::result_type , void >::value
        , member_type * const
        >::type member
    )
    {
      using fType = typename Type::functor_type ;
      static_cast<fType*>(task)->operator()( *member );
    }

  template< typename Type >
  KOKKOS_INLINE_FUNCTION static
  void apply_functor
    ( Type * const task
    , typename std::enable_if
        < ! std::is_same< typename Type::result_type , void >::value
        , member_type * const
        >::type member
    )
    {
      using fType = typename Type::functor_type ;
      static_cast<fType*>(task)->operator()( *member , task->m_result );
    }

  KOKKOS_FUNCTION static
  void apply( root_type * root , void * exec )
    {
      TaskBase    * const task   = static_cast< TaskBase * >( root );
      member_type * const member = reinterpret_cast< member_type * >( exec );

      TaskBase::template apply_functor( task , member );

      // Task may be serial or team.
      // If team then must synchronize before querying if respawn was requested.
      // If team then only one thread calls destructor.

      member->team_barrier();

      if ( 0 == member->team_rank() && !(task->requested_respawn()) ) {
        // Did not respawn, destroy the functor to free memory.
        static_cast<functor_type*>(task)->~functor_type();
        // Cannot destroy and deallocate the task until its dependences
        // have been processed.
      }
    }

  // Constructor for runnable task
  KOKKOS_INLINE_FUNCTION
  constexpr TaskBase( function_type arg_apply
                    , queue_type  * arg_queue
                    , root_type   * arg_dependence
                    , int           arg_ref_count
                    , int           arg_alloc_size
                    , int           arg_task_type
                    , int           arg_priority
                    , FunctorType && arg_functor
                    )
    : base_type( arg_apply 
               , arg_queue
               , arg_dependence
               , arg_ref_count
               , arg_alloc_size
               , arg_task_type
               , arg_priority
               )
    , functor_type( arg_functor )
    {}

  KOKKOS_INLINE_FUNCTION
  ~TaskBase() {}
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_TASKQUEUE_HPP */
