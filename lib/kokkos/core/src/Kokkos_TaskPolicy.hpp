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

#ifndef KOKKOS_TASKPOLICY_HPP
#define KOKKOS_TASKPOLICY_HPP

//----------------------------------------------------------------------------

#include <Kokkos_Core_fwd.hpp>

// If compiling with CUDA then must be using CUDA 8 or better
// and use relocateable device code to enable the task policy.
// nvcc relocatable device code option: --relocatable-device-code=true

#if ( defined( KOKKOS_COMPILER_NVCC ) )
  #if ( 8000 <= CUDA_VERSION ) && \
      defined( KOKKOS_CUDA_USE_RELOCATABLE_DEVICE_CODE )

  #define KOKKOS_ENABLE_TASKPOLICY

  #endif
#else

#define KOKKOS_ENABLE_TASKPOLICY

#endif


#if defined( KOKKOS_ENABLE_TASKPOLICY )

//----------------------------------------------------------------------------

#include <Kokkos_MemoryPool.hpp>
#include <impl/Kokkos_Tags.hpp>
#include <impl/Kokkos_TaskQueue.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

enum TaskType { TaskTeam   = Impl::TaskBase<void,void,void>::TaskTeam
              , TaskSingle = Impl::TaskBase<void,void,void>::TaskSingle };

enum TaskPriority { TaskHighPriority    = 0
                  , TaskRegularPriority = 1
                  , TaskLowPriority     = 2 };

template< typename Space >
class TaskPolicy ;

template< typename Space >
void wait( TaskPolicy< Space > const & );

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/*\brief  Implementation data for task data management, access, and execution.
 *
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
template< typename Space , typename ResultType , typename FunctorType >
class TaskBase ;

template< typename Space >
class TaskExec ;

}} // namespace Kokkos::Impl

//----------------------------------------------------------------------------

namespace Kokkos {

/**
 *
 *  Future< space >  // value_type == void
 *  Future< value >  // space == Default
 *  Future< value , space >
 *
 */
template< typename Arg1 /* = void */ , typename Arg2 /* = void */ >
class Future {
private:

  template< typename > friend class TaskPolicy ;
  template< typename , typename > friend class Future ;
  template< typename , typename , typename > friend class Impl::TaskBase ;

  enum { Arg1_is_space  = Kokkos::Impl::is_space< Arg1 >::value };
  enum { Arg2_is_space  = Kokkos::Impl::is_space< Arg2 >::value };
  enum { Arg1_is_value  = ! Arg1_is_space &&
                          ! std::is_same< Arg1 , void >::value };
  enum { Arg2_is_value  = ! Arg2_is_space &&
                          ! std::is_same< Arg2 , void >::value };

  static_assert( ! ( Arg1_is_space && Arg2_is_space )
               , "Future cannot be given two spaces" );

  static_assert( ! ( Arg1_is_value && Arg2_is_value )
               , "Future cannot be given two value types" );

  using ValueType =
    typename std::conditional< Arg1_is_value , Arg1 ,
    typename std::conditional< Arg2_is_value , Arg2 , void
    >::type >::type ;

  using Space =
    typename std::conditional< Arg1_is_space , Arg1 ,
    typename std::conditional< Arg2_is_space , Arg2 , void
    >::type >::type ;

  using task_base  = Impl::TaskBase< Space , ValueType , void > ;
  using queue_type = Impl::TaskQueue< Space > ;

  task_base * m_task ;

  KOKKOS_INLINE_FUNCTION explicit
  Future( task_base * task ) : m_task(0)
    { if ( task ) queue_type::assign( & m_task , task ); }

  //----------------------------------------

public:

  using execution_space = typename Space::execution_space ;
  using value_type      = ValueType ;

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_task ; }

  KOKKOS_INLINE_FUNCTION
  int reference_count() const
    { return 0 != m_task ? m_task->reference_count() : 0 ; }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  ~Future() { if ( m_task ) queue_type::assign( & m_task , (task_base*)0 ); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  constexpr Future() noexcept : m_task(0) {}

  KOKKOS_INLINE_FUNCTION
  Future( Future && rhs )
    : m_task( rhs.m_task ) { rhs.m_task = 0 ; }

  KOKKOS_INLINE_FUNCTION
  Future( const Future & rhs )
    : m_task(0)
    { if ( rhs.m_task ) queue_type::assign( & m_task , rhs.m_task ); }

  KOKKOS_INLINE_FUNCTION
  Future & operator = ( Future && rhs )
    {
      if ( m_task ) queue_type::assign( & m_task , (task_base*)0 );
      m_task = rhs.m_task ;
      rhs.m_task = 0 ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION
  Future & operator = ( const Future & rhs )
    {
      if ( m_task || rhs.m_task ) queue_type::assign( & m_task , rhs.m_task );
      return *this ;
    }

  //----------------------------------------

  template< class A1 , class A2 >
  KOKKOS_INLINE_FUNCTION
  Future( Future<A1,A2> && rhs )
    : m_task( rhs.m_task )
    {
      static_assert
        ( std::is_same< Space , void >::value ||
          std::is_same< Space , typename Future<A1,A2>::Space >::value
        , "Assigned Futures must have the same space" );

      static_assert
        ( std::is_same< value_type , void >::value ||
          std::is_same< value_type , typename Future<A1,A2>::value_type >::value
        , "Assigned Futures must have the same value_type" );

      rhs.m_task = 0 ;
    }

  template< class A1 , class A2 >
  KOKKOS_INLINE_FUNCTION
  Future( const Future<A1,A2> & rhs )
    : m_task(0)
    {
      static_assert
        ( std::is_same< Space , void >::value ||
          std::is_same< Space , typename Future<A1,A2>::Space >::value
        , "Assigned Futures must have the same space" );

      static_assert
        ( std::is_same< value_type , void >::value ||
          std::is_same< value_type , typename Future<A1,A2>::value_type >::value
        , "Assigned Futures must have the same value_type" );

      if ( rhs.m_task ) queue_type::assign( & m_task , rhs.m_task );
    }

  template< class A1 , class A2 >
  KOKKOS_INLINE_FUNCTION
  Future & operator = ( const Future<A1,A2> & rhs )
    {
      static_assert
        ( std::is_same< Space , void >::value ||
          std::is_same< Space , typename Future<A1,A2>::Space >::value
        , "Assigned Futures must have the same space" );

      static_assert
        ( std::is_same< value_type , void >::value ||
          std::is_same< value_type , typename Future<A1,A2>::value_type >::value
        , "Assigned Futures must have the same value_type" );

      if ( m_task || rhs.m_task ) queue_type::assign( & m_task , rhs.m_task );
      return *this ;
    }

  template< class A1 , class A2 >
  KOKKOS_INLINE_FUNCTION
  Future & operator = ( Future<A1,A2> && rhs )
    {
      static_assert
        ( std::is_same< Space , void >::value ||
          std::is_same< Space , typename Future<A1,A2>::Space >::value
        , "Assigned Futures must have the same space" );

      static_assert
        ( std::is_same< value_type , void >::value ||
          std::is_same< value_type , typename Future<A1,A2>::value_type >::value
        , "Assigned Futures must have the same value_type" );

      if ( m_task ) queue_type::assign( & m_task , (task_base*) 0 );
      m_task = rhs.m_task ;
      rhs.m_task = 0 ;
      return *this ;
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  typename task_base::get_return_type
  get() const
    {
      if ( 0 == m_task ) {
        Kokkos::abort( "Kokkos:::Future::get ERROR: is_null()");
      }
      return m_task->get();
    }
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< typename ExecSpace >
class TaskPolicy
{
private:

  using track_type = Kokkos::Experimental::Impl::SharedAllocationTracker ;
  using queue_type = Kokkos::Impl::TaskQueue< ExecSpace > ;
  using task_base  = Impl::TaskBase< ExecSpace , void , void > ;

  track_type   m_track ;
  queue_type * m_queue ;

  //----------------------------------------
  // Process optional arguments to spawn and respawn functions

  KOKKOS_INLINE_FUNCTION static
  void assign( task_base * const ) {}

  // TaskTeam or TaskSingle
  template< typename ... Options >
  KOKKOS_INLINE_FUNCTION static
  void assign( task_base * const task
             , TaskType const & arg
             , Options const & ... opts )
    {
      task->m_task_type = arg ;
      assign( task , opts ... );
    }

  // TaskHighPriority or TaskRegularPriority or TaskLowPriority
  template< typename ... Options >
  KOKKOS_INLINE_FUNCTION static
  void assign( task_base * const task
             , TaskPriority const & arg
             , Options const & ... opts )
    {
      task->m_priority = arg ;
      assign( task , opts ... );
    }

  // Future for a dependence
  template< typename A1 , typename A2 , typename ... Options >
  KOKKOS_INLINE_FUNCTION static
  void assign( task_base * const task
             , Future< A1 , A2 > const & arg 
             , Options const & ... opts )
    {
      // Assign dependence to task->m_next
      // which will be processed within subsequent call to schedule.
      // Error if the dependence is reset.

      if ( 0 != Kokkos::atomic_exchange(& task->m_next, arg.m_task) ) {
        Kokkos::abort("TaskPolicy ERROR: resetting task dependence");
      }

      if ( 0 != arg.m_task ) {
        // The future may be destroyed upon returning from this call
        // so increment reference count to track this assignment.
        Kokkos::atomic_fetch_add( &(arg.m_task->m_ref_count) , 1 );
      }

      assign( task , opts ... );
    }

  //----------------------------------------

public:

  using execution_policy = TaskPolicy ;
  using execution_space  = ExecSpace ;
  using memory_space     = typename queue_type::memory_space ;
  using member_type      = Kokkos::Impl::TaskExec< ExecSpace > ;

  KOKKOS_INLINE_FUNCTION
  TaskPolicy() : m_track(), m_queue(0) {}

  KOKKOS_INLINE_FUNCTION
  TaskPolicy( TaskPolicy && rhs ) = default ;

  KOKKOS_INLINE_FUNCTION
  TaskPolicy( TaskPolicy const & rhs ) = default ;

  KOKKOS_INLINE_FUNCTION
  TaskPolicy & operator = ( TaskPolicy && rhs ) = default ;

  KOKKOS_INLINE_FUNCTION
  TaskPolicy & operator = ( TaskPolicy const & rhs ) = default ;

  TaskPolicy( memory_space const & arg_memory_space
            , unsigned const arg_memory_pool_capacity
            , unsigned const arg_memory_pool_log2_superblock = 12 )
    : m_track()
    , m_queue(0)
    {
      typedef Kokkos::Experimental::Impl::SharedAllocationRecord
        < memory_space , typename queue_type::Destroy >
          record_type ;

      record_type * record =
        record_type::allocate( arg_memory_space
                             , "TaskQueue"
                             , sizeof(queue_type)
                             );

      m_queue = new( record->data() )
        queue_type( arg_memory_space
                  , arg_memory_pool_capacity
                  , arg_memory_pool_log2_superblock );

      record->m_destroy.m_queue = m_queue ;

      m_track.assign_allocated_record_to_uninitialized( record );
    }

  //----------------------------------------
  /**\brief  Allocation size for a spawned task */
  template< typename FunctorType >
  KOKKOS_FUNCTION
  size_t spawn_allocation_size() const
    {
      using task_type  = Impl::TaskBase< execution_space
                                       , typename FunctorType::value_type
                                       , FunctorType > ;

      return m_queue->allocate_block_size( sizeof(task_type) );
    }

  /**\brief  Allocation size for a when_all aggregate */
  KOKKOS_FUNCTION
  size_t when_all_allocation_size( int narg ) const
    {
      using task_base  = Kokkos::Impl::TaskBase< ExecSpace , void , void > ;

      return m_queue->allocate_block_size( sizeof(task_base) + narg * sizeof(task_base*) );
    }

  //----------------------------------------

  /**\brief  A task spawns a task with options
   *
   *  1) High, Normal, or Low priority
   *  2) With or without dependence
   *  3) Team or Serial
   */
  template< typename FunctorType , typename ... Options >
  KOKKOS_FUNCTION
  Future< typename FunctorType::value_type , ExecSpace >
  task_spawn( FunctorType const & arg_functor 
            , Options const & ... arg_options
            ) const
    {
      using value_type  = typename FunctorType::value_type ;
      using future_type = Future< value_type , execution_space > ;
      using task_type   = Impl::TaskBase< execution_space
                                        , value_type
                                        , FunctorType > ;

      //----------------------------------------
      // Give single-thread back-ends an opportunity to clear
      // queue of ready tasks before allocating a new task

      m_queue->iff_single_thread_recursive_execute();

      //----------------------------------------

      future_type f ;

      // Allocate task from memory pool
      f.m_task =
        reinterpret_cast< task_type * >(m_queue->allocate(sizeof(task_type)));

      if ( f.m_task ) {

        // Placement new construction
        new ( f.m_task ) task_type( arg_functor );

        // Reference count starts at two
        // +1 for matching decrement when task is complete
        // +1 for future
        f.m_task->m_queue      = m_queue ;
        f.m_task->m_ref_count  = 2 ;
        f.m_task->m_alloc_size = sizeof(task_type);

        assign( f.m_task , arg_options... );

        // Spawning from within the execution space so the
        // apply function pointer is guaranteed to be valid
        f.m_task->m_apply = task_type::apply ;

        m_queue->schedule( f.m_task );
        // this task may be updated or executed at any moment
      }

      return f ;
    }

  /**\brief  The host process spawns a task with options
   *
   *  1) High, Normal, or Low priority
   *  2) With or without dependence
   *  3) Team or Serial
   */
  template< typename FunctorType , typename ... Options >
  inline
  Future< typename FunctorType::value_type , ExecSpace >
  host_spawn( FunctorType const & arg_functor 
            , Options const & ... arg_options
            ) const
    {
      using value_type  = typename FunctorType::value_type ;
      using future_type = Future< value_type , execution_space > ;
      using task_type   = Impl::TaskBase< execution_space
                                        , value_type
                                        , FunctorType > ;

      future_type f ;

      // Allocate task from memory pool
      f.m_task = 
        reinterpret_cast<task_type*>( m_queue->allocate(sizeof(task_type)) );

      if ( f.m_task ) {

        // Placement new construction
        new( f.m_task ) task_type( arg_functor );

        // Reference count starts at two:
        // +1 to match decrement when task completes
        // +1 for the future
        f.m_task->m_queue      = m_queue ;
        f.m_task->m_ref_count  = 2 ;
        f.m_task->m_alloc_size = sizeof(task_type);

        assign( f.m_task , arg_options... );

        // Potentially spawning outside execution space so the
        // apply function pointer must be obtained from execution space.
        // Required for Cuda execution space function pointer.
        queue_type::specialization::template
          proc_set_apply< FunctorType >( & f.m_task->m_apply );

        m_queue->schedule( f.m_task );
      }
      return f ;
    }

  /**\brief  Return a future that is complete
   *         when all input futures are complete.
   */
  template< typename A1 , typename A2 >
  KOKKOS_FUNCTION
  Future< ExecSpace >
  when_all( int narg , Future< A1 , A2 > const * const arg ) const
    {
      static_assert
        ( std::is_same< execution_space
                      , typename Future< A1 , A2 >::execution_space
                      >::value
        , "Future must have same execution space" );

      using future_type = Future< ExecSpace > ;
      using task_base   = Kokkos::Impl::TaskBase< ExecSpace , void , void > ;

      future_type f ;

      size_t const size  = sizeof(task_base) + narg * sizeof(task_base*);

      f.m_task =
        reinterpret_cast< task_base * >( m_queue->allocate( size ) );

      if ( f.m_task ) {

        new( f.m_task ) task_base();

        // Reference count starts at two:
        // +1 to match decrement when task completes
        // +1 for the future
        f.m_task->m_queue      = m_queue ;
        f.m_task->m_ref_count  = 2 ;
        f.m_task->m_alloc_size = size ;
        f.m_task->m_dep_count  = narg ;
        f.m_task->m_task_type  = task_base::Aggregate ;

        task_base ** const dep = f.m_task->aggregate_dependences();

        // Assign dependences to increment their reference count
        // The futures may be destroyed upon returning from this call
        // so increment reference count to track this assignment.

        for ( int i = 0 ; i < narg ; ++i ) {
          task_base * const t = dep[i] = arg[i].m_task ;
          if ( 0 != t ) {
            Kokkos::atomic_fetch_add( &(t->m_ref_count) , 1 );
          }
        }

        m_queue->schedule( f.m_task );
        // this when_all may be processed at any moment
      }

      return f ;
    }

  /**\brief  An executing task respawns itself with options
   *
   *  1) High, Normal, or Low priority
   *  2) With or without dependence
   */
  template< class FunctorType , typename ... Options >
  KOKKOS_FUNCTION
  void respawn( FunctorType * task_self
              , Options const & ... arg_options ) const
    {
      using value_type  = typename FunctorType::value_type ;
      using task_type   = Impl::TaskBase< execution_space
                                        , value_type
                                        , FunctorType > ;

      task_base * const zero = (task_base *) 0 ;
      task_base * const lock = (task_base *) task_base::LockTag ;
      task_type * const task = static_cast< task_type * >( task_self );

      // Precondition:
      //   task is in Executing state
      //   therefore  m_next == LockTag
      //
      // Change to m_next == 0 for no dependence

      if ( lock != Kokkos::atomic_exchange( & task->m_next, zero ) ) {
        Kokkos::abort("TaskPolicy::respawn ERROR: already respawned");
      }

      assign( task , arg_options... );

      // Postcondition:
      //   task is in Executing-Respawn state
      //   therefore  m_next == dependece or 0
    }

  //----------------------------------------

  template< typename S >
  friend
  void Kokkos::wait( Kokkos::TaskPolicy< S > const & );

  //----------------------------------------

  inline
  int allocation_capacity() const noexcept
    { return m_queue->m_memory.get_mem_size(); }

  KOKKOS_INLINE_FUNCTION
  int allocated_task_count() const noexcept
    { return m_queue->m_count_alloc ; }

  KOKKOS_INLINE_FUNCTION
  int allocated_task_count_max() const noexcept
    { return m_queue->m_max_alloc ; }

  KOKKOS_INLINE_FUNCTION
  long allocated_task_count_accum() const noexcept
    { return m_queue->m_accum_alloc ; }

};

template< typename ExecSpace >
inline
void wait( TaskPolicy< ExecSpace > const & policy )
{ policy.m_queue->execute(); }

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct FutureValueTypeIsVoidError {};

template < class ExecSpace , class ResultType , class FunctorType >
class TaskMember ;

} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

/**\brief  States of a task */
enum TaskState
  { TASK_STATE_NULL         = 0  ///<  Does not exist
  , TASK_STATE_CONSTRUCTING = 1  ///<  Is under construction
  , TASK_STATE_WAITING      = 2  ///<  Is waiting for execution
  , TASK_STATE_EXECUTING    = 4  ///<  Is executing
  , TASK_STATE_COMPLETE     = 8  ///<  Execution is complete
  };

/**\brief  Tag for Future<Latch,Space>
 */
struct Latch {};

/**
 *
 *  Future< space >  // value_type == void
 *  Future< value >  // space == Default
 *  Future< value , space >
 *
 */
template< class Arg1 = void , class Arg2 = void >
class Future {
private:

  template< class , class , class > friend class Impl::TaskMember ;
  template< class > friend class TaskPolicy ;
  template< class , class > friend class Future ;

  // Argument #2, if not void, must be the space.
  enum { Arg1_is_space  = Kokkos::Impl::is_execution_space< Arg1 >::value };
  enum { Arg2_is_space  = Kokkos::Impl::is_execution_space< Arg2 >::value };
  enum { Arg2_is_void   = std::is_same< Arg2 , void >::value };

  struct ErrorNoExecutionSpace {};

  enum { Opt1  =   Arg1_is_space && Arg2_is_void
       , Opt2  = ! Arg1_is_space && Arg2_is_void
       , Opt3  = ! Arg1_is_space && Arg2_is_space
       , OptOK = Kokkos::Impl::StaticAssert< Opt1 || Opt2 || Opt3 , ErrorNoExecutionSpace >::value
       };

  typedef typename
    Kokkos::Impl::if_c< Opt2 || Opt3 , Arg1 , void >::type
      ValueType ;

  typedef typename
    Kokkos::Impl::if_c< Opt1 , Arg1 , typename
    Kokkos::Impl::if_c< Opt2 , Kokkos::DefaultExecutionSpace , typename
    Kokkos::Impl::if_c< Opt3 , Arg2 , void
    >::type >::type >::type
      ExecutionSpace ;

  typedef Impl::TaskMember< ExecutionSpace , void , void >       TaskRoot ;
  typedef Impl::TaskMember< ExecutionSpace , ValueType , void >  TaskValue ;

  TaskRoot * m_task ;

  KOKKOS_INLINE_FUNCTION explicit
  Future( TaskRoot * task )
    : m_task(0)
    { TaskRoot::assign( & m_task , TaskRoot::template verify_type< ValueType >( task ) ); }

  //----------------------------------------

public:

  typedef ValueType       value_type;
  typedef ExecutionSpace  execution_space ;

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  TaskState get_task_state() const
    { return 0 != m_task ? m_task->get_state() : TASK_STATE_NULL ; }

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_task ; }

  KOKKOS_INLINE_FUNCTION
  int reference_count() const
    { return 0 != m_task ? m_task->reference_count() : 0 ; }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  ~Future() { TaskRoot::assign( & m_task , 0 ); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  Future() : m_task(0) {}

  KOKKOS_INLINE_FUNCTION
  Future( const Future & rhs )
    : m_task(0)
    { TaskRoot::assign( & m_task , rhs.m_task ); }

  KOKKOS_INLINE_FUNCTION
  Future & operator = ( const Future & rhs )
    { TaskRoot::assign( & m_task , rhs.m_task ); return *this ; }

  //----------------------------------------

  template< class A1 , class A2 >
  KOKKOS_INLINE_FUNCTION
  Future( const Future<A1,A2> & rhs )
    : m_task(0)
    { TaskRoot::assign( & m_task , TaskRoot::template verify_type< value_type >( rhs.m_task ) ); }

  template< class A1 , class A2 >
  KOKKOS_INLINE_FUNCTION
  Future & operator = ( const Future<A1,A2> & rhs )
    { TaskRoot::assign( & m_task , TaskRoot::template verify_type< value_type >( rhs.m_task ) ); return *this ; }

  //----------------------------------------

  typedef typename TaskValue::get_result_type get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const
    {
      if ( 0 == m_task ) {
        Kokkos::abort( "Kokkos::Experimental::Future::get ERROR: is_null()");
      }
      return static_cast<TaskValue*>( m_task )->get();  
    }

  //----------------------------------------
};

template< class Arg2 >
class Future< Latch , Arg2 > {
private:

  template< class , class , class > friend class Impl::TaskMember ;
  template< class > friend class TaskPolicy ;
  template< class , class > friend class Future ;

  // Argument #2, if not void, must be the space.
  enum { Arg2_is_space  = Kokkos::Impl::is_execution_space< Arg2 >::value };
  enum { Arg2_is_void   = std::is_same< Arg2 , void >::value };

  static_assert( Arg2_is_space || Arg2_is_void 
               , "Future template argument #2 must be a space" );

  typedef typename
    std::conditional< Arg2_is_space , Arg2 , Kokkos::DefaultExecutionSpace >
     ::type ExecutionSpace ;

  typedef Impl::TaskMember< ExecutionSpace , void , void >  TaskRoot ;

  TaskRoot * m_task ;

  KOKKOS_INLINE_FUNCTION explicit
  Future( TaskRoot * task )
    : m_task(0)
    { TaskRoot::assign( & m_task , task ); }

  //----------------------------------------

public:

  typedef void            value_type;
  typedef ExecutionSpace  execution_space ;

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  void add( const int k ) const
    { if ( 0 != m_task ) m_task->latch_add(k); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  TaskState get_task_state() const
    { return 0 != m_task ? m_task->get_state() : TASK_STATE_NULL ; }

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_task ; }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  ~Future() { TaskRoot::assign( & m_task , 0 ); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  Future() : m_task(0) {}

  KOKKOS_INLINE_FUNCTION
  Future( const Future & rhs )
    : m_task(0)
    { TaskRoot::assign( & m_task , rhs.m_task ); }

  KOKKOS_INLINE_FUNCTION
  Future & operator = ( const Future & rhs )
    { TaskRoot::assign( & m_task , rhs.m_task ); return *this ; }

  //----------------------------------------

  typedef void get_result_type ;

  KOKKOS_INLINE_FUNCTION
  void get() const {}

  //----------------------------------------

};

namespace Impl {

template< class T >
struct is_future : public std::false_type {};

template< class Arg0 , class Arg1 >
struct is_future< Kokkos::Experimental::Future<Arg0,Arg1> >
  : public std::true_type {};

} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

/** \brief  If the argument is an execution space then a serial task in that space */
template< class Arg0 = Kokkos::DefaultExecutionSpace >
class TaskPolicy {
public:

  typedef typename Arg0::execution_space  execution_space ;

  //----------------------------------------

  TaskPolicy
    ( const unsigned arg_task_max_count
    , const unsigned arg_task_max_size
    , const unsigned arg_task_default_dependence_capacity = 4
    , const unsigned arg_task_team_size = 0 /* choose default */
    );

  TaskPolicy() = default ;
  TaskPolicy( TaskPolicy && rhs ) = default ;
  TaskPolicy( const TaskPolicy & rhs ) = default ;
  TaskPolicy & operator = ( TaskPolicy && rhs ) = default ;
  TaskPolicy & operator = ( const TaskPolicy & rhs ) = default ;

  //----------------------------------------
  /** \brief  Create a serial task with storage for dependences.
   *
   *  Postcondition: Task is in the 'constructing' state.
   */
  template< class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  create( const FunctorType & functor
        , const unsigned      dependence_capacity /* = default */ );

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< typename FunctorType::value_type , execution_space >
  create_team( const FunctorType & functor
             , const unsigned dependence_capacity /* = default */ );

  /** \brief  Set dependence that 'after' cannot start execution
   *          until 'before' has completed.
   *
   *  Precondition: The 'after' task must be in then 'Constructing' state.
   */
  template< class TA , class TB >
  void add_dependence( const Future<TA,execution_space> & after
                     , const Future<TB,execution_space> & before ) const ;

  /** \brief  Spawn a task in the 'Constructing' state
   *
   *  Precondition:  Task is in the 'constructing' state.
   *  Postcondition: Task is waiting, executing, or complete.
   */
  template< class T >
  const Future<T,execution_space> &
  spawn( const Future<T,execution_space> & ) const ;

  //----------------------------------------
  /** \brief  Query dependence of an executing task */

  template< class FunctorType >
  Future< execution_space >
  get_dependence( FunctorType * , const int ) const ;

  //----------------------------------------
  /** \brief  Clear current dependences of an executing task
   *          in preparation for setting new dependences and
   *          respawning.
   *
   * Precondition: The functor must be a task in the executing state.
   */
  template< class FunctorType >
  void clear_dependence( FunctorType * ) const ;

  /** \brief  Set dependence that 'after' cannot resume execution
   *          until 'before' has completed.
   *
   *  The 'after' functor must be in the executing state
   */
  template< class FunctorType , class TB >
  void add_dependence( FunctorType * after
                     , const Future<TB,execution_space> & before ) const ;

  /** \brief  Respawn (reschedule) an executing task to be called again
   *          after all dependences have completed.
   */
  template< class FunctorType >
  void respawn( FunctorType * ) const ;
};

//----------------------------------------------------------------------------
/** \brief  Create and spawn a single-thread task */
template< class ExecSpace , class FunctorType >
inline
Future< typename FunctorType::value_type , ExecSpace >
spawn( TaskPolicy<ExecSpace> & policy , const FunctorType & functor )
{ return policy.spawn( policy.create( functor ) ); }

/** \brief  Create and spawn a single-thread task with dependences */
template< class ExecSpace , class FunctorType , class Arg0 , class Arg1 >
inline
Future< typename FunctorType::value_type , ExecSpace >
spawn( TaskPolicy<ExecSpace>   & policy
     , const FunctorType       & functor
     , const Future<Arg0,Arg1> & before_0
     , const Future<Arg0,Arg1> & before_1 )
{
  Future< typename FunctorType::value_type , ExecSpace > f ;
  f = policy.create( functor , 2 );
  policy.add_dependence( f , before_0 );
  policy.add_dependence( f , before_1 );
  policy.spawn( f );
  return f ;
}

//----------------------------------------------------------------------------
/** \brief  Create and spawn a parallel_for task */
template< class ExecSpace , class ParallelPolicyType , class FunctorType >
inline
Future< typename FunctorType::value_type , ExecSpace >
spawn_foreach( TaskPolicy<ExecSpace>     & task_policy
             , const ParallelPolicyType  & parallel_policy
             , const FunctorType         & functor )
{ return task_policy.spawn( task_policy.create_foreach( parallel_policy , functor ) ); }

/** \brief  Create and spawn a parallel_reduce task */
template< class ExecSpace , class ParallelPolicyType , class FunctorType >
inline
Future< typename FunctorType::value_type , ExecSpace >
spawn_reduce( TaskPolicy<ExecSpace>     & task_policy
            , const ParallelPolicyType  & parallel_policy
            , const FunctorType         & functor )
{ return task_policy.spawn( task_policy.create_reduce( parallel_policy , functor ) ); }

//----------------------------------------------------------------------------
/** \brief  Respawn a task functor with dependences */
template< class ExecSpace , class FunctorType , class Arg0 , class Arg1 >
inline
void respawn( TaskPolicy<ExecSpace>   & policy
            , FunctorType *             functor
            , const Future<Arg0,Arg1> & before_0
            , const Future<Arg0,Arg1> & before_1
            )
{
  policy.clear_dependence( functor );
  policy.add_dependence( functor , before_0 );
  policy.add_dependence( functor , before_1 );
  policy.respawn( functor );
}

//----------------------------------------------------------------------------

template< class ExecSpace >
void wait( TaskPolicy< ExecSpace > & );

} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKPOLICY ) */
#endif /* #ifndef KOKKOS_TASKPOLICY_HPP */

