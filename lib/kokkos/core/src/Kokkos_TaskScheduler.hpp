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

#ifndef KOKKOS_TASKSCHEDULER_HPP
#define KOKKOS_TASKSCHEDULER_HPP

//----------------------------------------------------------------------------

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_TASKDAG )

#include <Kokkos_Core_fwd.hpp>
//----------------------------------------------------------------------------

#include <Kokkos_MemoryPool.hpp>
#include <impl/Kokkos_Tags.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

// Forward declarations used in Impl::TaskQueue

template< typename Arg1 = void , typename Arg2 = void >
class Future ;

template< typename Space >
class TaskScheduler ;

template< typename Space >
void wait( TaskScheduler< Space > const & );

template< typename Space >
struct is_scheduler : public std::false_type {};

template< typename Space >
struct is_scheduler< TaskScheduler< Space > > : public std::true_type {};

} // namespace Kokkos

#include <impl/Kokkos_TaskQueue.hpp>

//----------------------------------------------------------------------------
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

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

/**
 *
 *  Future< space >  // value_type == void
 *  Future< value >  // space == Default
 *  Future< value , space >
 *
 */
template< typename Arg1 , typename Arg2 >
class Future {
private:

  template< typename > friend class TaskScheduler ;
  template< typename , typename > friend class Future ;
  template< typename , typename , typename > friend class Impl::TaskBase ;

  enum { Arg1_is_space  = Kokkos::is_space< Arg1 >::value };
  enum { Arg2_is_space  = Kokkos::is_space< Arg2 >::value };
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

  using task_base  = Impl::TaskBase< void , void , void > ;
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
  void clear()
    { if ( m_task ) queue_type::assign( & m_task , (task_base*)0 ); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  ~Future() { clear(); }

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
      clear();
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

      clear();
      m_task = rhs.m_task ;
      rhs.m_task = 0 ;
      return *this ;
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  int is_ready() const noexcept
    { return ( 0 == m_task ) || ( ((task_base*) task_base::LockTag) == m_task->m_wait ); }

  KOKKOS_INLINE_FUNCTION
  const typename Impl::TaskResult< ValueType >::reference_type
  get() const
    {
      if ( 0 == m_task ) {
        Kokkos::abort( "Kokkos:::Future::get ERROR: is_null()");
      }
      return Impl::TaskResult< ValueType >::get( m_task );
    }
};

// Is a Future with the given execution space
template< typename , typename ExecSpace = void >
struct is_future : public std::false_type {};

template< typename Arg1 , typename Arg2 , typename ExecSpace >
struct is_future< Future<Arg1,Arg2> , ExecSpace >
  : public std::integral_constant
      < bool ,
      ( std::is_same< ExecSpace , void >::value ||
        std::is_same< ExecSpace
                    , typename Future<Arg1,Arg2>::execution_space >::value )
      > {};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

enum class TaskPriority : int { High    = 0
                              , Regular = 1
                              , Low     = 2 };

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template< int TaskEnum , typename DepFutureType >
struct TaskPolicyData
{
  using execution_space = typename DepFutureType::execution_space ;
  using scheduler_type  = TaskScheduler< execution_space > ;

  enum : int { m_task_type = TaskEnum };

  scheduler_type const * m_scheduler ;
  DepFutureType  const   m_dependence ;
  int                    m_priority ;

  TaskPolicyData() = delete ;
  TaskPolicyData( TaskPolicyData && ) = default ;
  TaskPolicyData( TaskPolicyData const & ) = default ;
  TaskPolicyData & operator = ( TaskPolicyData && ) = default ;
  TaskPolicyData & operator = ( TaskPolicyData const & ) = default ;

  KOKKOS_INLINE_FUNCTION
  TaskPolicyData( DepFutureType        const & arg_future
                , Kokkos::TaskPriority const & arg_priority )
    : m_scheduler( 0 )
    , m_dependence( arg_future )
    , m_priority( static_cast<int>( arg_priority ) )
    {}

  KOKKOS_INLINE_FUNCTION
  TaskPolicyData( scheduler_type       const & arg_scheduler
                , Kokkos::TaskPriority const & arg_priority )
    : m_scheduler( & arg_scheduler )
    , m_dependence()
    , m_priority( static_cast<int>( arg_priority ) )
    {}

  KOKKOS_INLINE_FUNCTION
  TaskPolicyData( scheduler_type       const & arg_scheduler
                , DepFutureType        const & arg_future
                , Kokkos::TaskPriority const & arg_priority )
    : m_scheduler( & arg_scheduler )
    , m_dependence( arg_future )
    , m_priority( static_cast<int>( arg_priority ) )
    {}
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template< typename ExecSpace >
class TaskScheduler
{
private:

  using track_type = Kokkos::Impl::SharedAllocationTracker ;
  using queue_type = Kokkos::Impl::TaskQueue< ExecSpace > ;
  using task_base  = Impl::TaskBase< void , void , void > ;

  track_type   m_track ;
  queue_type * m_queue ;

  //----------------------------------------

public:

  using execution_space  = ExecSpace ;
  using memory_space     = typename queue_type::memory_space ;
  using memory_pool      = typename queue_type::memory_pool ;
  using member_type      =
    typename Kokkos::Impl::TaskQueueSpecialization< ExecSpace >::member_type ;

  KOKKOS_INLINE_FUNCTION
  TaskScheduler() : m_track(), m_queue(0) {}

  KOKKOS_INLINE_FUNCTION
  TaskScheduler( TaskScheduler && rhs ) = default ;

  KOKKOS_INLINE_FUNCTION
  TaskScheduler( TaskScheduler const & rhs ) = default ;

  KOKKOS_INLINE_FUNCTION
  TaskScheduler & operator = ( TaskScheduler && rhs ) = default ;

  KOKKOS_INLINE_FUNCTION
  TaskScheduler & operator = ( TaskScheduler const & rhs ) = default ;

  TaskScheduler( memory_pool const & arg_memory_pool )
    : m_track()
    , m_queue(0)
    {
      typedef Kokkos::Impl::SharedAllocationRecord
        < memory_space , typename queue_type::Destroy >
          record_type ;

      record_type * record =
        record_type::allocate( memory_space()
                             , "TaskQueue"
                             , sizeof(queue_type)
                             );

      m_queue = new( record->data() ) queue_type( arg_memory_pool );

      record->m_destroy.m_queue = m_queue ;

      m_track.assign_allocated_record_to_uninitialized( record );
    }

  TaskScheduler( memory_space const & arg_memory_space
               , size_t const mempool_capacity
               , unsigned const mempool_min_block_size  // = 1u << 6
               , unsigned const mempool_max_block_size  // = 1u << 10
               , unsigned const mempool_superblock_size // = 1u << 12
               )
    : TaskScheduler( memory_pool( arg_memory_space
                                , mempool_capacity
                                , mempool_min_block_size
                                , mempool_max_block_size
                                , mempool_superblock_size ) )
    {}

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  memory_pool * memory() const noexcept
    { return m_queue ? &( m_queue->m_memory ) : (memory_pool*) 0 ; }

  //----------------------------------------
  /**\brief  Allocation size for a spawned task */
  template< typename FunctorType >
  KOKKOS_FUNCTION
  size_t spawn_allocation_size() const
    { return m_queue->template spawn_allocation_size< FunctorType >(); }

  /**\brief  Allocation size for a when_all aggregate */
  KOKKOS_FUNCTION
  size_t when_all_allocation_size( int narg ) const
    { return m_queue->when_all_allocation_size( narg ); }

  //----------------------------------------

  template< int TaskEnum , typename DepFutureType , typename FunctorType >
  KOKKOS_FUNCTION static
  Kokkos::Future< typename FunctorType::value_type , execution_space >
  spawn( Impl::TaskPolicyData<TaskEnum,DepFutureType> const & arg_policy
       , typename task_base::function_type                    arg_function
       , FunctorType                                       && arg_functor
       )
    {
      using value_type  = typename FunctorType::value_type ;
      using future_type = Future< value_type , execution_space > ;
      using task_type   = Impl::TaskBase< execution_space
                                        , value_type
                                        , FunctorType > ;

      queue_type * const queue =
        arg_policy.m_scheduler ? arg_policy.m_scheduler->m_queue : (
        arg_policy.m_dependence.m_task
          ? static_cast<queue_type*>(arg_policy.m_dependence.m_task->m_queue)
          : (queue_type*) 0 );

      if ( 0 == queue ) {
        Kokkos::abort("Kokkos spawn requires scheduler or non-null Future");
      }

      if ( arg_policy.m_dependence.m_task != 0 &&
           arg_policy.m_dependence.m_task->m_queue != queue ) {
        Kokkos::abort("Kokkos spawn given incompatible scheduler and Future");
      }

      //----------------------------------------
      // Give single-thread back-ends an opportunity to clear
      // queue of ready tasks before allocating a new task

      queue->iff_single_thread_recursive_execute();

      //----------------------------------------

      future_type f ;

      // Allocate task from memory pool

      const size_t alloc_size =
        queue->template spawn_allocation_size< FunctorType >();

      f.m_task =
        reinterpret_cast< task_type * >(queue->allocate(alloc_size) );

      if ( f.m_task ) {

        // Placement new construction
        // Reference count starts at two:
        //   +1 for the matching decrement when task is complete
        //   +1 for the future
        new ( f.m_task ) task_type( std::move(arg_functor) );

        f.m_task->m_apply      = arg_function ;
        f.m_task->m_queue      = queue ;
        f.m_task->m_next       = arg_policy.m_dependence.m_task ;
        f.m_task->m_ref_count  = 2 ;
        f.m_task->m_alloc_size = alloc_size ;
        f.m_task->m_task_type  = arg_policy.m_task_type ;
        f.m_task->m_priority   = arg_policy.m_priority ;

        Kokkos::memory_fence();

        // The dependence (if any) is processed immediately
        // within the schedule function, as such the dependence's
        // reference count does not need to be incremented for
        // the assignment.

        queue->schedule_runnable( f.m_task );
        // This task may be updated or executed at any moment,
        // even during the call to 'schedule'.
      }

      return f ;
    }

  template< typename FunctorType , typename A1 , typename A2 >
  KOKKOS_FUNCTION static
  void
  respawn( FunctorType         * arg_self
         , Future<A1,A2> const & arg_dependence
         , TaskPriority  const & arg_priority
         )
    {
      // Precondition: task is in Executing state

      using value_type  = typename FunctorType::value_type ;
      using task_type   = Impl::TaskBase< execution_space
                                        , value_type
                                        , FunctorType > ;

      task_type * const task = static_cast< task_type * >( arg_self );

      task->m_priority = static_cast<int>(arg_priority);

      task->add_dependence( arg_dependence.m_task );

      // Postcondition: task is in Executing-Respawn state
    }

  template< typename FunctorType >
  KOKKOS_FUNCTION static
  void
  respawn( FunctorType         * arg_self
         , TaskScheduler const &
         , TaskPriority  const & arg_priority
         )
    {
      // Precondition: task is in Executing state

      using value_type  = typename FunctorType::value_type ;
      using task_type   = Impl::TaskBase< execution_space
                                        , value_type
                                        , FunctorType > ;

      task_type * const task = static_cast< task_type * >( arg_self );

      task->m_priority = static_cast<int>(arg_priority);

      task->add_dependence( (task_base*) 0 );

      // Postcondition: task is in Executing-Respawn state
    }

  //----------------------------------------
  /**\brief  Return a future that is complete
   *         when all input futures are complete.
   */
  template< typename A1 , typename A2 >
  KOKKOS_FUNCTION static
  Future< execution_space >
  when_all( Future< A1 , A2 > const arg[] , int narg )
    {
      using future_type = Future< execution_space > ;
      using task_base   = Kokkos::Impl::TaskBase< void , void , void > ;

      future_type f ;

      if ( narg ) {

        queue_type * queue = 0 ;

        for ( int i = 0 ; i < narg ; ++i ) {
          task_base * const t = arg[i].m_task ;
          if ( 0 != t ) {
            // Increment reference count to track subsequent assignment.
            Kokkos::atomic_increment( &(t->m_ref_count) );
            if ( queue == 0 ) {
              queue = static_cast< queue_type * >( t->m_queue );
            }
            else if ( queue != static_cast< queue_type * >( t->m_queue ) ) {
              Kokkos::abort("Kokkos when_all Futures must be in the same scheduler" );
            }
          }
        }

        if ( queue != 0 ) {

          size_t const alloc_size = queue->when_all_allocation_size( narg );

          f.m_task =
            reinterpret_cast< task_base * >( queue->allocate( alloc_size ) );

          if ( f.m_task ) {

            // Reference count starts at two:
            // +1 to match decrement when task completes
            // +1 for the future

            new( f.m_task ) task_base();

            f.m_task->m_queue      = queue ;
            f.m_task->m_ref_count  = 2 ;
            f.m_task->m_alloc_size = alloc_size ;
            f.m_task->m_dep_count  = narg ;
            f.m_task->m_task_type  = task_base::Aggregate ;

            // Assign dependences, reference counts were already incremented

            task_base * volatile * const dep =
              f.m_task->aggregate_dependences();

            for ( int i = 0 ; i < narg ; ++i ) { dep[i] = arg[i].m_task ; }

            Kokkos::memory_fence();

            queue->schedule_aggregate( f.m_task );
            // this when_all may be processed at any moment
          }
        }
      }

      return f ;
    }

  template < class F >
  KOKKOS_FUNCTION
  Future< execution_space >
  when_all( int narg , F const func )
    {
      using input_type  = decltype( func(0) );
      using future_type = Future< execution_space > ;
      using task_base   = Kokkos::Impl::TaskBase< void , void , void > ;

      static_assert( is_future< input_type >::value
                   , "Functor must return a Kokkos::Future" );

      future_type f ;

      if ( 0 == narg ) return f ;

      size_t const alloc_size = m_queue->when_all_allocation_size( narg );

      f.m_task =
        reinterpret_cast< task_base * >( m_queue->allocate( alloc_size ) );

      if ( f.m_task ) {

        // Reference count starts at two:
        // +1 to match decrement when task completes
        // +1 for the future

        new( f.m_task ) task_base();

        f.m_task->m_queue      = m_queue ;
        f.m_task->m_ref_count  = 2 ;
        f.m_task->m_alloc_size = alloc_size ;
        f.m_task->m_dep_count  = narg ;
        f.m_task->m_task_type  = task_base::Aggregate ;

        // Assign dependences, reference counts were already incremented

        task_base * volatile * const dep =
          f.m_task->aggregate_dependences();

        for ( int i = 0 ; i < narg ; ++i ) {
          const input_type arg_f = func(i);
          if ( 0 != arg_f.m_task ) {

            if ( m_queue != static_cast< queue_type * >( arg_f.m_task->m_queue ) ) {
              Kokkos::abort("Kokkos when_all Futures must be in the same scheduler" );
            }
            // Increment reference count to track subsequent assignment.
            Kokkos::atomic_increment( &(arg_f.m_task->m_ref_count) );
            dep[i] = arg_f.m_task ;
          }
        }

        Kokkos::memory_fence();

        m_queue->schedule_aggregate( f.m_task );
        // this when_all may be processed at any moment
      }
      return f ;
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  int allocation_capacity() const noexcept
    { return m_queue->m_memory.capacity(); }

  KOKKOS_INLINE_FUNCTION
  int allocated_task_count() const noexcept
    { return m_queue->m_count_alloc ; }

  KOKKOS_INLINE_FUNCTION
  int allocated_task_count_max() const noexcept
    { return m_queue->m_max_alloc ; }

  KOKKOS_INLINE_FUNCTION
  long allocated_task_count_accum() const noexcept
    { return m_queue->m_accum_alloc ; }

  //----------------------------------------

  template< typename S >
  friend
  void Kokkos::wait( Kokkos::TaskScheduler< S > const & );

};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//----------------------------------------------------------------------------
// Construct a TaskTeam execution policy

template< typename T >
Kokkos::Impl::TaskPolicyData
  < Kokkos::Impl::TaskBase<void,void,void>::TaskTeam
  , typename std::conditional< Kokkos::is_future< T >::value , T ,
    typename Kokkos::Future< typename T::execution_space > >::type
  >
KOKKOS_INLINE_FUNCTION
TaskTeam( T            const & arg
        , TaskPriority const & arg_priority = TaskPriority::Regular
        )
{
  static_assert( Kokkos::is_future<T>::value ||
                 Kokkos::is_scheduler<T>::value
               , "Kokkos TaskTeam argument must be Future or TaskScheduler" );

  return
    Kokkos::Impl::TaskPolicyData
      < Kokkos::Impl::TaskBase<void,void,void>::TaskTeam
      , typename std::conditional< Kokkos::is_future< T >::value , T ,
        typename Kokkos::Future< typename T::execution_space > >::type
      >( arg , arg_priority );
}

template< typename E , typename F >
Kokkos::Impl::
  TaskPolicyData< Kokkos::Impl::TaskBase<void,void,void>::TaskTeam , F >
KOKKOS_INLINE_FUNCTION
TaskTeam( TaskScheduler<E> const & arg_scheduler
        , F                const & arg_future
        , typename std::enable_if< Kokkos::is_future<F>::value ,
            TaskPriority >::type const & arg_priority = TaskPriority::Regular
        )
{
  return
    Kokkos::Impl::TaskPolicyData
      < Kokkos::Impl::TaskBase<void,void,void>::TaskTeam , F >
        ( arg_scheduler , arg_future , arg_priority );
}

// Construct a TaskSingle execution policy

template< typename T >
Kokkos::Impl::TaskPolicyData
  < Kokkos::Impl::TaskBase<void,void,void>::TaskSingle
  , typename std::conditional< Kokkos::is_future< T >::value , T ,
    typename Kokkos::Future< typename T::execution_space > >::type
  >
KOKKOS_INLINE_FUNCTION
TaskSingle( T            const & arg
          , TaskPriority const & arg_priority = TaskPriority::Regular
          )
{
  static_assert( Kokkos::is_future<T>::value ||
                 Kokkos::is_scheduler<T>::value
               , "Kokkos TaskSingle argument must be Future or TaskScheduler" );

  return
    Kokkos::Impl::TaskPolicyData
      < Kokkos::Impl::TaskBase<void,void,void>::TaskSingle
      , typename std::conditional< Kokkos::is_future< T >::value , T ,
        typename Kokkos::Future< typename T::execution_space > >::type
      >( arg , arg_priority );
}

template< typename E , typename F >
Kokkos::Impl::
  TaskPolicyData< Kokkos::Impl::TaskBase<void,void,void>::TaskSingle , F >
KOKKOS_INLINE_FUNCTION
TaskSingle( TaskScheduler<E> const & arg_scheduler
          , F                const & arg_future
          , typename std::enable_if< Kokkos::is_future<F>::value ,
              TaskPriority >::type const & arg_priority = TaskPriority::Regular
          )
{
  return
    Kokkos::Impl::TaskPolicyData
      < Kokkos::Impl::TaskBase<void,void,void>::TaskSingle , F >
        ( arg_scheduler , arg_future , arg_priority );
}

//----------------------------------------------------------------------------

/**\brief  A host control thread spawns a task with options
 *
 *  1) Team or Serial
 *  2) With scheduler or dependence
 *  3) High, Normal, or Low priority
 */
template< int TaskEnum
        , typename DepFutureType
        , typename FunctorType >
Future< typename FunctorType::value_type
      , typename DepFutureType::execution_space >
host_spawn( Impl::TaskPolicyData<TaskEnum,DepFutureType> const & arg_policy
          , FunctorType                                       && arg_functor
          )
{
  using exec_space = typename DepFutureType::execution_space ;
  using scheduler  = TaskScheduler< exec_space > ;

  typedef Impl::TaskBase< exec_space
                        , typename FunctorType::value_type
                        , FunctorType
                        > task_type ;

  static_assert( TaskEnum == task_type::TaskTeam ||
                 TaskEnum == task_type::TaskSingle
               , "Kokkos host_spawn requires TaskTeam or TaskSingle" );

  // May be spawning a Cuda task, must use the specialization
  // to query on-device function pointer.
  typename task_type::function_type const ptr =
    Kokkos::Impl::TaskQueueSpecialization< exec_space >::
      template get_function_pointer< task_type >();

  return scheduler::spawn( arg_policy , ptr , std::move(arg_functor) );
}

/**\brief  A task spawns a task with options
 *
 *  1) Team or Serial
 *  2) With scheduler or dependence
 *  3) High, Normal, or Low priority
 */
template< int TaskEnum
        , typename DepFutureType
        , typename FunctorType >
Future< typename FunctorType::value_type
      , typename DepFutureType::execution_space >
KOKKOS_INLINE_FUNCTION
task_spawn( Impl::TaskPolicyData<TaskEnum,DepFutureType> const & arg_policy
          , FunctorType                                       && arg_functor
          )
{
  using exec_space = typename DepFutureType::execution_space ;
  using scheduler  = TaskScheduler< exec_space > ;

  typedef Impl::TaskBase< exec_space
                        , typename FunctorType::value_type
                        , FunctorType
                        > task_type ;

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST ) && \
    defined( KOKKOS_ENABLE_CUDA )

  static_assert( ! std::is_same< Kokkos::Cuda , exec_space >::value
               , "Error calling Kokkos::task_spawn for Cuda space within Host code" );

#endif

  static_assert( TaskEnum == task_type::TaskTeam ||
                 TaskEnum == task_type::TaskSingle
               , "Kokkos host_spawn requires TaskTeam or TaskSingle" );

  typename task_type::function_type const ptr = task_type::apply ;

  return scheduler::spawn( arg_policy , ptr , std::move(arg_functor) );
}

/**\brief  A task respawns itself with options
 *
 *  1) With scheduler or dependence
 *  2) High, Normal, or Low priority
 */
template< typename FunctorType , typename T >
void
KOKKOS_INLINE_FUNCTION
respawn( FunctorType         * arg_self
       , T             const & arg
       , TaskPriority  const & arg_priority = TaskPriority::Regular
       )
{
  static_assert( Kokkos::is_future<T>::value ||
                 Kokkos::is_scheduler<T>::value
               , "Kokkos respawn argument must be Future or TaskScheduler" );

  TaskScheduler< typename T::execution_space >::
    respawn( arg_self , arg , arg_priority );
}

//----------------------------------------------------------------------------

template< typename A1 , typename A2 >
KOKKOS_INLINE_FUNCTION
Future< typename Future< A1 , A2 >::execution_space >
when_all( Future< A1 , A2 > const arg[]
        , int                     narg
        )
{
  return TaskScheduler< typename Future<A1,A2>::execution_space >::
    when_all( arg , narg );
}

//----------------------------------------------------------------------------
// Wait for all runnable tasks to complete

template< typename ExecSpace >
inline
void wait( TaskScheduler< ExecSpace > const & scheduler )
{ scheduler.m_queue->execute(); }

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_TASKSCHEDULER_HPP */

