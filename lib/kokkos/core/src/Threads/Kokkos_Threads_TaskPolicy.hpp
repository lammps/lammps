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

#ifndef KOKKOS_THREADS_TASKPOLICY_HPP
#define KOKKOS_THREADS_TASKPOLICY_HPP


#include <Kokkos_Threads.hpp>
#include <Kokkos_TaskPolicy.hpp>

#if defined( KOKKOS_HAVE_PTHREAD )

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct ThreadsTaskPolicyQueue ;

/** \brief  Base class for all Kokkos::Threads tasks */
template<>
class TaskMember< Kokkos::Threads , void , void > {
public:

  template < class > friend class Kokkos::Experimental::TaskPolicy ;
  friend struct ThreadsTaskPolicyQueue ;

  typedef TaskMember * (* function_verify_type) ( TaskMember * );
  typedef void         (* function_single_type) ( TaskMember * );
  typedef void         (* function_team_type)   ( TaskMember * , Kokkos::Impl::ThreadsExecTeamMember & );

private:


  ThreadsTaskPolicyQueue * m_policy ;
  TaskMember * volatile  * m_queue ;
  function_verify_type     m_verify ;
  function_team_type       m_team ;         ///< Apply function
  function_single_type     m_serial ;       ///< Apply function
  TaskMember **            m_dep ;          ///< Dependences
  TaskMember *             m_wait ;         ///< Head of linked list of tasks waiting on this task
  TaskMember *             m_next ;         ///< Member of linked list of tasks
  int                      m_dep_capacity ; ///< Capacity of dependences
  int                      m_dep_size ;     ///< Actual count of dependences
  int                      m_size_alloc ;
  int                      m_shmem_size ;
  int                      m_ref_count ;    ///< Reference count
  int                      m_state ;        ///< State of the task


  TaskMember( TaskMember && ) = delete ;
  TaskMember( const TaskMember & ) = delete ;
  TaskMember & operator = ( TaskMember && ) = delete ;
  TaskMember & operator = ( const TaskMember & ) = delete ;

protected:

  TaskMember()
    : m_policy(0)
    , m_verify(0)
    , m_team(0)
    , m_serial(0)
    , m_dep(0)
    , m_wait(0)
    , m_next(0)
    , m_dep_capacity(0)
    , m_dep_size(0)
    , m_size_alloc(0)
    , m_shmem_size(0)
    , m_ref_count(0)
    , m_state( TASK_STATE_CONSTRUCTING )
    {}

public:

  ~TaskMember();

  KOKKOS_INLINE_FUNCTION
  int reference_count() const
    { return *((volatile int *) & m_ref_count ); }

  template< typename ResultType >
  KOKKOS_FUNCTION static
  TaskMember * verify_type( TaskMember * t )
    {
      enum { check_type = ! std::is_same< ResultType , void >::value };

      if ( check_type && t != 0 ) {

        // Verify that t->m_verify is this function
        const function_verify_type self = & TaskMember::template verify_type< ResultType > ;

        if ( t->m_verify != self ) {
          t = 0 ;
          Kokkos::abort("TaskPolicy< Threads > verify_result_type" );
        }
      }
      return t ;
    }

  //----------------------------------------
  /*  Inheritence Requirements on task types:
   *
   *    class TaskMember< Threads , DerivedType::value_type , FunctorType >
   *      : public TaskMember< Threads , DerivedType::value_type , void >
   *      , public Functor
   *      { ... };
   *
   *  If value_type != void
   *    class TaskMember< Threads , value_type , void >
   *      : public TaskMember< Threads , void , void >
   *
   */
  //----------------------------------------

  template< class DerivedTaskType , class Tag >
  KOKKOS_FUNCTION static
  void apply_single(
    typename std::enable_if
      <( std::is_same<Tag,void>::value &&
         std::is_same< typename DerivedTaskType::result_type , void >::value
       ), TaskMember * >::type t )
    {
      {
        typedef typename DerivedTaskType::functor_type  functor_type ;

        functor_type * const f = 
          static_cast< functor_type * >( static_cast< DerivedTaskType * >(t) );

        f->apply();

        if ( t->m_state == int(Kokkos::Experimental::TASK_STATE_EXECUTING) ) {
          f->~functor_type();
        }
      }
    }

  template< class DerivedTaskType , class Tag >
  KOKKOS_FUNCTION static
  void apply_single(
    typename std::enable_if
      <( std::is_same< Tag , void >::value &&
         ! std::is_same< typename DerivedTaskType::result_type , void >::value
       ), TaskMember * >::type t )
    {
      {
        typedef typename DerivedTaskType::functor_type  functor_type ;

        DerivedTaskType * const self = static_cast< DerivedTaskType * >(t);
        functor_type    * const f    = static_cast< functor_type * >( self );

        f->apply( self->m_result );

        if ( t->m_state == int(Kokkos::Experimental::TASK_STATE_EXECUTING) ) {
          f->~functor_type();
        }
      }
    }

  //----------------------------------------

  template< class DerivedTaskType , class Tag >
  KOKKOS_FUNCTION static
  void apply_team(
    typename std::enable_if
      <( std::is_same<Tag,void>::value &&
         std::is_same<typename DerivedTaskType::result_type,void>::value
       ), TaskMember * >::type t
    , Kokkos::Impl::ThreadsExecTeamMember & member
    )
    {
      typedef typename DerivedTaskType::functor_type  functor_type ;

      functor_type * const f =
        static_cast< functor_type * >( static_cast< DerivedTaskType * >(t) );
    
      f->apply( member );

      // Synchronize for possible functor destruction and
      // completion of team task.
      if ( member.team_fan_in() ) {
        if ( t->m_state == int(Kokkos::Experimental::TASK_STATE_EXECUTING) ) {
          f->~functor_type();
        }
      }
    }

  template< class DerivedTaskType , class Tag >
  KOKKOS_FUNCTION static
  void apply_team(
    typename std::enable_if
      <( std::is_same<Tag,void>::value &&
         ! std::is_same<typename DerivedTaskType::result_type,void>::value
       ), TaskMember * >::type t
    , Kokkos::Impl::ThreadsExecTeamMember & member
    )
    {
      typedef typename DerivedTaskType::functor_type  functor_type ;

      DerivedTaskType * const self = static_cast< DerivedTaskType * >(t);
      functor_type    * const f    = static_cast< functor_type * >( self );
    
      f->apply( member , self->m_result );

      // Synchronize for possible functor destruction and
      // completion of team task.
      if ( member.team_fan_in() ) {
        if ( t->m_state == int(Kokkos::Experimental::TASK_STATE_EXECUTING) ) {
          f->~functor_type();
        }
      }
    }

  //----------------------------------------

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  static
  void assign( TaskMember ** const lhs , TaskMember * const rhs );
#else
  KOKKOS_INLINE_FUNCTION static
  void assign( TaskMember ** const lhs , TaskMember * const rhs ) {}
#endif

  TaskMember * get_dependence( int i ) const ;

  KOKKOS_INLINE_FUNCTION
  int get_dependence() const { return m_dep_size ; }

  void clear_dependence();

  void latch_add( const int k );

  //----------------------------------------

  typedef FutureValueTypeIsVoidError get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const { return get_result_type() ; }

  inline static
  void construct_result( TaskMember * const ) {}

  KOKKOS_INLINE_FUNCTION
  Kokkos::Experimental::TaskState get_state() const { return Kokkos::Experimental::TaskState( m_state ); }

};

/** \brief  A Future< Kokkos::Threads , ResultType > will cast
 *          from  TaskMember< Kokkos::Threads , void , void >
 *          to    TaskMember< Kokkos::Threads , ResultType , void >
 *          to query the result.
 */
template< class ResultType >
class TaskMember< Kokkos::Threads , ResultType , void >
  : public TaskMember< Kokkos::Threads , void , void >
{
public:

  typedef ResultType result_type ;

  result_type  m_result ;

  typedef const result_type & get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const { return m_result ; }

  inline static
  void construct_result( TaskMember * const ptr )
    { 
      new((void*)(& ptr->m_result)) result_type();
    }

  inline
  TaskMember() : TaskMember< Kokkos::Threads , void , void >(), m_result() {}

  TaskMember( TaskMember && ) = delete ;
  TaskMember( const TaskMember & ) = delete ;
  TaskMember & operator = ( TaskMember && ) = delete ;
  TaskMember & operator = ( const TaskMember & ) = delete ;
};

/** \brief  Callback functions will cast
 *          from  TaskMember< Kokkos::Threads , void , void >
 *          to    TaskMember< Kokkos::Threads , ResultType , FunctorType >
 *          to execute work functions.
 */
template< class ResultType , class FunctorType >
class TaskMember< Kokkos::Threads , ResultType , FunctorType >
  : public TaskMember< Kokkos::Threads , ResultType , void >
  , public FunctorType
{
public:
  typedef ResultType   result_type ;
  typedef FunctorType  functor_type ;

  inline
  TaskMember( const functor_type & arg_functor )
    : TaskMember< Kokkos::Threads , ResultType , void >()
    , functor_type( arg_functor )
    {}

  inline static
  void copy_construct( TaskMember * const ptr
                     , const functor_type & arg_functor )
    {
      typedef TaskMember< Kokkos::Threads , ResultType , void > base_type ;

      new((void*)static_cast<FunctorType*>(ptr)) functor_type( arg_functor );

      base_type::construct_result( static_cast<base_type*>( ptr ) );
    }

  TaskMember() = delete ;
  TaskMember( TaskMember && ) = delete ;
  TaskMember( const TaskMember & ) = delete ;
  TaskMember & operator = ( TaskMember && ) = delete ;
  TaskMember & operator = ( const TaskMember & ) = delete ;
};

//----------------------------------------------------------------------------

struct ThreadsTaskPolicyQueue {

  enum { NPRIORITY = 3 };

  typedef Kokkos::Experimental::MemoryPool< Kokkos::Threads >
    memory_space ;

  typedef Kokkos::Experimental::Impl::TaskMember< Kokkos::Threads, void, void >
    task_root_type ;

  memory_space     m_space ;
  task_root_type * m_team[ NPRIORITY ];
  task_root_type * m_serial[ NPRIORITY ];
  int              m_team_size ;    ///< Fixed size of a task-team
  int              m_default_dependence_capacity ;
  int     volatile m_count_ready ;  ///< Ready plus executing tasks
  int     volatile m_count_alloc ;  ///< Total allocated tasks

  // Execute tasks until all non-waiting tasks are complete.
  static void driver( Kokkos::Impl::ThreadsExec & exec
                    , const void * arg );

  task_root_type * allocate_task
   ( const unsigned arg_sizeof_task
   , const unsigned arg_dep_capacity
   , const unsigned arg_team_shmem
   );

  void deallocate_task( void * , unsigned );
  void schedule_task( task_root_type * const
                    , const bool initial_spawn = true );
  void reschedule_task( task_root_type * const );
  void add_dependence( task_root_type * const after
                     , task_root_type * const before );

  // When a task finishes executing update its dependences
  // and either deallocate the task if complete
  // or reschedule the task if respawned.
  void complete_executed_task( task_root_type * );

  // Pop a task from a ready queue
  static task_root_type *
    pop_ready_task( task_root_type * volatile * const queue );

  ThreadsTaskPolicyQueue() = delete ;
  ThreadsTaskPolicyQueue( ThreadsTaskPolicyQueue && ) = delete ;
  ThreadsTaskPolicyQueue( const ThreadsTaskPolicyQueue & ) = delete ;
  ThreadsTaskPolicyQueue & operator = ( ThreadsTaskPolicyQueue && ) = delete ;
  ThreadsTaskPolicyQueue & operator = ( const ThreadsTaskPolicyQueue & ) = delete ;

  ~ThreadsTaskPolicyQueue();

  ThreadsTaskPolicyQueue
    ( const unsigned arg_task_max_count
    , const unsigned arg_task_max_size
    , const unsigned arg_task_default_dependence_capacity
    , const unsigned arg_task_team_size
    );

  // Callback to destroy the shared memory tracked queue.
  struct Destroy {
    ThreadsTaskPolicyQueue * m_policy ;
    void destroy_shared_allocation();
  };
};

} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

void wait( TaskPolicy< Kokkos::Threads > & );

template<>
class TaskPolicy< Kokkos::Threads >
{
public:

  typedef Kokkos::Threads                      execution_space ;
  typedef TaskPolicy                           execution_policy ;
  typedef Kokkos::Impl::ThreadsExecTeamMember  member_type ;

private:

  typedef Impl::TaskMember< Kokkos::Threads , void , void >  task_root_type ;
  typedef Kokkos::Experimental::MemoryPool< Kokkos::Threads > memory_space ;

  typedef Kokkos::Experimental::Impl::SharedAllocationTracker track_type ;

  track_type                      m_track ;
  Impl::ThreadsTaskPolicyQueue  * m_policy ;

  template< class FunctorType >
  static inline
  const task_root_type * get_task_root( const FunctorType * f )
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType > task_type ;
      return static_cast< const task_root_type * >( static_cast< const task_type * >(f) );
    }

  template< class FunctorType >
  static inline
  task_root_type * get_task_root( FunctorType * f )
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType > task_type ;
      return static_cast< task_root_type * >( static_cast< task_type * >(f) );
    }

  /** \brief  Allocate and construct a task.
   *
   *  Allocate space for DerivedTaskType followed by TaskMember*[ dependence_capacity ]
   */
  template< class DerivedTaskType , class Tag >
  task_root_type *
  create( const typename DerivedTaskType::functor_type &  arg_functor
        , const task_root_type::function_single_type      arg_apply_single
        , const task_root_type::function_team_type        arg_apply_team
        , const unsigned                                  arg_team_shmem
        , const unsigned                                  arg_dependence_capacity
        )
    {
      task_root_type * const t =
        m_policy->allocate_task( sizeof(DerivedTaskType)
                               , arg_dependence_capacity
                               , arg_team_shmem
                               );
      if ( t != 0 ) {

        DerivedTaskType * const task = static_cast<DerivedTaskType*>(t);

        DerivedTaskType::copy_construct( task , arg_functor );

        task->task_root_type::m_verify  = & task_root_type::template verify_type< typename DerivedTaskType::value_type > ;
        task->task_root_type::m_team    = arg_apply_team ;
        task->task_root_type::m_serial  = arg_apply_single ;

        // Do not proceed until initialization is written to memory
        Kokkos::memory_fence();
      }
      return t ;
    }

public:

  // Valid team sizes are 1,
  // Threads::pool_size(1) == threads per numa, or
  // Threads::pool_size(2) == threads per core

  TaskPolicy
    ( const unsigned arg_task_max_count
    , const unsigned arg_task_max_size
    , const unsigned arg_task_default_dependence_capacity = 4
    , const unsigned arg_task_team_size = 0 /* choose default */
    );

  KOKKOS_FUNCTION TaskPolicy() = default ;
  KOKKOS_FUNCTION TaskPolicy( TaskPolicy && rhs ) = default ;
  KOKKOS_FUNCTION TaskPolicy( const TaskPolicy & rhs ) = default ;
  KOKKOS_FUNCTION TaskPolicy & operator = ( TaskPolicy && rhs ) = default ;
  KOKKOS_FUNCTION TaskPolicy & operator = ( const TaskPolicy & rhs ) = default ;

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  int allocated_task_count() const { return m_policy->m_count_alloc ; }

  //----------------------------------------
  // Create serial-thread task

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< typename FunctorType::value_type , execution_space >
  task_create( const FunctorType & functor
             , const unsigned dependence_capacity = ~0u )
    {
      typedef typename FunctorType::value_type  value_type ;
      typedef Impl::TaskMember< execution_space , value_type , FunctorType >  task_type ;

      return Future< value_type , execution_space >(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        TaskPolicy::create< task_type , void >
          ( functor
          , & task_root_type::template apply_single< task_type , void >
          , task_root_type::function_team_type(0)
          , 0
          , dependence_capacity
          )
#endif
        );
    }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< typename FunctorType::value_type , execution_space >
  proc_create( const FunctorType & functor
             , const unsigned dependence_capacity = ~0u )
    { return task_create( functor , dependence_capacity ); }

  // Create thread-team task

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< typename FunctorType::value_type , execution_space >
  task_create_team( const FunctorType & functor
                  , const unsigned dependence_capacity = ~0u )
    {
      typedef typename FunctorType::value_type  value_type ;
      typedef Impl::TaskMember< execution_space , value_type , FunctorType >  task_type ;

      return Future< value_type , execution_space >(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        TaskPolicy::create< task_type , void >
          ( functor
          , task_root_type::function_single_type(0)
          , & task_root_type::template apply_team< task_type , void >
          , Kokkos::Impl::FunctorTeamShmemSize< FunctorType >::
              value( functor , m_policy->m_team_size )
          , dependence_capacity
          )
#endif
        );
    }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< typename FunctorType::value_type , execution_space >
  proc_create_team( const FunctorType & functor
                  , const unsigned dependence_capacity = ~0u )
    { return task_create_team( functor , dependence_capacity ); }

  template< class A1 , class A2 , class A3 , class A4 >
  KOKKOS_INLINE_FUNCTION
  void add_dependence( const Future<A1,A2> & after
                     , const Future<A3,A4> & before
                     , typename std::enable_if
                        < std::is_same< typename Future<A1,A2>::execution_space , execution_space >::value
                          &&
                          std::is_same< typename Future<A3,A4>::execution_space , execution_space >::value
                        >::type * = 0
                      ) const
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      m_policy->add_dependence( after.m_task , before.m_task );
#endif
    }

  //----------------------------------------

  Future< Latch , execution_space >
  KOKKOS_INLINE_FUNCTION
  create_latch( const int N ) const
    {
      task_root_type * const task =
        m_policy->allocate_task( sizeof(task_root_type) , 0 , 0 );
      task->m_dep_size = N ; // Using m_dep_size for latch counter
      task->m_state = TASK_STATE_WAITING ;
      return Future< Latch , execution_space >( task );
    }

  //----------------------------------------

  template< class FunctorType , class A3 , class A4 >
  KOKKOS_INLINE_FUNCTION
  void add_dependence( FunctorType * task_functor
                     , const Future<A3,A4> & before
                     , typename std::enable_if
                        < std::is_same< typename Future<A3,A4>::execution_space , execution_space >::value
                        >::type * = 0
                      ) const
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      m_policy->add_dependence( get_task_root(task_functor) , before.m_task );
#endif
    }

  template< class ValueType >
  const Future< ValueType , execution_space > &
    spawn( const Future< ValueType , execution_space > & f
         , const bool priority = false ) const
      {
        if ( f.m_task ) {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
          f.m_task->m_queue =
            ( f.m_task->m_team != 0
            ? & ( m_policy->m_team[   priority ? 0 : 1 ] )
            : & ( m_policy->m_serial[ priority ? 0 : 1 ] ) );
          m_policy->schedule_task( f.m_task );
#endif
        }
        return f ;
      }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  void respawn( FunctorType * task_functor 
              , const bool priority = false ) const
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      task_root_type * const t = get_task_root(task_functor);
      t->m_queue =
        ( t->m_team != 0 ? & ( m_policy->m_team[   priority ? 0 : 1 ] )
                         : & ( m_policy->m_serial[ priority ? 0 : 1 ] ) );
      m_policy->reschedule_task( t );
#endif
    }

  // When a create method fails by returning a null Future
  // the task that called the create method may respawn
  // with a dependence on memory becoming available.
  // This is a race as more than one task may be respawned
  // with this need.

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  void respawn_needing_memory( FunctorType * task_functor ) const
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      task_root_type * const t = get_task_root(task_functor);
      t->m_queue =
        ( t->m_team != 0 ? & ( m_policy->m_team[   2 ] )
                         : & ( m_policy->m_serial[ 2 ] ) );
      m_policy->reschedule_task( t );
#endif
    }

  //----------------------------------------
  // Functions for an executing task functor to query dependences,
  // set new dependences, and respawn itself.

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< void , execution_space >
  get_dependence( const FunctorType * task_functor , int i ) const
    {
      return Future<void,execution_space>(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        get_task_root(task_functor)->get_dependence(i)
#endif
        );
    }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  int get_dependence( const FunctorType * task_functor ) const
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { return get_task_root(task_functor)->get_dependence(); }
#else
    { return 0 ; }
#endif

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  void clear_dependence( FunctorType * task_functor ) const
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { get_task_root(task_functor)->clear_dependence(); }
#else
    {}
#endif

  //----------------------------------------

  static member_type & member_single();

  friend void wait( TaskPolicy< Kokkos::Threads > & );
};

} /* namespace Experimental */
} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_HAVE_PTHREAD ) */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_THREADS_TASKPOLICY_HPP */


