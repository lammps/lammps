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

#ifndef KOKKOS_CUDA_TASKPOLICY_HPP
#define KOKKOS_CUDA_TASKPOLICY_HPP

#include <Kokkos_Core_fwd.hpp>

#if defined( KOKKOS_HAVE_CUDA ) && \
    defined( KOKKOS_CUDA_USE_RELOCATABLE_DEVICE_CODE )

#define KOKKOS_ENABLE_CUDA_TASK_POLICY

/* The TaskPolicy< Cuda > capability requires nvcc using the option:
 *    --relocatable-device-code=true
 */

#include <Kokkos_Cuda.hpp>
#include <Kokkos_TaskPolicy.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct CudaTaskPolicyQueue ;

/** \brief  Base class for all Kokkos::Cuda tasks */
template<>
class TaskMember< Kokkos::Cuda , void , void > {
public:

  template< class > friend class Kokkos::Experimental::TaskPolicy ;
  friend struct CudaTaskPolicyQueue ;

  typedef void (* function_single_type) ( TaskMember * );
  typedef void (* function_team_type)   ( TaskMember * , Kokkos::Impl::CudaTeamMember & );

private:

  friend struct CudaTaskPolicyQueue ;

  CudaTaskPolicyQueue   * m_policy ;
  TaskMember * volatile * m_queue ;
  function_team_type      m_team ;    ///< Apply function on CUDA
  function_single_type    m_serial ;  ///< Apply function on CUDA
  TaskMember **           m_dep ;     ///< Dependences
  TaskMember *            m_wait ;    ///< Linked list of tasks waiting on this task
  TaskMember *            m_next ;    ///< Linked list of tasks waiting on a different task
  int                    m_dep_capacity ; ///< Capacity of dependences
  int                    m_dep_size ;     ///< Actual count of dependences
  int                    m_size_alloc ;
  int                    m_shmem_size ;
  int                    m_ref_count ;    ///< Reference count
  int                    m_state ;        ///< State of the task


  TaskMember( TaskMember && ) = delete ;
  TaskMember( const TaskMember & ) = delete ;
  TaskMember & operator = ( TaskMember && ) = delete ;
  TaskMember & operator = ( const TaskMember & ) = delete ;

protected:

  KOKKOS_INLINE_FUNCTION
  TaskMember()
    : m_policy(0)
    , m_queue(0)
    , m_team(0)
    , m_serial(0)
    , m_dep(0)
    , m_wait(0)
    , m_next(0)
    , m_size_alloc(0)
    , m_dep_capacity(0)
    , m_dep_size(0)
    , m_shmem_size(0)
    , m_ref_count(0)
    , m_state( TASK_STATE_CONSTRUCTING )
    {}

public:

  KOKKOS_FUNCTION
  ~TaskMember();

  KOKKOS_INLINE_FUNCTION
  int reference_count() const
    { return *((volatile int *) & m_ref_count ); }

  // Cannot use the function pointer to verify the type
  // since the function pointer is not unique between
  // Host and Cuda. Don't run verificaton for Cuda. 
  // Assume testing on Host-only back-end will catch such errors.

  template< typename ResultType >
  KOKKOS_INLINE_FUNCTION static
  TaskMember * verify_type( TaskMember * t ) { return t ; }

  //----------------------------------------
  /*  Inheritence Requirements on task types:
   *
   *    class DerivedTaskType
   *      : public TaskMember< Cuda , DerivedType::value_type , FunctorType >
   *      { ... };
   *
   *    class TaskMember< Cuda , DerivedType::value_type , FunctorType >
   *      : public TaskMember< Cuda , DerivedType::value_type , void >
   *      , public Functor
   *      { ... };
   *
   *  If value_type != void
   *    class TaskMember< Cuda , value_type , void >
   *      : public TaskMember< Cuda , void , void >
   *
   *  Allocate space for DerivedTaskType followed by TaskMember*[ dependence_capacity ]
   *
   */
  //----------------------------------------
  // If after the 'apply' the task's state is waiting 
  // then it will be rescheduled and called again.
  // Otherwise the functor must be destroyed.

  template< class DerivedTaskType , class Tag >
  __device__ static
  void apply_single(
    typename std::enable_if
      <( std::is_same< Tag , void >::value &&
        std::is_same< typename DerivedTaskType::result_type , void >::value
       ), TaskMember * >::type t )
    {
      typedef typename DerivedTaskType::functor_type  functor_type ;

      functor_type * const f =
        static_cast< functor_type * >( static_cast< DerivedTaskType * >(t) );

      f->apply();

      if ( t->m_state == int(Kokkos::Experimental::TASK_STATE_EXECUTING) ) {
        f->~functor_type();
      }
    }

  template< class DerivedTaskType , class Tag >
  __device__ static
  void apply_single(
    typename std::enable_if
      <( std::is_same< Tag , void >::value &&
        ! std::is_same< typename DerivedTaskType::result_type , void >::value
       ), TaskMember * >::type t )
    {
      typedef typename DerivedTaskType::functor_type  functor_type ;

      DerivedTaskType * const self = static_cast< DerivedTaskType * >(t);
      functor_type    * const f    = static_cast< functor_type * >( self );

      f->apply( self->m_result );

      if ( t->m_state == int(Kokkos::Experimental::TASK_STATE_EXECUTING) ) {
        f->~functor_type();
      }
    }

  template< class DerivedTaskType , class Tag >
  __device__
  void set_apply_single()
    {
      m_serial = & TaskMember::template apply_single<DerivedTaskType,Tag> ;
    }

  //----------------------------------------

  template< class DerivedTaskType , class Tag >
  __device__ static
  void apply_team(
    typename std::enable_if
      <( std::is_same<Tag,void>::value &&
         std::is_same<typename DerivedTaskType::result_type,void>::value
       ), TaskMember * >::type t
    , Kokkos::Impl::CudaTeamMember & member
    )
    {
      typedef typename DerivedTaskType::functor_type functor_type ;

      functor_type * const f =
        static_cast< functor_type * >( static_cast< DerivedTaskType * >(t) );

      f->apply( member );

      __syncthreads(); // Wait for team to finish calling function

      if ( threadIdx.x == 0 &&
           threadIdx.y == 0 &&
           t->m_state == int(Kokkos::Experimental::TASK_STATE_EXECUTING) ) {
        f->~functor_type();
      }
    }

  template< class DerivedTaskType , class Tag >
  __device__ static
  void apply_team(
    typename std::enable_if
      <( std::is_same<Tag,void>::value &&
         ! std::is_same<typename DerivedTaskType::result_type,void>::value
       ), TaskMember * >::type t
    , Kokkos::Impl::CudaTeamMember & member
    )
    {
      typedef typename DerivedTaskType::functor_type  functor_type ;

      DerivedTaskType * const self = static_cast< DerivedTaskType * >(t);
      functor_type    * const f    = static_cast< functor_type * >( self );

      f->apply( member , self->m_result );

      __syncthreads(); // Wait for team to finish calling function

      if ( threadIdx.x == 0 &&
           threadIdx.y == 0 &&
           t->m_state == int(Kokkos::Experimental::TASK_STATE_EXECUTING) ) {
        f->~functor_type();
      }
    }

  template< class DerivedTaskType , class Tag >
  __device__
  void set_apply_team()
    {
      m_team = & TaskMember::template apply_team<DerivedTaskType,Tag> ;
    }

  //----------------------------------------

  KOKKOS_FUNCTION static
  void assign( TaskMember ** const lhs , TaskMember * const rhs );

  __device__
  TaskMember * get_dependence( int i ) const ;

  __device__
  int get_dependence() const ;

  KOKKOS_FUNCTION void clear_dependence();

  __device__
  void latch_add( const int k );

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION static
  void construct_result( TaskMember * const ) {}

  typedef FutureValueTypeIsVoidError get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const { return get_result_type() ; }

  KOKKOS_INLINE_FUNCTION
  Kokkos::Experimental::TaskState get_state() const { return Kokkos::Experimental::TaskState( m_state ); }

};

/** \brief  A Future< Kokkos::Cuda , ResultType > will cast
 *          from  TaskMember< Kokkos::Cuda , void , void >
 *          to    TaskMember< Kokkos::Cuda , ResultType , void >
 *          to query the result.
 */
template< class ResultType >
class TaskMember< Kokkos::Cuda , ResultType , void >
  : public TaskMember< Kokkos::Cuda , void , void >
{
public:

  typedef ResultType result_type ;

  result_type  m_result ;

  typedef const result_type & get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const { return m_result ; }

  KOKKOS_INLINE_FUNCTION static
  void construct_result( TaskMember * const ptr )
    {
      new((void*)(& ptr->m_result)) result_type();
    }

  TaskMember() = delete ;
  TaskMember( TaskMember && ) = delete ;
  TaskMember( const TaskMember & ) = delete ;
  TaskMember & operator = ( TaskMember && ) = delete ;
  TaskMember & operator = ( const TaskMember & ) = delete ;
};

/** \brief  Callback functions will cast
 *          from  TaskMember< Kokkos::Cuda , void , void >
 *          to    TaskMember< Kokkos::Cuda , ResultType , FunctorType >
 *          to execute work functions.
 */
template< class ResultType , class FunctorType >
class TaskMember< Kokkos::Cuda , ResultType , FunctorType >
  : public TaskMember< Kokkos::Cuda , ResultType , void >
  , public FunctorType
{
public:
  typedef ResultType   result_type ;
  typedef FunctorType  functor_type ;

  KOKKOS_INLINE_FUNCTION static
  void copy_construct( TaskMember * const ptr
                     , const functor_type & arg_functor )
    {
      typedef TaskMember< Kokkos::Cuda , ResultType , void > base_type ;

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

namespace {

template< class DerivedTaskType , class Tag >
__global__
void cuda_set_apply_single( DerivedTaskType * task )
{
  typedef Kokkos::Experimental::Impl::TaskMember< Kokkos::Cuda , void , void >
    task_root_type ;

  task->task_root_type::template set_apply_single< DerivedTaskType , Tag >();
}

template< class DerivedTaskType , class Tag >
__global__
void cuda_set_apply_team( DerivedTaskType * task )
{
  typedef Kokkos::Experimental::Impl::TaskMember< Kokkos::Cuda , void , void >
    task_root_type ;

  task->task_root_type::template set_apply_team< DerivedTaskType , Tag >();
}

} /* namespace */
} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct CudaTaskPolicyQueue {

  enum { NPRIORITY = 3 };

  // Must use UVM so that tasks can be created in both
  // Host and Cuda space.

  typedef Kokkos::Experimental::MemoryPool< Kokkos::CudaUVMSpace >
    memory_space ;

  typedef Kokkos::Experimental::Impl::TaskMember< Kokkos::Cuda , void , void >
    task_root_type ;

  memory_space     m_space ;
  task_root_type * m_team[ NPRIORITY ] ;
  task_root_type * m_serial[ NPRIORITY ];
  int              m_team_size ;
  int              m_default_dependence_capacity ;
  int volatile     m_count_ready ; ///< Ready plus executing tasks

  // Execute tasks until all non-waiting tasks are complete
  __device__
  void driver();

  __device__ static
  task_root_type * pop_ready_task( task_root_type * volatile * const queue );

  // When a task finishes executing.
  __device__
  void complete_executed_task( task_root_type * );

  KOKKOS_FUNCTION void schedule_task( task_root_type * const 
                                    , const bool initial_spawn = true );
  KOKKOS_FUNCTION void reschedule_task( task_root_type * const );
  KOKKOS_FUNCTION
  void add_dependence( task_root_type * const after
                     , task_root_type * const before );


  CudaTaskPolicyQueue() = delete ;
  CudaTaskPolicyQueue( CudaTaskPolicyQueue && ) = delete ;
  CudaTaskPolicyQueue( const CudaTaskPolicyQueue & ) = delete ;
  CudaTaskPolicyQueue & operator = ( CudaTaskPolicyQueue && ) = delete ;
  CudaTaskPolicyQueue & operator = ( const CudaTaskPolicyQueue & ) = delete ;


  ~CudaTaskPolicyQueue();

  // Construct only on the Host
  CudaTaskPolicyQueue
    ( const unsigned arg_task_max_count
    , const unsigned arg_task_max_size
    , const unsigned arg_task_default_dependence_capacity
    , const unsigned arg_task_team_size
    );

  struct Destroy {
    CudaTaskPolicyQueue * m_policy ;
    void destroy_shared_allocation();
  };

  //----------------------------------------
  /** \brief  Allocate and construct a task.
   *
   *  Allocate space for DerivedTaskType followed
   *  by TaskMember*[ dependence_capacity ]
   */
  KOKKOS_FUNCTION
  task_root_type *
  allocate_task( const unsigned arg_sizeof_task
               , const unsigned arg_dep_capacity
               , const unsigned arg_team_shmem = 0 );

  KOKKOS_FUNCTION void deallocate_task( task_root_type * const );
};

} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

void wait( TaskPolicy< Kokkos::Cuda > & );

template<>
class TaskPolicy< Kokkos::Cuda >
{
public:

  typedef Kokkos::Cuda                  execution_space ;
  typedef TaskPolicy                    execution_policy ;
  typedef Kokkos::Impl::CudaTeamMember  member_type ;

private:

  typedef Impl::TaskMember< Kokkos::Cuda , void , void >  task_root_type ;
  typedef Kokkos::Experimental::MemoryPool< Kokkos::CudaUVMSpace > memory_space ;
  typedef Kokkos::Experimental::Impl::SharedAllocationTracker track_type ;

  track_type                   m_track ;
  Impl::CudaTaskPolicyQueue  * m_policy ;

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION static
  const task_root_type * get_task_root( const FunctorType * f )
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType > task_type ;
      return static_cast< const task_root_type * >( static_cast< const task_type * >(f) );
    }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION static
  task_root_type * get_task_root( FunctorType * f )
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType > task_type ;
      return static_cast< task_root_type * >( static_cast< task_type * >(f) );
    }

public:

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

  KOKKOS_FUNCTION
  int allocated_task_count() const { return 0 ; }

  //----------------------------------------
  // Create serial-thread task
  // Main process and tasks must use different functions
  // to work around CUDA limitation where __host__ __device__
  // functions are not allowed to invoke templated __global__ functions.

  template< class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  proc_create( const FunctorType & arg_functor
             , const unsigned      arg_dep_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type  value_type ;

      typedef Impl::TaskMember< execution_space , value_type , FunctorType >
        task_type ;

      task_type * const task =
        static_cast<task_type*>(
          m_policy->allocate_task( sizeof(task_type) , arg_dep_capacity ) );

      if ( task ) {
        // The root part of the class has been constructed.
        // Must now construct the functor and result specific part.

        task_type::copy_construct( task , arg_functor );

        // Setting the apply pointer on the device requires code
        // executing on the GPU.  This function is called on the
        // host process so a kernel must be run.

        // Launching a kernel will cause the allocated task in
        // UVM memory to be copied to the GPU.
        // Synchronize to guarantee non-concurrent access
        // between host and device.

        CUDA_SAFE_CALL( cudaDeviceSynchronize() );

        Impl::cuda_set_apply_single<task_type,void><<<1,1>>>( task );

        CUDA_SAFE_CALL( cudaGetLastError() );
        CUDA_SAFE_CALL( cudaDeviceSynchronize() );
      }

      return Future< value_type , execution_space >( task );
    }

  template< class FunctorType >
  __device__
  Future< typename FunctorType::value_type , execution_space >
  task_create( const FunctorType & arg_functor
             , const unsigned      arg_dep_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type  value_type ;

      typedef Impl::TaskMember< execution_space , value_type , FunctorType >
        task_type ;

      task_type * const task =
        static_cast<task_type*>(
          m_policy->allocate_task( sizeof(task_type) , arg_dep_capacity ) );

      if ( task ) {
        // The root part of the class has been constructed.
        // Must now construct the functor and result specific part.

        task_type::copy_construct( task , arg_functor );

        // Setting the apply pointer on the device requires code
        // executing on the GPU.  If this function is called on the
        // Host then a kernel must be run.

        task->task_root_type::template set_apply_single< task_type , void >();
      }

      return Future< value_type , execution_space >( task );
    }

  //----------------------------------------
  // Create thread-team task
  // Main process and tasks must use different functions
  // to work around CUDA limitation where __host__ __device__
  // functions are not allowed to invoke templated __global__ functions.

  template< class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  proc_create_team( const FunctorType & arg_functor
                  , const unsigned      arg_dep_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type  value_type ;

      typedef Impl::TaskMember< execution_space , value_type , FunctorType >
        task_type ;

      const unsigned team_shmem_size =
        Kokkos::Impl::FunctorTeamShmemSize< FunctorType >::value
           ( arg_functor , m_policy->m_team_size );

      task_type * const task =
        static_cast<task_type*>(
          m_policy->allocate_task( sizeof(task_type) , arg_dep_capacity , team_shmem_size ) );

      if ( task ) {
        // The root part of the class has been constructed.
        // Must now construct the functor and result specific part.

        task_type::copy_construct( task , arg_functor );

        // Setting the apply pointer on the device requires code
        // executing on the GPU.  This function is called on the
        // host process so a kernel must be run.

        // Launching a kernel will cause the allocated task in
        // UVM memory to be copied to the GPU.
        // Synchronize to guarantee non-concurrent access
        // between host and device.

        CUDA_SAFE_CALL( cudaDeviceSynchronize() );

        Impl::cuda_set_apply_team<task_type,void><<<1,1>>>( task );

        CUDA_SAFE_CALL( cudaGetLastError() );
        CUDA_SAFE_CALL( cudaDeviceSynchronize() );
      }

      return Future< value_type , execution_space >( task );
    }

  template< class FunctorType >
  __device__
  Future< typename FunctorType::value_type , execution_space >
  task_create_team( const FunctorType & arg_functor
                  , const unsigned      arg_dep_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type  value_type ;

      typedef Impl::TaskMember< execution_space , value_type , FunctorType >
        task_type ;

      const unsigned team_shmem_size =
        Kokkos::Impl::FunctorTeamShmemSize< FunctorType >::value
           ( arg_functor , m_policy->m_team_size );

      task_type * const task =
        static_cast<task_type*>(
          m_policy->allocate_task( sizeof(task_type) , arg_dep_capacity , team_shmem_size ) );

      if ( task ) {
        // The root part of the class has been constructed.
        // Must now construct the functor and result specific part.

        task_type::copy_construct( task , arg_functor );

        // Setting the apply pointer on the device requires code
        // executing on the GPU.  If this function is called on the
        // Host then a kernel must be run.

        task->task_root_type::template set_apply_team< task_type , void >();
      }

      return Future< value_type , execution_space >( task );
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
    { m_policy->add_dependence( after.m_task , before.m_task ); }

  template< class FunctorType , class A3 , class A4 >
  KOKKOS_INLINE_FUNCTION
  void add_dependence( FunctorType * task_functor
                     , const Future<A3,A4> & before
                     , typename std::enable_if
                        < std::is_same< typename Future<A3,A4>::execution_space , execution_space >::value
                        >::type * = 0
                      ) const
    { m_policy->add_dependence( get_task_root(task_functor) , before.m_task ); }


  template< class ValueType >
  KOKKOS_INLINE_FUNCTION
  const Future< ValueType , execution_space > &
    spawn( const Future< ValueType , execution_space > & f 
         , const bool priority = false ) const
      {
        if ( f.m_task ) {
          f.m_task->m_queue =
            ( f.m_task->m_team != 0
            ? & ( m_policy->m_team[   priority ? 0 : 1 ] )
            : & ( m_policy->m_serial[ priority ? 0 : 1 ] ) );
          m_policy->schedule_task( f.m_task );
        }
        return f ;
      }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  void respawn( FunctorType * task_functor 
              , const bool priority = false ) const
    {
      task_root_type * const t = get_task_root(task_functor);
      t->m_queue =
        ( t->m_team != 0 ? & ( m_policy->m_team[   priority ? 0 : 1 ] )
                         : & ( m_policy->m_serial[ priority ? 0 : 1 ] ) );
      m_policy->reschedule_task( t );
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
      task_root_type * const t = get_task_root(task_functor);
      t->m_queue =
        ( t->m_team != 0 ? & ( m_policy->m_team[   2 ] )
                         : & ( m_policy->m_serial[ 2 ] ) );
      m_policy->reschedule_task( t );
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
        get_task_root(task_functor)->get_dependence(i)
      );
    }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  int get_dependence( const FunctorType * task_functor ) const
    { return get_task_root(task_functor)->get_dependence(); }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  void clear_dependence( FunctorType * task_functor ) const
    { get_task_root(task_functor)->clear_dependence(); }

  //----------------------------------------

  __device__
  static member_type member_single()
    {
      return
        member_type( 0 /* shared memory */
                   , 0 /* shared memory begin */
                   , 0 /* shared memory size */
                   , 0 /* league rank */
                   , 1 /* league size */ );
    }

  friend void wait( TaskPolicy< Kokkos::Cuda > & );
};

} /* namespace Experimental */
} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_HAVE_CUDA ) && defined( KOKKOS_CUDA_USE_RELOCATABLE_DEVICE_CODE ) */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CUDA_TASKPOLICY_HPP */


