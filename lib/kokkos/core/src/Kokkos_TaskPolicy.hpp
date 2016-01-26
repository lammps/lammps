
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

#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Tags.hpp>
#include <impl/Kokkos_StaticAssert.hpp>
#include <impl/Kokkos_AllocationTracker.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

struct FutureValueTypeIsVoidError {};

template < class ExecSpace , class ResultType , class FunctorType >
class TaskMember ;

template< class ExecPolicy , class ResultType , class FunctorType >
class TaskForEach ;

template< class ExecPolicy , class ResultType , class FunctorType >
class TaskReduce ;

template< class ExecPolicy , class ResultType , class FunctorType >
struct TaskScan ;

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
  enum { Arg2_is_void   = Kokkos::Impl::is_same< Arg2 , void >::value };

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

public:

  typedef ValueType       value_type;
  typedef ExecutionSpace  execution_space ;

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  TaskState get_task_state() const
    { return 0 != m_task ? m_task->get_state() : TASK_STATE_NULL ; }

  //----------------------------------------

  explicit
  Future( TaskRoot * task )
    : m_task(0)
    { TaskRoot::assign( & m_task , TaskRoot::template verify_type< value_type >( task ) ); }

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
    { return static_cast<TaskValue*>( m_task )->get(); }
};

namespace Impl {

template< class T >
struct is_future : public Kokkos::Impl::bool_< false > {};

template< class Arg0 , class Arg1 >
struct is_future< Kokkos::Experimental::Future<Arg0,Arg1> > : public Kokkos::Impl::bool_< true > {};

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
  /** \brief  Create a serial task with storage for dependences.
   *
   *  Postcondition: Task is in the 'constructing' state.
   */
  template< class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  create( const FunctorType & functor
        , const unsigned      dependence_capacity /* = default */ ) const ;

  /** \brief  Create a foreach task with storage for dependences. */
  template< class ExecPolicy , class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  create_foreach( const ExecPolicy  & policy
                , const FunctorType & functor
                , const unsigned      dependence_capacity /* = default */ ) const ;

  /** \brief  Create a reduce task with storage for dependences. */
  template< class ExecPolicy , class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  create_reduce( const ExecPolicy  & policy
               , const FunctorType & functor
               , const unsigned      dependence_capacity /* = default */ ) const ;

  /** \brief  Create a scan task with storage for dependences. */
  template< class ExecPolicy , class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  create_scan( const ExecPolicy  & policy
             , const FunctorType & functor
             , const unsigned      dependence_capacity /* = default */ ) const ;

  /** \brief  Set dependence that 'after' cannot start execution
   *          until 'before' has completed.
   *
   *  Precondition: The 'after' task must be in then 'Constructing' state.
   */
  template< class TA , class TB >
  void set_dependence( const Future<TA,execution_space> & after
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

  /** \brief  Set dependence that 'after' cannot start execution
   *          until 'before' has completed.
   *
   *  The 'after' functor must be in the executing state
   */
  template< class FunctorType , class TB >
  void set_dependence( FunctorType * after
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

#endif /* #define KOKKOS_TASKPOLICY_HPP */

