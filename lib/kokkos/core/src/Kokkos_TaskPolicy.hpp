
/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Tags.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct FutureValueTypeIsVoidError {};

template< class ExecSpace >
class TaskManager ;

template < class Policy , class ResultType = void , class Functor = void >
class TaskMember ;

template< class ExecPolicy >
struct TaskDepends { typedef typename ExecPolicy::execution_space  execution_space ; };

template< class ExecPolicy >
struct TaskForEach { typedef typename ExecPolicy::execution_space  execution_space ; };

template< class ExecPolicy >
struct TaskReduce { typedef typename ExecPolicy::execution_space  execution_space ; };

template< class ExecPolicy >
struct TaskScan { typedef typename ExecPolicy::execution_space  execution_space ; };

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {

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

  template< class > friend class Impl::TaskManager ;
  template< class , class > friend class Future ;

  // Argument #2, if not void, must be the space.
  enum { Arg1_is_space  = Impl::is_execution_space< Arg1 >::value };
  enum { Arg2_is_space  = Impl::is_execution_space< Arg2 >::value };
  enum { Arg2_is_void   = Impl::is_same< Arg2 , void >::value };

  struct ErrorNoExecutionSpace {};

  enum { Opt1  =   Arg1_is_space && Arg2_is_void
       , Opt2  = ! Arg1_is_space && Arg2_is_void
       , Opt3  = ! Arg1_is_space && Arg2_is_space 
       , OptOK = Impl::StaticAssert< Opt1 || Opt2 || Opt3 , ErrorNoExecutionSpace >::value
       };

  typedef typename
    Impl::if_c< Opt2 || Opt3 , Arg1 , void >::type
      ValueType ;

  typedef typename
    Impl::if_c< Opt1 , Arg1 , typename
    Impl::if_c< Opt2 , Kokkos::DefaultExecutionSpace , typename
    Impl::if_c< Opt3 , Arg2 , void
    >::type >::type >::type
      ExecutionSpace ;

  typedef Impl::TaskManager< ExecutionSpace >              TaskManager ;
  typedef Impl::TaskMember<  ExecutionSpace >              TaskRoot ;
  typedef Impl::TaskMember<  ExecutionSpace , ValueType >  TaskValue ;

  TaskRoot * m_task ;

public:

  typedef ValueType       value_type;
  typedef ExecutionSpace  execution_space ;

  //----------------------------------------

  explicit
  Future( TaskRoot * task )
    : m_task(0)
    { TaskManager::assign( & m_task , TaskValue::verify_type( task ) ); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  ~Future() { TaskManager::assign( & m_task , 0 ); }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  Future() : m_task(0) {}

  KOKKOS_INLINE_FUNCTION
  Future( const Future & rhs )
    : m_task(0)
    { TaskManager::assign( & m_task , rhs.m_task ); }

  KOKKOS_INLINE_FUNCTION
  Future & operator = ( const Future & rhs )
    { TaskManager::assign( & m_task , rhs.m_task ); return *this ; }

  //----------------------------------------

  template< class A1 , class A2 >
  KOKKOS_INLINE_FUNCTION
  Future( const Future<A1,A2> & rhs )
    : m_task(0)
    { TaskManager::assign( & m_task , TaskValue::verify_type( rhs.m_task ) ); }

  template< class A1 , class A2 >
  KOKKOS_INLINE_FUNCTION
  Future & operator = ( const Future<A1,A2> & rhs )
    { TaskManager::assign( & m_task , TaskValue::verify_type( rhs.m_task ) ); return *this ; }

  //----------------------------------------

  typedef typename TaskValue::get_result_type get_result_type ;

  KOKKOS_INLINE_FUNCTION
  typename TaskValue::get_result_type get() const
    { return static_cast<TaskValue*>( m_task )->get(); }
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {

template< class Policy = Kokkos::DefaultExecutionSpace >
class TaskPolicy {
public:

  typedef typename Policy::execution_space  execution_space ;

  template< class A1 , class A2 >
  void wait( const Future<A1,A2> & ) const ;

  template< class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  spawn( const FunctorType & ) const ;

  template< class FunctorType >
  void respawn( FunctorType * ) const ;

  template< class FunctorType >
  Future< void , execution_space >
  get_dependence( FunctorType * ) const ;

  template< class ValueType >
  TaskPolicy< void /* ... */ >
  depends( const Future< ValueType , execution_space > * const , const int );

  template< class ExecPolicy >
  TaskPolicy< void /* ... */ > foreach( const ExecPolicy & );

  template< class ExecPolicy >
  TaskPolicy< void /* ... */ > reduce( const ExecPolicy & );

  template< class ExecPolicy >
  TaskPolicy< void /* ... */ > scan( const ExecPolicy & );
};

// spawn( M.depends(n,d).foreach(K) , functor );
// M.depends(n,d).foreach(K).spawn( functor );

template< class PolicyType , class FunctorType >
Future< typename FunctorType::value_type
      , typename PolicyType::execution_space >
inline
spawn( const TaskPolicy< PolicyType > & policy
     , const FunctorType              & functor )
{ return policy.spawn( functor ); }

template< class PolicyType , class A1 , class A2 >
void wait( const TaskPolicy< PolicyType > & policy 
         , const Future<A1,A2>            & future
         , typename Impl::enable_if<
             Impl::is_same< typename PolicyType::execution_space
                          , typename Future<A1,A2>::execution_space >::value
          >::type * = 0 )
{ policy.wait( future ); }

} /* namespace Kokkos */

//----------------------------------------------------------------------------

#endif /* #define KOKKOS_TASKPOLICY_HPP */

