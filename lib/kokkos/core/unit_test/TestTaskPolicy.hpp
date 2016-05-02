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


#ifndef KOKKOS_UNITTEST_TASKPOLICY_HPP
#define KOKKOS_UNITTEST_TASKPOLICY_HPP

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <Kokkos_TaskPolicy.hpp>

namespace TestTaskPolicy {

//----------------------------------------------------------------------------

template< class ExecSpace >
struct FibChild {

  typedef long value_type ;

  Kokkos::Experimental::TaskPolicy<ExecSpace> policy ;
  Kokkos::Experimental::Future<long,ExecSpace> fib_1 ;
  Kokkos::Experimental::Future<long,ExecSpace> fib_2 ;
  const value_type n ;
  int has_nested ;

  KOKKOS_INLINE_FUNCTION
  FibChild( const Kokkos::Experimental::TaskPolicy<ExecSpace> & arg_policy
          , const value_type arg_n )
    : policy(arg_policy)
    , fib_1() , fib_2()
    , n( arg_n ), has_nested(0) {}

  KOKKOS_INLINE_FUNCTION
  void apply( value_type & result )
    {
      typedef Kokkos::Experimental::Future<long,ExecSpace> future_type ;

      if ( n < 2 ) {

        has_nested = -1 ;

        result = n ;
      }
      else {
        if ( has_nested == 0 ) {
          // Spawn new children and respawn myself to sum their results:
          // Spawn lower value at higher priority as it has a shorter
          // path to completion.
          if ( fib_2.is_null() ) {
            fib_2 = policy.task_create( FibChild(policy,n-2) );
          }

          if ( ! fib_2.is_null() && fib_1.is_null() ) {
            fib_1 = policy.task_create( FibChild(policy,n-1) );
          }

          if ( ! fib_1.is_null() ) {
            has_nested = 2 ;

            policy.spawn( fib_2 , true /* high priority */ );
            policy.spawn( fib_1 );
            policy.add_dependence( this , fib_1 );
            policy.add_dependence( this , fib_2 );
            policy.respawn( this );
          }
          else {
            // Release task memory before spawning the task,
            // after spawning memory cannot be released.
            fib_2 = future_type();
            // Respawn when more memory is available
            policy.respawn_needing_memory( this );
          }
        }
        else if ( has_nested == 2 ) {

          has_nested = -1 ;

          result = fib_1.get() + fib_2.get();

if ( false ) {
  printf("FibChild %ld = fib(%ld), task_count(%d)\n"
        , long(n), long(result), policy.allocated_task_count());
}

        }
        else {
          printf("FibChild(%ld) execution error\n",(long)n);
          Kokkos::abort("FibChild execution error");
        }
      }
    }
};

template< class ExecSpace >
struct FibChild2 {

  typedef long value_type ;

  Kokkos::Experimental::TaskPolicy<ExecSpace> policy ;
  Kokkos::Experimental::Future<long,ExecSpace> fib_a ;
  Kokkos::Experimental::Future<long,ExecSpace> fib_b ;
  const value_type n ;
  int has_nested ;

  KOKKOS_INLINE_FUNCTION
  FibChild2( const Kokkos::Experimental::TaskPolicy<ExecSpace> & arg_policy
           , const value_type arg_n )
    : policy(arg_policy)
    , n( arg_n ), has_nested(0) {}

  KOKKOS_INLINE_FUNCTION
  void apply( value_type & result )
    {
      if ( 0 == has_nested ) {
        if ( n < 2 ) {

          has_nested = -1 ;

          result = n ;
        }
        else if ( n < 4 ) {
          // Spawn new children and respawn myself to sum their results:
          // result = Fib(n-1) + Fib(n-2)
          has_nested = 2 ;

          // Spawn lower value at higher priority as it has a shorter
          // path to completion.

          policy.clear_dependence( this );
          fib_a = policy.spawn( policy.task_create( FibChild2(policy,n-1) ) );
          fib_b = policy.spawn( policy.task_create( FibChild2(policy,n-2) ) , true );
          policy.add_dependence( this , fib_a );
          policy.add_dependence( this , fib_b );
          policy.respawn( this );
        }
        else {
          // Spawn new children and respawn myself to sum their results:
          // result = Fib(n-1) + Fib(n-2)
          // result = ( Fib(n-2) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = ( ( Fib(n-3) + Fib(n-4) ) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = 3 * Fib(n-3) + 2 * Fib(n-4)
          has_nested = 4 ;

          // Spawn lower value at higher priority as it has a shorter
          // path to completion.

          policy.clear_dependence( this );
          fib_a = policy.spawn( policy.task_create( FibChild2(policy,n-3) ) );
          fib_b = policy.spawn( policy.task_create( FibChild2(policy,n-4) ) , true );
          policy.add_dependence( this , fib_a );
          policy.add_dependence( this , fib_b );
          policy.respawn( this );
        }
     }
     else if ( 2 == has_nested || 4 == has_nested ) {
        result = ( has_nested == 2 ) ? fib_a.get() + fib_b.get()
                                     : 3 * fib_a.get() + 2 * fib_b.get() ;

        has_nested = -1 ;
      }
      else {
        printf("FibChild2(%ld) execution error\n",(long)n);
        Kokkos::abort("FibChild2 execution error");
      }
    }
};

namespace {

long eval_fib( long n )
{
  if ( 2 <= n ) {
    std::vector<long> fib(n+1);

    fib[0] = 0 ;
    fib[1] = 1 ;

    for ( long i = 2 ; i <= n ; ++i ) { fib[i] = fib[i-2] + fib[i-1]; }

    n = fib[n] ;
  }

  return n ;
}

}

template< class ExecSpace >
void test_fib( long n , const unsigned task_max_count = 1024 )
{
  const unsigned task_max_size   = 256 ;
  const unsigned task_dependence = 4 ;

  Kokkos::Experimental::TaskPolicy<ExecSpace>
    policy( task_max_count
          , task_max_size
          , task_dependence );

  Kokkos::Experimental::Future<long,ExecSpace> f =
    policy.spawn( policy.proc_create( FibChild<ExecSpace>(policy,n) ) );

  Kokkos::Experimental::wait( policy );

  if ( f.get() != eval_fib(n) ) {
    std::cout << "Fib(" << n << ") = " << f.get();
    std::cout << " != " << eval_fib(n);
    std::cout << std::endl ;
  }
}

template< class ExecSpace >
void test_fib2( long n , const unsigned task_max_count = 1024 )
{
  const unsigned task_max_size   = 256 ;
  const unsigned task_dependence = 4 ;

  Kokkos::Experimental::TaskPolicy<ExecSpace>
    policy( task_max_count
          , task_max_size
          , task_dependence );

  Kokkos::Experimental::Future<long,ExecSpace> f =
    policy.spawn( policy.proc_create( FibChild2<ExecSpace>(policy,n) ) );

  Kokkos::Experimental::wait( policy );

  if ( f.get() != eval_fib(n) ) {
    std::cout << "Fib2(" << n << ") = " << f.get();
    std::cout << " != " << eval_fib(n);
    std::cout << std::endl ;
  }
}

//----------------------------------------------------------------------------

template< class ExecSpace >
struct Norm2 {

  typedef double value_type ;

  const double * const m_x ;

  Norm2( const double * x ) : m_x(x) {}

  inline
  void init( double & val ) const { val = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void operator()( int i , double & val ) const { val += m_x[i] * m_x[i] ; }

  void apply( double & dst ) const { dst = std::sqrt( dst ); }
};

template< class ExecSpace >
void test_norm2( const int n )
{
  const unsigned task_max_count  = 1024 ;
  const unsigned task_max_size   = 256 ;
  const unsigned task_dependence = 4 ;

  Kokkos::Experimental::TaskPolicy<ExecSpace>
    policy( task_max_count
          , task_max_size
          , task_dependence );

  double * const x = new double[n];

  for ( int i = 0 ; i < n ; ++i ) x[i] = 1 ;

  Kokkos::RangePolicy<ExecSpace> r(0,n);

  Kokkos::Experimental::Future<double,ExecSpace> f =
    Kokkos::Experimental::spawn_reduce( policy , r , Norm2<ExecSpace>(x) );

  Kokkos::Experimental::wait( policy );

#if defined(PRINT)
  std::cout << "Norm2: " << f.get() << std::endl ;
#endif

  delete[] x ;
}

//----------------------------------------------------------------------------

template< class Space >
struct TaskDep {

  typedef int value_type ;
  typedef Kokkos::Experimental::TaskPolicy< Space > policy_type ;

  const policy_type policy ;
  const int         input ;

  TaskDep( const policy_type & arg_p , const int arg_i )
    : policy( arg_p ), input( arg_i ) {}

  KOKKOS_INLINE_FUNCTION
  void apply( int & val )
  {
    val = input ;
    const int num = policy.get_dependence( this );

    for ( int i = 0 ; i < num ; ++i ) {
      Kokkos::Experimental::Future<int,Space> f = policy.get_dependence( this , i );
      val += f.get();
    }
  }
};


template< class Space >
void test_task_dep( const int n )
{
  enum { NTEST = 64 };

  const unsigned task_max_count  = 1024 ;
  const unsigned task_max_size   = 64 ;
  const unsigned task_dependence = 4 ;

  Kokkos::Experimental::TaskPolicy<Space>
    policy( task_max_count
          , task_max_size
          , task_dependence );

  Kokkos::Experimental::Future<int,Space> f[ NTEST ];

  for ( int i = 0 ; i < NTEST ; ++i ) {
    // Create task in the "constructing" state with capacity for 'n+1' dependences
    f[i] = policy.proc_create( TaskDep<Space>(policy,0) , n + 1 );

    if ( f[i].get_task_state() != Kokkos::Experimental::TASK_STATE_CONSTRUCTING ) {
      Kokkos::Impl::throw_runtime_exception("get_task_state() != Kokkos::Experimental::TASK_STATE_CONSTRUCTING");
    }

    // Only use 'n' dependences

    for ( int j = 0 ; j < n ; ++j ) {

      Kokkos::Experimental::Future<int,Space> nested =
        policy.proc_create( TaskDep<Space>(policy,j+1) );

      policy.spawn( nested );

      // Add dependence to a "constructing" task
      policy.add_dependence( f[i] , nested );
    }

    // Spawn task from the "constructing" to the "waiting" state
    policy.spawn( f[i] );
  }

  const int answer = n % 2 ? n * ( ( n + 1 ) / 2 ) : ( n / 2 ) * ( n + 1 );

  Kokkos::Experimental::wait( policy );

  int error = 0 ;
  for ( int i = 0 ; i < NTEST ; ++i ) {
    if ( f[i].get_task_state() != Kokkos::Experimental::TASK_STATE_COMPLETE ) {
      Kokkos::Impl::throw_runtime_exception("get_task_state() != Kokkos::Experimental::TASK_STATE_COMPLETE");
    }
    if ( answer != f[i].get() && 0 == error ) {
      std::cout << "test_task_dep(" << n << ") ERROR at[" << i << "]"
                << " answer(" << answer << ") != result(" << f[i].get() << ")" << std::endl ;
    }
  }
}

//----------------------------------------------------------------------------

template< class ExecSpace >
struct TaskTeam {

  enum { SPAN = 8 };

  typedef void value_type ;
  typedef Kokkos::Experimental::TaskPolicy<ExecSpace>  policy_type ;
  typedef Kokkos::Experimental::Future<void,ExecSpace> future_type ;
  typedef Kokkos::View<long*,ExecSpace>                view_type ;

  policy_type  policy ;
  future_type  future ;

  view_type  result ;
  const long nvalue ;

  KOKKOS_INLINE_FUNCTION
  TaskTeam( const policy_type & arg_policy
          , const view_type   & arg_result
          , const long          arg_nvalue )
    : policy(arg_policy)
    , future()
    , result( arg_result )
    , nvalue( arg_nvalue )
    {}

  KOKKOS_INLINE_FUNCTION
  void apply( const typename policy_type::member_type & member )
    {
      const long end   = nvalue + 1 ;
      const long begin = 0 < end - SPAN ? end - SPAN : 0 ;

      if ( 0 < begin && future.get_task_state() == Kokkos::Experimental::TASK_STATE_NULL ) {
        if ( member.team_rank() == 0 ) {
          future = policy.spawn( policy.task_create_team( TaskTeam( policy , result , begin - 1 ) ) );
          policy.clear_dependence( this );
          policy.add_dependence( this , future );
          policy.respawn( this );
        }
        return ;
      }

      Kokkos::parallel_for( Kokkos::TeamThreadRange(member,begin,end)
                          , [&]( int i ) { result[i] = i + 1 ; }
                          );
    }
};

template< class ExecSpace >
struct TaskTeamValue {

  enum { SPAN = 8 };

  typedef long value_type ;
  typedef Kokkos::Experimental::TaskPolicy<ExecSpace>         policy_type ;
  typedef Kokkos::Experimental::Future<value_type,ExecSpace>  future_type ;
  typedef Kokkos::View<long*,ExecSpace>                       view_type ;

  policy_type  policy ;
  future_type  future ;

  view_type  result ;
  const long nvalue ;

  KOKKOS_INLINE_FUNCTION
  TaskTeamValue( const policy_type & arg_policy
               , const view_type   & arg_result
               , const long          arg_nvalue )
    : policy(arg_policy)
    , future()
    , result( arg_result )
    , nvalue( arg_nvalue )
    {}

  KOKKOS_INLINE_FUNCTION
  void apply( const typename policy_type::member_type & member , value_type & final )
    {
      const long end   = nvalue + 1 ;
      const long begin = 0 < end - SPAN ? end - SPAN : 0 ;

      if ( 0 < begin && future.is_null() ) {
        if ( member.team_rank() == 0 ) {

          future = policy.task_create_team( TaskTeamValue( policy , result , begin - 1 ) );

          policy.spawn( future );
          policy.add_dependence( this , future );
          policy.respawn( this );
        }
        return ;
      }

      Kokkos::parallel_for( Kokkos::TeamThreadRange(member,begin,end)
                          , [&]( int i ) { result[i] = i + 1 ; }
                          );

      if ( member.team_rank() == 0 ) {
        final = result[nvalue] ;
      }

      Kokkos::memory_fence();
    }
};

template< class ExecSpace >
void test_task_team( long n )
{
  typedef TaskTeam< ExecSpace >            task_type ;
  typedef TaskTeamValue< ExecSpace >       task_value_type ;
  typedef typename task_type::view_type    view_type ;
  typedef typename task_type::policy_type  policy_type ;

  typedef typename task_type::future_type        future_type ;
  typedef typename task_value_type::future_type  future_value_type ;

  const unsigned task_max_count  = 1024 ;
  const unsigned task_max_size   = 256 ;
  const unsigned task_dependence = 4 ;

  policy_type
    policy( task_max_count
          , task_max_size
          , task_dependence );

  view_type    result("result",n+1);

  typename view_type::HostMirror
    host_result = Kokkos::create_mirror_view( result );

  future_type f = policy.proc_create_team( task_type( policy , result , n ) );

  ASSERT_FALSE( f.is_null() );

  policy.spawn( f );

  Kokkos::Experimental::wait( policy );

  Kokkos::deep_copy( host_result , result );

  for ( long i = 0 ; i <= n ; ++i ) {
    const long answer = i + 1 ;
    if ( host_result(i) != answer ) {
      std::cerr << "test_task_team void ERROR result(" << i << ") = "
                << host_result(i) << " != " << answer << std::endl ;
    }
  }

  future_value_type fv = policy.proc_create_team( task_value_type( policy , result , n ) );

  ASSERT_FALSE( fv.is_null() );

  policy.spawn( fv );

  Kokkos::Experimental::wait( policy );

  Kokkos::deep_copy( host_result , result );

  if ( fv.get() != n + 1 ) {
    std::cerr << "test_task_team value ERROR future = "
              << fv.get() << " != " << n + 1 << std::endl ;
  }
  for ( long i = 0 ; i <= n ; ++i ) {
    const long answer = i + 1 ;
    if ( host_result(i) != answer ) {
      std::cerr << "test_task_team value ERROR result(" << i << ") = "
                << host_result(i) << " != " << answer << std::endl ;
    }
  }
}

//----------------------------------------------------------------------------

template< class ExecSpace >
struct TaskLatchAdd {

  typedef void value_type ;
  typedef Kokkos::Experimental::Future< Kokkos::Experimental::Latch , ExecSpace >  future_type ;

  future_type     latch ;
  volatile int *  count ;

  KOKKOS_INLINE_FUNCTION
  TaskLatchAdd( const future_type & arg_latch 
              , volatile int * const arg_count )
    : latch( arg_latch )
    , count( arg_count )
    {}

  KOKKOS_INLINE_FUNCTION
  void apply()
    {
      Kokkos::atomic_fetch_add( count , 1 );
      latch.add(1);
    }
};

template< class ExecSpace >
struct TaskLatchRun {

  typedef void value_type ;
  typedef Kokkos::Experimental::TaskPolicy< ExecSpace >      policy_type ;
  typedef Kokkos::Experimental::Future< Kokkos::Experimental::Latch , ExecSpace >  future_type ;

  policy_type policy ;
  int total ;
  volatile int count ;

  KOKKOS_INLINE_FUNCTION
  TaskLatchRun( const policy_type & arg_policy , const int arg_total )
    : policy(arg_policy), total(arg_total), count(0) {}

  KOKKOS_INLINE_FUNCTION
  void apply()
    {
      if ( 0 == count && 0 < total ) {
        future_type latch = policy.create_latch( total );

        for ( int i = 0 ; i < total ; ++i ) {
          auto f = policy.task_create( TaskLatchAdd<ExecSpace>(latch,&count) , 0 );
          if ( f.is_null() ) {
            Kokkos::abort("TaskLatchAdd allocation FAILED" );
          }

          if ( policy.spawn( f ).is_null() ) {
            Kokkos::abort("TaskLatcAdd spawning FAILED" );
          }
        }

        policy.add_dependence( this , latch );
        policy.respawn( this );
      }
      else if ( count != total ) {
        printf("TaskLatchRun FAILED %d != %d\n",count,total);
      }
    }
};


template< class ExecSpace >
void test_latch( int n )
{
  typedef TaskLatchRun< ExecSpace >        task_type ;
  typedef typename task_type::policy_type  policy_type ;

  // Primary + latch + n*LatchAdd
  const unsigned task_max_count  = n + 2 ;
  const unsigned task_max_size   = sizeof(task_type);
  const unsigned task_dependence = 4 ;

  policy_type
    policy( task_max_count
          , task_max_size
          , task_dependence );


  policy.spawn( policy.proc_create( TaskLatchRun<ExecSpace>(policy,n) ) );

  wait( policy );
}

//----------------------------------------------------------------------------

} // namespace TestTaskPolicy

#endif /* #ifndef KOKKOS_UNITTEST_TASKPOLICY_HPP */


