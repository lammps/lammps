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
  const value_type n ;
  int has_nested ;

  FibChild( const Kokkos::Experimental::TaskPolicy<ExecSpace> & arg_policy
          , const value_type arg_n )
    : policy(arg_policy,2) /* default dependence capacity = 2 */
    , n( arg_n ), has_nested(0) {}

  inline
  void apply( value_type & result )
    {
      if ( n < 2 ) {

        has_nested = -1 ;

        result = n ;
      }
      else {
        if ( has_nested == 0 ) {
          // Spawn new children and respawn myself to sum their results:
          has_nested = 2 ;

          Kokkos::Experimental::respawn
            ( policy
            , this
            , Kokkos::Experimental::spawn( policy , FibChild(policy,n-1) )
            , Kokkos::Experimental::spawn( policy , FibChild(policy,n-2) )
            );

        }
        else if ( has_nested == 2 ) {

          has_nested = -1 ;

          const Kokkos::Experimental::Future<long,ExecSpace> fib_1 = policy.get_dependence(this,0);
          const Kokkos::Experimental::Future<long,ExecSpace> fib_2 = policy.get_dependence(this,1);

          result = fib_1.get() + fib_2.get();
        }
        else {
          fprintf(stderr,"FibChild(%ld) execution error\n",(long)n);
          fflush(stderr);
        }
      }
    }
};

template< class ExecSpace >
struct FibChild2 {

  typedef long value_type ;

  Kokkos::Experimental::TaskPolicy<ExecSpace> policy ;
  const value_type n ;
  int has_nested ;

  FibChild2( const Kokkos::Experimental::TaskPolicy<ExecSpace> & arg_policy
           , const value_type arg_n )
    : policy(arg_policy,2) /* default dependence capacity = 2 */
    , n( arg_n ), has_nested(0) {}

  inline
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
          // Kokkos::respawn implements the following steps:
          policy.clear_dependence( this );
          policy.add_dependence( this , Kokkos::Experimental::spawn( policy , FibChild2(policy,n-1) ) );
          policy.add_dependence( this , Kokkos::Experimental::spawn( policy , FibChild2(policy,n-2) ) );
          policy.respawn( this );
        }
        else {
          // Spawn new children and respawn myself to sum their results:
          // result = Fib(n-1) + Fib(n-2)
          // result = ( Fib(n-2) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = ( ( Fib(n-3) + Fib(n-4) ) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = 3 * Fib(n-3) + 2 * Fib(n-4)
          has_nested = 4 ;
          // Kokkos::Experimental::respawn implements the following steps:
          policy.clear_dependence( this );
          policy.add_dependence( this , Kokkos::Experimental::spawn( policy , FibChild2(policy,n-3) ) );
          policy.add_dependence( this , Kokkos::Experimental::spawn( policy , FibChild2(policy,n-4) ) );
          policy.respawn( this );
        }
     }
     else if ( 2 == has_nested || 4 == has_nested ) {
        const Kokkos::Experimental::Future<long,ExecSpace> fib_a = policy.get_dependence(this,0);
        const Kokkos::Experimental::Future<long,ExecSpace> fib_b = policy.get_dependence(this,1);

        result = ( has_nested == 2 ) ? fib_a.get() + fib_b.get()
                                     : 3 * fib_a.get() + 2 * fib_b.get() ;

        has_nested = -1 ;
      }
      else {
        fprintf(stderr,"FibChild2(%ld) execution error\n",(long)n);
        fflush(stderr);
      }
    }
};

namespace {

long eval_fib( long n )
{
  if ( n < 2 ) return n ;

  std::vector<long> fib(n+1);

  fib[0] = 0 ;
  fib[1] = 1 ;

  for ( long i = 2 ; i <= n ; ++i ) { fib[i] = fib[i-2] + fib[i-1]; }

  return fib[n];
}

}

template< class ExecSpace >
void test_fib( long n )
{
  Kokkos::Experimental::TaskPolicy<ExecSpace> policy(2);

  Kokkos::Experimental::Future<long,ExecSpace> f = Kokkos::Experimental::spawn( policy , FibChild<ExecSpace>(policy,n) );

  Kokkos::Experimental::wait( policy );

  if ( f.get() != eval_fib(n) ) {
    std::cout << "Fib(" << n << ") = " << f.get();
    std::cout << " != " << eval_fib(n);
    std::cout << std::endl ;
  }
}

template< class ExecSpace >
void test_fib2( long n )
{
  Kokkos::Experimental::TaskPolicy<ExecSpace> policy(2); // default dependence capacity

  Kokkos::Experimental::Future<long,ExecSpace> f = Kokkos::Experimental::spawn( policy , FibChild2<ExecSpace>(policy,n) );

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

  inline
  void operator()( int i , double & val ) const { val += m_x[i] * m_x[i] ; }

  void apply( double & dst ) const { dst = std::sqrt( dst ); }
};

template< class ExecSpace >
void test_norm2( const int n )
{
  Kokkos::Experimental::TaskPolicy< ExecSpace > policy ;

  double * const x = new double[n];

  for ( int i = 0 ; i < n ; ++i ) x[i] = 1 ;

  Kokkos::RangePolicy<ExecSpace> r(0,n);

  Kokkos::Experimental::Future<double,ExecSpace> f = Kokkos::Experimental::spawn_reduce( policy , r , Norm2<ExecSpace>(x) );

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

  Kokkos::Experimental::TaskPolicy< Space > policy ;

  Kokkos::Experimental::Future<int,Space> f[ NTEST ];

  for ( int i = 0 ; i < NTEST ; ++i ) {
    // Create task in the "constructing" state with capacity for 'n+1' dependences
    f[i] = policy.create( TaskDep<Space>(policy,0) , n + 1 );

    if ( f[i].get_task_state() != Kokkos::Experimental::TASK_STATE_CONSTRUCTING ) {
      Kokkos::Impl::throw_runtime_exception("get_task_state() != Kokkos::Experimental::TASK_STATE_CONSTRUCTING");
    }

    // Only use 'n' dependences

    for ( int j = 0 ; j < n ; ++j ) {

      Kokkos::Experimental::Future<int,Space> nested = policy.create( TaskDep<Space>(policy,j+1) );

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
  typedef Kokkos::Experimental::Future<ExecSpace>      future_type ;
  typedef Kokkos::View<long*,ExecSpace>                view_type ;

  policy_type  policy ;
  future_type  future ;

  view_type  result ;
  const long nvalue ;

  TaskTeam( const policy_type & arg_policy
          , const view_type   & arg_result
          , const long          arg_nvalue )
    : policy(arg_policy)
    , future()
    , result( arg_result )
    , nvalue( arg_nvalue )
    {}

  inline
  void apply( const typename policy_type::member_type & member )
    {
      const long end   = nvalue + 1 ;
      const long begin = 0 < end - SPAN ? end - SPAN : 0 ;

      if ( 0 < begin && future.get_task_state() == Kokkos::Experimental::TASK_STATE_NULL ) {
        if ( member.team_rank() == 0 ) {
          future = policy.spawn( policy.create_team( TaskTeam( policy , result , begin - 1 ) ) );
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

  TaskTeamValue( const policy_type & arg_policy
               , const view_type   & arg_result
               , const long          arg_nvalue )
    : policy(arg_policy)
    , future()
    , result( arg_result )
    , nvalue( arg_nvalue )
    {}

  inline
  void apply( const typename policy_type::member_type & member , value_type & final )
    {
      const long end   = nvalue + 1 ;
      const long begin = 0 < end - SPAN ? end - SPAN : 0 ;

      if ( 0 < begin && future.get_task_state() == Kokkos::Experimental::TASK_STATE_NULL ) {
        if ( member.team_rank() == 0 ) {
          future = policy.spawn( policy.create_team( TaskTeamValue( policy , result , begin - 1 ) ) );
          policy.clear_dependence( this );
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

  policy_type  policy ;
  view_type    result("result",n+1);

  future_type f = policy.spawn( policy.create_team( task_type( policy , result , n ) ) );

  Kokkos::Experimental::wait( policy );

  for ( long i = 0 ; i <= n ; ++i ) {
    const long answer = i + 1 ;
    if ( result(i) != answer ) {
      std::cerr << "test_task_team void ERROR result(" << i << ") = " << result(i) << " != " << answer << std::endl ;
    }
  }

  future_value_type fv = policy.spawn( policy.create_team( task_value_type( policy , result , n ) ) );

  Kokkos::Experimental::wait( policy );

  if ( fv.get() != n + 1 ) {
    std::cerr << "test_task_team value ERROR future = " << fv.get() << " != " << n + 1 << std::endl ;
  }
  for ( long i = 0 ; i <= n ; ++i ) {
    const long answer = i + 1 ;
    if ( result(i) != answer ) {
      std::cerr << "test_task_team value ERROR result(" << i << ") = " << result(i) << " != " << answer << std::endl ;
    }
  }
}

//----------------------------------------------------------------------------

} // namespace TestTaskPolicy

#endif /* #ifndef KOKKOS_UNITTEST_TASKPOLICY_HPP */


