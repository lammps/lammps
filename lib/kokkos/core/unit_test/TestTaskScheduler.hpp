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

#ifndef KOKKOS_UNITTEST_TASKSCHEDULER_HPP
#define KOKKOS_UNITTEST_TASKSCHEDULER_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_TASKDAG )
#include <Kokkos_Core.hpp>
#include <cstdio>
#include <iostream>
#include <cmath>


namespace TestTaskScheduler {

namespace {

inline
long eval_fib( long n )
{
  constexpr long mask = 0x03;

  long fib[4] = { 0, 1, 1, 2 };

  for ( long i = 2; i <= n; ++i ) {
    fib[ i & mask ] = fib[ ( i - 1 ) & mask ] + fib[ ( i - 2 ) & mask ];
  }

  return fib[ n & mask ];
}

}

template< typename Space >
struct TestFib
{
  typedef Kokkos::TaskScheduler< Space >  sched_type;
  typedef Kokkos::Future< long, Space >   future_type;
  typedef long                            value_type;

  sched_type  sched;
  future_type fib_m1;
  future_type fib_m2;
  const value_type n;

  KOKKOS_INLINE_FUNCTION
  TestFib( const sched_type & arg_sched, const value_type arg_n )
    : sched( arg_sched ), fib_m1(), fib_m2(), n( arg_n ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( typename sched_type::member_type &, value_type & result )
  {
#if 0
    printf( "\nTestFib(%ld) %d %d\n", n, int( !fib_m1.is_null() ), int( !fib_m2.is_null() ) );
#endif

    if ( n < 2 ) {
      result = n;
    }
    else if ( !fib_m2.is_null() && !fib_m1.is_null() ) {
      result = fib_m1.get() + fib_m2.get();
    }
    else {
      // Spawn new children and respawn myself to sum their results.
      // Spawn lower value at higher priority as it has a shorter
      // path to completion.

      fib_m2 = Kokkos::task_spawn( Kokkos::TaskSingle( sched, Kokkos::TaskPriority::High )
                                 , TestFib( sched, n - 2 ) );

      fib_m1 = Kokkos::task_spawn( Kokkos::TaskSingle( sched )
                                 , TestFib( sched, n - 1 ) );

      Kokkos::Future< Space > dep[] = { fib_m1, fib_m2 };
      Kokkos::Future< Space > fib_all = Kokkos::when_all( dep, 2 );

      if ( !fib_m2.is_null() && !fib_m1.is_null() && !fib_all.is_null() ) {
        // High priority to retire this branch.
        Kokkos::respawn( this, fib_all, Kokkos::TaskPriority::High );
      }
      else {
#if 1
        printf( "TestFib(%ld) insufficient memory alloc_capacity(%d) task_max(%d) task_accum(%ld)\n"
               , n
               , sched.allocation_capacity()
               , sched.allocated_task_count_max()
               , sched.allocated_task_count_accum()
               );
#endif

        Kokkos::abort( "TestFib insufficient memory" );

      }
    }
  }

  static void run( int i, size_t MemoryCapacity = 16000 )
  {
    typedef typename sched_type::memory_space memory_space;

    enum { MinBlockSize   =   64 };
    enum { MaxBlockSize   = 1024 };
    enum { SuperBlockSize = 4096 };

    sched_type root_sched( memory_space()
                         , MemoryCapacity
                         , MinBlockSize
                         , std::min(size_t(MaxBlockSize),MemoryCapacity)
                         , std::min(size_t(SuperBlockSize),MemoryCapacity) );

    future_type f = Kokkos::host_spawn( Kokkos::TaskSingle( root_sched )
                                      , TestFib( root_sched, i ) );

    Kokkos::wait( root_sched );

    ASSERT_EQ( eval_fib( i ), f.get() );

#if 0
    fprintf( stdout, "\nTestFib::run(%d) spawn_size(%d) when_all_size(%d) alloc_capacity(%d) task_max(%d) task_accum(%ld)\n"
           , i
           , int(root_sched.template spawn_allocation_size<TestFib>())
           , int(root_sched.when_all_allocation_size(2))
           , root_sched.allocation_capacity()
           , root_sched.allocated_task_count_max()
           , root_sched.allocated_task_count_accum()
           );
    fflush( stdout );
#endif
  }
};

} // namespace TestTaskScheduler

//----------------------------------------------------------------------------

namespace TestTaskScheduler {

template< class Space >
struct TestTaskSpawn {
  typedef Kokkos::TaskScheduler< Space >  sched_type;
  typedef Kokkos::Future< Space >         future_type;
  typedef void                            value_type;

  sched_type   m_sched ;
  future_type  m_future ;

  KOKKOS_INLINE_FUNCTION
  TestTaskSpawn( const sched_type & arg_sched
               , const future_type & arg_future
               )
    : m_sched( arg_sched )
    , m_future( arg_future )
    {}

  KOKKOS_INLINE_FUNCTION
  void operator()( typename sched_type::member_type & )
  {
    if ( ! m_future.is_null() ) {
      Kokkos::task_spawn( Kokkos::TaskSingle( m_sched ) , TestTaskSpawn( m_sched , future_type() ) );
    }
  }

  static void run()
  {
    typedef typename sched_type::memory_space memory_space;

    enum { MemoryCapacity = 16000 };
    enum { MinBlockSize   =   64 };
    enum { MaxBlockSize   = 1024 };
    enum { SuperBlockSize = 4096 };

    sched_type sched( memory_space()
                    , MemoryCapacity
                    , MinBlockSize
                    , MaxBlockSize
                    , SuperBlockSize );

    auto f = Kokkos::host_spawn( Kokkos::TaskSingle( sched ), TestTaskSpawn( sched, future_type() ) );
    Kokkos::host_spawn( Kokkos::TaskSingle( f ), TestTaskSpawn( sched, f ) );

    Kokkos::wait( sched );
  }
};

template< class Space >
struct TestTaskDependence {
  typedef Kokkos::TaskScheduler< Space >  sched_type;
  typedef Kokkos::Future< Space >         future_type;
  typedef Kokkos::View< long, Space >     accum_type;
  typedef void                            value_type;

  sched_type  m_sched;
  accum_type  m_accum;
  long        m_count;

  KOKKOS_INLINE_FUNCTION
  TestTaskDependence( long n
                    , const sched_type & arg_sched
                    , const accum_type & arg_accum )
    : m_sched( arg_sched )
    , m_accum( arg_accum )
    , m_count( n ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( typename sched_type::member_type & )
  {
    enum { CHUNK = 8 };
    const int n = CHUNK < m_count ? CHUNK : m_count;

    if ( 1 < m_count ) {

      const int increment = ( m_count + n - 1 ) / n;

      future_type f =
        m_sched.when_all( n , [this,increment]( int i ) {
          const long inc   = increment ;
          const long begin = i * inc ;
          const long count = begin + inc < m_count ? inc : m_count - begin ;

          return Kokkos::task_spawn
            ( Kokkos::TaskSingle( m_sched )
            , TestTaskDependence( count, m_sched, m_accum ) );
        });

      m_count = 0;

      Kokkos::respawn( this, f );
    }
    else if ( 1 == m_count ) {
      Kokkos::atomic_increment( & m_accum() );
    }
  }

  static void run( int n )
  {
    typedef typename sched_type::memory_space memory_space;

    enum { MemoryCapacity = 16000 };
    enum { MinBlockSize   =   64 };
    enum { MaxBlockSize   = 1024 };
    enum { SuperBlockSize = 4096 };

    sched_type sched( memory_space()
                    , MemoryCapacity
                    , MinBlockSize
                    , MaxBlockSize
                    , SuperBlockSize );

    accum_type accum( "accum" );

    typename accum_type::HostMirror host_accum = Kokkos::create_mirror_view( accum );

    Kokkos::host_spawn( Kokkos::TaskSingle( sched ), TestTaskDependence( n, sched, accum ) );

    Kokkos::wait( sched );

    Kokkos::deep_copy( host_accum, accum );

    ASSERT_EQ( host_accum(), n );
  }
};

} // namespace TestTaskScheduler

//----------------------------------------------------------------------------

namespace TestTaskScheduler {

template< class ExecSpace >
struct TestTaskTeam {
  //enum { SPAN = 8 };
  enum { SPAN = 33 };
  //enum { SPAN = 1 };

  typedef void                                value_type;
  typedef Kokkos::TaskScheduler< ExecSpace >  sched_type;
  typedef Kokkos::Future< ExecSpace >         future_type;
  typedef Kokkos::View< long*, ExecSpace >    view_type;

  sched_type   sched;
  future_type  future;

  view_type   parfor_result;
  view_type   parreduce_check;
  view_type   parscan_result;
  view_type   parscan_check;
  const long  nvalue;

  KOKKOS_INLINE_FUNCTION
  TestTaskTeam( const sched_type & arg_sched
              , const view_type  & arg_parfor_result
              , const view_type  & arg_parreduce_check
              , const view_type  & arg_parscan_result
              , const view_type  & arg_parscan_check
              , const long         arg_nvalue )
    : sched( arg_sched )
    , future()
    , parfor_result( arg_parfor_result )
    , parreduce_check( arg_parreduce_check )
    , parscan_result( arg_parscan_result )
    , parscan_check( arg_parscan_check )
    , nvalue( arg_nvalue ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( typename sched_type::member_type & member )
  {
    const long end   = nvalue + 1;
    const long begin = 0 < end - SPAN ? end - SPAN : 0;

    if ( 0 < begin && future.is_null() ) {
      if ( member.team_rank() == 0 ) {
        future = Kokkos::task_spawn( Kokkos::TaskTeam( sched )
                                   , TestTaskTeam( sched
                                                 , parfor_result
                                                 , parreduce_check
                                                 , parscan_result
                                                 , parscan_check
                                                 , begin - 1 )
                                   );

        #ifndef __HCC_ACCELERATOR__
        assert( !future.is_null() );
        #endif

        Kokkos::respawn( this, future );
      }

      return;
    }

    Kokkos::parallel_for( Kokkos::TeamThreadRange( member, begin, end )
                        , [&] ( int i ) { parfor_result[i] = i; }
                        );

    // Test parallel_reduce without join.

    long tot = 0;
    long expected = ( begin + end - 1 ) * ( end - begin ) * 0.5;

    Kokkos::parallel_reduce( Kokkos::TeamThreadRange( member, begin, end )
                           , [&] ( int i, long & res ) { res += parfor_result[i]; }
                           , tot
                           );

    Kokkos::parallel_for( Kokkos::TeamThreadRange( member, begin, end )
                        , [&] ( int i ) { parreduce_check[i] = expected - tot; }
                        );

    // Test parallel_reduce with join.

    tot = 0;
    Kokkos::parallel_reduce( Kokkos::TeamThreadRange( member, begin, end )
                           , [&] ( int i, long & res ) { res += parfor_result[i]; }
                           , Kokkos::Experimental::Sum<long>( tot )
                           );

    Kokkos::parallel_for( Kokkos::TeamThreadRange( member, begin, end )
                        , [&] ( int i ) { parreduce_check[i] += expected - tot; }
                        );

    // Test parallel_scan.

    // Exclusive scan.
    Kokkos::parallel_scan<long>( Kokkos::TeamThreadRange( member, begin, end )
                               , [&] ( int i, long & val, const bool final )
    {
      if ( final ) { parscan_result[i] = val; }

      val += i;
    });

    // Wait for 'parscan_result' before testing it.
    member.team_barrier();

    if ( member.team_rank() == 0 ) {
      for ( long i = begin; i < end; ++i ) {
        parscan_check[i] = ( i * ( i - 1 ) - begin * ( begin - 1 ) ) * 0.5 - parscan_result[i];
      }
    }

    // Don't overwrite 'parscan_result' until it has been tested.
    member.team_barrier();

    // Inclusive scan.
    Kokkos::parallel_scan<long>( Kokkos::TeamThreadRange( member, begin, end )
                               , [&] ( int i, long & val, const bool final )
    {
      val += i;

      if ( final ) { parscan_result[i] = val; }
    });

    // Wait for 'parscan_result' before testing it.
    member.team_barrier();

    if ( member.team_rank() == 0 ) {
      for ( long i = begin; i < end; ++i ) {
        parscan_check[i] += ( i * ( i + 1 ) - begin * ( begin - 1 ) ) * 0.5 - parscan_result[i];
      }
    }

    // ThreadVectorRange check.
/*
    long result = 0;
    expected = ( begin + end - 1 ) * ( end - begin ) * 0.5;
    Kokkos::parallel_reduce( Kokkos::TeamThreadRange( member, 0, 1 )
                           , [&] ( const int i, long & outerUpdate )
    {
      long sum_j = 0.0;

      Kokkos::parallel_reduce( Kokkos::ThreadVectorRange( member, end - begin )
                             , [&] ( const int j, long & innerUpdate )
      {
        innerUpdate += begin + j;
      }, sum_j );

      outerUpdate += sum_j;
    }, result );

    Kokkos::parallel_for( Kokkos::TeamThreadRange( member, begin, end )
                        , [&] ( int i )
    {
      parreduce_check[i] += result - expected;
    });
*/

  }

  static void run( long n )
  {
    const unsigned memory_capacity = 400000;

    enum { MinBlockSize   =   64 };
    enum { MaxBlockSize   = 1024 };
    enum { SuperBlockSize = 4096 };

    sched_type root_sched( typename sched_type::memory_space()
                         , memory_capacity
                         , MinBlockSize
                         , MaxBlockSize
                         , SuperBlockSize );

    view_type root_parfor_result( "parfor_result", n + 1 );
    view_type root_parreduce_check( "parreduce_check", n + 1 );
    view_type root_parscan_result( "parscan_result", n + 1 );
    view_type root_parscan_check( "parscan_check", n + 1 );

    typename view_type::HostMirror
      host_parfor_result = Kokkos::create_mirror_view( root_parfor_result );
    typename view_type::HostMirror
      host_parreduce_check = Kokkos::create_mirror_view( root_parreduce_check );
    typename view_type::HostMirror
      host_parscan_result = Kokkos::create_mirror_view( root_parscan_result );
    typename view_type::HostMirror
      host_parscan_check = Kokkos::create_mirror_view( root_parscan_check );

    future_type f = Kokkos::host_spawn( Kokkos::TaskTeam( root_sched )
                                      , TestTaskTeam( root_sched
                                                    , root_parfor_result
                                                    , root_parreduce_check
                                                    , root_parscan_result
                                                    , root_parscan_check
                                                    , n )
                                      );

    Kokkos::wait( root_sched );

    Kokkos::deep_copy( host_parfor_result, root_parfor_result );
    Kokkos::deep_copy( host_parreduce_check, root_parreduce_check );
    Kokkos::deep_copy( host_parscan_result, root_parscan_result );
    Kokkos::deep_copy( host_parscan_check, root_parscan_check );

    long error_count = 0 ;

    for ( long i = 0; i <= n; ++i ) {
      const long answer = i;

      if ( host_parfor_result( i ) != answer ) {
        ++error_count ;
        std::cerr << "TestTaskTeam::run ERROR parallel_for result(" << i << ") = "
                  << host_parfor_result( i ) << " != " << answer << std::endl;
      }

      if ( host_parreduce_check( i ) != 0 ) {
        ++error_count ;
        std::cerr << "TestTaskTeam::run ERROR parallel_reduce check(" << i << ") = "
                  << host_parreduce_check( i ) << " != 0" << std::endl;
      }

      if ( host_parscan_check( i ) != 0 ) {
        ++error_count ;
        std::cerr << "TestTaskTeam::run ERROR parallel_scan check(" << i << ") = "
                  << host_parscan_check( i ) << " != 0" << std::endl;
      }
    }

    ASSERT_EQ( 0L , error_count );
  }
};

template< class ExecSpace >
struct TestTaskTeamValue {
  enum { SPAN = 8 };

  typedef long                                     value_type;
  typedef Kokkos::TaskScheduler< ExecSpace >       sched_type;
  typedef Kokkos::Future< value_type, ExecSpace >  future_type;
  typedef Kokkos::View< long*, ExecSpace >         view_type;

  sched_type   sched;
  future_type  future;

  view_type   result;
  const long  nvalue;

  KOKKOS_INLINE_FUNCTION
  TestTaskTeamValue( const sched_type & arg_sched
                   , const view_type  & arg_result
                   , const long         arg_nvalue )
    : sched( arg_sched )
    , future()
    , result( arg_result )
    , nvalue( arg_nvalue ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( typename sched_type::member_type const & member
                 , value_type & final )
  {
    const long end   = nvalue + 1;
    const long begin = 0 < end - SPAN ? end - SPAN : 0;

    if ( 0 < begin && future.is_null() ) {
      if ( member.team_rank() == 0 ) {
        future = sched.task_spawn( TestTaskTeamValue( sched, result, begin - 1 )
                                 , Kokkos::TaskTeam );

        assert( !future.is_null() );

        sched.respawn( this , future );
      }

      return;
    }

    Kokkos::parallel_for( Kokkos::TeamThreadRange( member, begin, end )
                        , [&] ( int i ) { result[i] = i + 1; }
                        );

    if ( member.team_rank() == 0 ) {
      final = result[nvalue];
    }

    Kokkos::memory_fence();
  }

  static void run( long n )
  {
    const unsigned memory_capacity = 100000;

    enum { MinBlockSize   =   64 };
    enum { MaxBlockSize   = 1024 };
    enum { SuperBlockSize = 4096 };

    sched_type root_sched( typename sched_type::memory_space()
                         , memory_capacity
                         , MinBlockSize
                         , MaxBlockSize
                         , SuperBlockSize );

    view_type root_result( "result", n + 1 );

    typename view_type::HostMirror host_result = Kokkos::create_mirror_view( root_result );

    future_type fv = root_sched.host_spawn( TestTaskTeamValue( root_sched, root_result, n )
                                          , Kokkos::TaskTeam );

    Kokkos::wait( root_sched );

    Kokkos::deep_copy( host_result, root_result );

    if ( fv.get() != n + 1 ) {
      std::cerr << "TestTaskTeamValue ERROR future = "
                << fv.get() << " != " << n + 1 << std::endl;
    }

    for ( long i = 0; i <= n; ++i ) {
      const long answer = i + 1;

      if ( host_result( i ) != answer ) {
        std::cerr << "TestTaskTeamValue ERROR result(" << i << ") = "
                  << host_result( i ) << " != " << answer << std::endl;
      }
    }
  }
};

} // namespace TestTaskScheduler

namespace Test {

TEST_F( TEST_CATEGORY, task_fib )
{
  const int N = 27 ;
  for ( int i = 0; i < N; ++i ) {
    TestTaskScheduler::TestFib< TEST_EXECSPACE >::run( i , ( i + 1 ) * ( i + 1 ) * 2000 );
  }
}

TEST_F( TEST_CATEGORY, task_depend )
{
  for ( int i = 0; i < 25; ++i ) {
    TestTaskScheduler::TestTaskDependence< TEST_EXECSPACE >::run( i );
  }
}

TEST_F( TEST_CATEGORY, task_team )
{
  TestTaskScheduler::TestTaskTeam< TEST_EXECSPACE >::run( 1000 );
  //TestTaskScheduler::TestTaskTeamValue< TEST_EXECSPACE >::run( 1000 ); // Put back after testing.
}

}

#endif // #if defined( KOKKOS_ENABLE_TASKDAG )
#endif // #ifndef KOKKOS_UNITTEST_TASKSCHEDULER_HPP

