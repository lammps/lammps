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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <cstdio>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {

namespace {

template< class ExecSpace, class ScheduleType >
struct TestTeamPolicy {
  typedef typename Kokkos::TeamPolicy< ScheduleType,  ExecSpace >::member_type team_member;
  typedef Kokkos::View< int**, ExecSpace > view_type;

  view_type m_flags;

  TestTeamPolicy( const size_t league_size )
    : m_flags( Kokkos::ViewAllocateWithoutInitializing( "flags" ),
               Kokkos::TeamPolicy< ScheduleType,  ExecSpace >(1,1).team_size_max( *this, Kokkos::ParallelReduceTag() ),
               league_size ) {}

  struct VerifyInitTag {};

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & member ) const
  {
    const int tid = member.team_rank() + member.team_size() * member.league_rank();

    m_flags( member.team_rank(), member.league_rank() ) = tid;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const VerifyInitTag &, const team_member & member ) const
  {
    const int tid = member.team_rank() + member.team_size() * member.league_rank();

    if ( tid != m_flags( member.team_rank(), member.league_rank() ) ) {
      printf( "TestTeamPolicy member(%d,%d) error %d != %d\n",
               member.league_rank(), member.team_rank(),
               tid, m_flags( member.team_rank(), member.league_rank() ) );
    }
  }

  // Included for test_small_league_size.
  TestTeamPolicy() : m_flags() {}

  // Included for test_small_league_size.
  struct NoOpTag {};

  KOKKOS_INLINE_FUNCTION
  void operator()( const NoOpTag &, const team_member & member ) const {}


  static void test_small_league_size() {
    int bs = 8; // batch size (number of elements per batch)
    int ns = 16; // total number of "problems" to process

    // Calculate total scratch memory space size.
    const int level = 0;
    int mem_size = 960;
    const int num_teams = ns / bs;
    Kokkos::TeamPolicy< ExecSpace, NoOpTag > policy( num_teams, Kokkos::AUTO() );

    Kokkos::parallel_for( policy.set_scratch_size( level, Kokkos::PerTeam( mem_size ), Kokkos::PerThread( 0 ) ),
                          TestTeamPolicy() );
  }

  static void test_for( const size_t league_size )
  {
    TestTeamPolicy functor( league_size );
    typedef Kokkos::TeamPolicy< ScheduleType,  ExecSpace > policy_type;
    typedef Kokkos::TeamPolicy< ScheduleType,  ExecSpace, VerifyInitTag > policy_type_init;

    const int team_size = policy_type(league_size,1).team_size_max( functor, Kokkos::ParallelForTag() );
    const int team_size_init = policy_type_init(league_size,1).team_size_max( functor, Kokkos::ParallelForTag() );

    Kokkos::parallel_for( policy_type( league_size, team_size ), functor );
    Kokkos::parallel_for( policy_type_init( league_size, team_size_init ), functor );

    test_small_league_size();
  }

  struct ReduceTag {};

  typedef long value_type;

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & member, value_type & update ) const
  {
    update += member.team_rank() + member.team_size() * member.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const ReduceTag &, const team_member & member, value_type & update ) const
  {
    update += 1 + member.team_rank() + member.team_size() * member.league_rank();
  }

  static void test_reduce( const size_t league_size )
  {
    TestTeamPolicy functor( league_size );

    typedef Kokkos::TeamPolicy< ScheduleType,  ExecSpace > policy_type;
    typedef Kokkos::TeamPolicy< ScheduleType,  ExecSpace, ReduceTag > policy_type_reduce;

    const int team_size = policy_type_reduce(league_size,1).team_size_max( functor, Kokkos::ParallelReduceTag() );

    const long N = team_size * league_size;

    long total = 0;

    Kokkos::parallel_reduce( policy_type( league_size, team_size ), functor, total );
    ASSERT_EQ( size_t( ( N - 1 ) * ( N ) ) / 2, size_t( total ) );

    Kokkos::parallel_reduce( policy_type_reduce( league_size, team_size ), functor, total );
    ASSERT_EQ( ( size_t( N ) * size_t( N + 1 ) ) / 2, size_t( total ) );
  }
};

} // namespace

} // namespace Test

/*--------------------------------------------------------------------------*/

namespace Test {

template< typename ScalarType, class DeviceType, class ScheduleType >
class ReduceTeamFunctor
{
public:
  typedef DeviceType                                           execution_space;
  typedef Kokkos::TeamPolicy< ScheduleType, execution_space >  policy_type;
  typedef typename execution_space::size_type                  size_type;

  struct value_type {
    ScalarType value[3];
  };

  const size_type nwork;

  KOKKOS_INLINE_FUNCTION
  ReduceTeamFunctor( const size_type & arg_nwork ) : nwork( arg_nwork ) {}

  KOKKOS_INLINE_FUNCTION
  ReduceTeamFunctor( const ReduceTeamFunctor & rhs ) : nwork( rhs.nwork ) {}

  KOKKOS_INLINE_FUNCTION
  void init( value_type & dst ) const
  {
    dst.value[0] = 0;
    dst.value[1] = 0;
    dst.value[2] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst, const volatile value_type & src ) const
  {
    dst.value[0] += src.value[0];
    dst.value[1] += src.value[1];
    dst.value[2] += src.value[2];
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const typename policy_type::member_type ind, value_type & dst ) const
  {
    const int thread_rank = ind.team_rank() + ind.team_size() * ind.league_rank();
    const int thread_size = ind.team_size() * ind.league_size();
    const int chunk = ( nwork + thread_size - 1 ) / thread_size;

    size_type iwork = chunk * thread_rank;
    const size_type iwork_end = iwork + chunk < nwork ? iwork + chunk : nwork;

    for ( ; iwork < iwork_end; ++iwork ) {
      dst.value[0] += 1;
      dst.value[1] += iwork + 1;
      dst.value[2] += nwork - iwork;
    }
  }
};

} // namespace Test

namespace {

template< typename ScalarType, class DeviceType, class ScheduleType >
class TestReduceTeam
{
public:
  typedef DeviceType                                            execution_space;
  typedef Kokkos::TeamPolicy< ScheduleType,  execution_space >  policy_type;
  typedef typename execution_space::size_type                   size_type;

  TestReduceTeam( const size_type & nwork ) { run_test( nwork ); }

  void run_test( const size_type & nwork )
  {
    typedef Test::ReduceTeamFunctor< ScalarType, execution_space, ScheduleType> functor_type;
    typedef typename functor_type::value_type value_type;
    typedef Kokkos::View< value_type, Kokkos::HostSpace, Kokkos::MemoryUnmanaged > result_type;

    enum { Count = 3 };
    enum { Repeat = 100 };

    value_type result[ Repeat ];

    const unsigned long nw   = nwork;
    const unsigned long nsum = nw % 2 ? nw * ( ( nw + 1 ) / 2 )
                                      : ( nw / 2 ) * ( nw + 1 );

    policy_type team_exec( nw, 1 );

    const unsigned team_size   = team_exec.team_size_recommended( functor_type( nwork ), Kokkos::ParallelReduceTag() );
    const unsigned league_size = ( nwork + team_size - 1 ) / team_size;

    team_exec = policy_type( league_size, team_size );

    for ( unsigned i = 0; i < Repeat; ++i ) {
      result_type tmp( & result[i] );
      Kokkos::parallel_reduce( team_exec, functor_type( nwork ), tmp );
    }

    execution_space::fence();

    for ( unsigned i = 0; i < Repeat; ++i ) {
      for ( unsigned j = 0; j < Count; ++j ) {
        const unsigned long correct = 0 == j % 3 ? nw : nsum;
        ASSERT_EQ( (ScalarType) correct, result[i].value[j] );
      }
    }
  }
};

} // namespace

/*--------------------------------------------------------------------------*/

namespace Test {

template< class DeviceType, class ScheduleType >
class ScanTeamFunctor
{
public:
  typedef DeviceType                                            execution_space;
  typedef Kokkos::TeamPolicy< ScheduleType,  execution_space >  policy_type;
  typedef long int                                              value_type;

  Kokkos::View< value_type, execution_space > accum;
  Kokkos::View< value_type, execution_space > total;

  ScanTeamFunctor() : accum( "accum" ), total( "total" ) {}

  KOKKOS_INLINE_FUNCTION
  void init( value_type & error ) const { error = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( value_type volatile & error, value_type volatile const & input ) const
  { if ( input ) error = 1; }

  struct JoinMax {
    typedef long int value_type;

    KOKKOS_INLINE_FUNCTION
    void join( value_type volatile & dst, value_type volatile const & input ) const
    { if ( dst < input ) dst = input; }
  };

  KOKKOS_INLINE_FUNCTION
  void operator()( const typename policy_type::member_type ind, value_type & error ) const
  {
    if ( 0 == ind.league_rank() && 0 == ind.team_rank() ) {
      const long int thread_count = ind.league_size() * ind.team_size();
      total() = ( thread_count * ( thread_count + 1 ) ) / 2;
    }

    // Team max:
    int long m = (long int) ( ind.league_rank() + ind.team_rank() );
    ind.team_reduce(  Kokkos::Max<int long>(m) );

    if ( m != ind.league_rank() + ( ind.team_size() - 1 ) ) {
      printf( "ScanTeamFunctor[%d.%d of %d.%d] reduce_max_answer(%ld) != reduce_max(%ld)\n",
               ind.league_rank(), ind.team_rank(),
               ind.league_size(), ind.team_size(),
               (long int) ( ind.league_rank() + ( ind.team_size() - 1 ) ), m );
    }

    // Scan:
    const long int answer =
      ( ind.league_rank() + 1 ) * ind.team_rank() + ( ind.team_rank() * ( ind.team_rank() + 1 ) ) / 2;

    const long int result =
      ind.team_scan( ind.league_rank() + 1 + ind.team_rank() + 1 );

    const long int result2 =
      ind.team_scan( ind.league_rank() + 1 + ind.team_rank() + 1 );

    if ( answer != result || answer != result2 ) {
      printf( "ScanTeamFunctor[%d.%d of %d.%d] answer(%ld) != scan_first(%ld) or scan_second(%ld)\n",
              ind.league_rank(), ind.team_rank(),
              ind.league_size(), ind.team_size(),
              answer, result, result2 );

      error = 1;
    }

    const long int thread_rank = ind.team_rank() +
                                 ind.team_size() * ind.league_rank();
    ind.team_scan( 1 + thread_rank, accum.data() );
  }
};

template< class DeviceType, class ScheduleType >
class TestScanTeam
{
public:
  typedef DeviceType                                            execution_space;
  typedef long int                                              value_type;
  typedef Kokkos::TeamPolicy< ScheduleType,  execution_space >  policy_type;
  typedef Test::ScanTeamFunctor<DeviceType, ScheduleType>       functor_type;

  TestScanTeam( const size_t nteam ) { run_test( nteam ); }

  void run_test( const size_t nteam )
  {
    typedef Kokkos::View< long int, Kokkos::HostSpace, Kokkos::MemoryUnmanaged >  result_type;

    const unsigned REPEAT = 100000;
    unsigned Repeat;

    if ( nteam == 0 ) {
      Repeat = 1;
    }
    else {
      Repeat = ( REPEAT + nteam - 1 ) / nteam; // Error here.
    }

    functor_type functor;

    policy_type team_exec( nteam, 1);
    team_exec = policy_type(nteam, team_exec.team_size_max(functor, Kokkos::ParallelReduceTag()));

    for ( unsigned i = 0; i < Repeat; ++i ) {
      long int accum = 0;
      long int total = 0;
      long int error = 0;
      Kokkos::deep_copy( functor.accum, total );

      Kokkos::parallel_reduce( team_exec, functor, result_type( & error ) );
      DeviceType::fence();

      Kokkos::deep_copy( accum, functor.accum );
      Kokkos::deep_copy( total, functor.total );

      ASSERT_EQ( error, 0 );
      ASSERT_EQ( total, accum );
    }

    execution_space::fence();
  }
};

} // namespace Test

/*--------------------------------------------------------------------------*/

namespace Test {

template< class ExecSpace, class ScheduleType >
struct SharedTeamFunctor {

  typedef ExecSpace                                             execution_space;
  typedef int                                                   value_type;
  typedef Kokkos::TeamPolicy< ScheduleType,  execution_space >  policy_type;

  enum { SHARED_COUNT = 1000 };

  typedef typename ExecSpace::scratch_memory_space  shmem_space;

  // TBD: MemoryUnmanaged should be the default for shared memory space.
  typedef Kokkos::View< int*, shmem_space, Kokkos::MemoryUnmanaged > shared_int_array_type;

  // Tell how much shared memory will be required by this functor.
  inline
  unsigned team_shmem_size( int team_size ) const
  {
    return shared_int_array_type::shmem_size( SHARED_COUNT ) +
           shared_int_array_type::shmem_size( SHARED_COUNT );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const typename policy_type::member_type & ind, value_type & update ) const
  {
    const shared_int_array_type shared_A( ind.team_shmem(), SHARED_COUNT );
    const shared_int_array_type shared_B( ind.team_shmem(), SHARED_COUNT );

    if ( ( shared_A.data() == nullptr && SHARED_COUNT > 0 ) ||
         ( shared_B.data() == nullptr && SHARED_COUNT > 0 ) )
    {
      printf ("member( %d/%d , %d/%d ) Failed to allocate shared memory of size %lu\n"
             , ind.league_rank()
             , ind.league_size()
             , ind.team_rank()
             , ind.team_size()
             , static_cast<unsigned long>( SHARED_COUNT )
             );

      ++update; // Failure to allocate is an error.
    }
    else {
      for ( int i = ind.team_rank(); i < SHARED_COUNT; i += ind.team_size() ) {
        shared_A[i] = i + ind.league_rank();
        shared_B[i] = 2 * i + ind.league_rank();
      }

      ind.team_barrier();

      if ( ind.team_rank() + 1 == ind.team_size() ) {
        for ( int i = 0; i < SHARED_COUNT; ++i ) {
          if ( shared_A[i] != i + ind.league_rank() ) {
            ++update;
          }

          if ( shared_B[i] != 2 * i + ind.league_rank() ) {
            ++update;
          }
        }
      }
    }
  }
};

} // namespace Test

namespace {

template< class ExecSpace, class ScheduleType >
struct TestSharedTeam {
  TestSharedTeam() { run(); }

  void run()
  {
    typedef Test::SharedTeamFunctor<ExecSpace, ScheduleType> Functor;
    typedef Kokkos::View< typename Functor::value_type, Kokkos::HostSpace, Kokkos::MemoryUnmanaged > result_type;

    const size_t team_size = Kokkos::TeamPolicy< ScheduleType, ExecSpace >(8192, 1).team_size_max( Functor(),
        Kokkos::ParallelReduceTag() );

    Kokkos::TeamPolicy< ScheduleType, ExecSpace > team_exec( 8192 / team_size, team_size );

    typename Functor::value_type error_count = 0;

    Kokkos::parallel_reduce( team_exec, Functor(), result_type( & error_count ) );

    ASSERT_EQ( error_count, 0 );
  }
};

} // namespace

namespace Test {

#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
#if !defined(KOKKOS_ENABLE_CUDA) || ( 8000 <= CUDA_VERSION )
template< class MemorySpace, class ExecSpace, class ScheduleType >
struct TestLambdaSharedTeam {
  TestLambdaSharedTeam() { run(); }

  void run()
  {
    typedef Test::SharedTeamFunctor< ExecSpace, ScheduleType > Functor;
    //typedef Kokkos::View< typename Functor::value_type, Kokkos::HostSpace, Kokkos::MemoryUnmanaged > result_type;
    typedef Kokkos::View< typename Functor::value_type, MemorySpace, Kokkos::MemoryUnmanaged > result_type;

    typedef typename ExecSpace::scratch_memory_space shmem_space;

    // TBD: MemoryUnmanaged should be the default for shared memory space.
    typedef Kokkos::View< int*, shmem_space, Kokkos::MemoryUnmanaged > shared_int_array_type;

    const int SHARED_COUNT = 1000;
    int team_size = 1;

#ifdef KOKKOS_ENABLE_CUDA
    if ( std::is_same< ExecSpace, Kokkos::Cuda >::value ) team_size = 128;
#endif

    Kokkos::TeamPolicy< ScheduleType,  ExecSpace > team_exec( 8192 / team_size, team_size );
    team_exec = team_exec.set_scratch_size( 0, Kokkos::PerTeam( SHARED_COUNT * 2 * sizeof( int ) ) );

    typename Functor::value_type error_count = 0;

    Kokkos::parallel_reduce( team_exec, KOKKOS_LAMBDA
        ( const typename Kokkos::TeamPolicy< ScheduleType,  ExecSpace >::member_type & ind, int & update )
    {
      const shared_int_array_type shared_A( ind.team_shmem(), SHARED_COUNT );
      const shared_int_array_type shared_B( ind.team_shmem(), SHARED_COUNT );

      if ( ( shared_A.data () == nullptr && SHARED_COUNT > 0 ) ||
           ( shared_B.data () == nullptr && SHARED_COUNT > 0 ) )
      {
        printf( "Failed to allocate shared memory of size %lu\n",
                static_cast<unsigned long>( SHARED_COUNT ) );

        ++update; // Failure to allocate is an error.
      }
      else {
        for ( int i = ind.team_rank(); i < SHARED_COUNT; i += ind.team_size() ) {
          shared_A[i] = i + ind.league_rank();
          shared_B[i] = 2 * i + ind.league_rank();
        }

        ind.team_barrier();

        if ( ind.team_rank() + 1 == ind.team_size() ) {
          for ( int i = 0; i < SHARED_COUNT; ++i ) {
            if ( shared_A[i] != i + ind.league_rank() ) {
              ++update;
            }

            if ( shared_B[i] != 2 * i + ind.league_rank() ) {
              ++update;
            }
          }
        }
      }
    }, result_type( & error_count ) );

    ASSERT_EQ( error_count, 0 );
  }
};
#endif
#endif

} // namespace Test

namespace Test {

template< class ExecSpace, class ScheduleType >
struct ScratchTeamFunctor {

  typedef ExecSpace                                            execution_space;
  typedef int                                                  value_type;
  typedef Kokkos::TeamPolicy< ScheduleType, execution_space >  policy_type;

  enum { SHARED_TEAM_COUNT = 100 };
  enum { SHARED_THREAD_COUNT = 10 };

  typedef typename ExecSpace::scratch_memory_space shmem_space;

  // TBD: MemoryUnmanaged should be the default for shared memory space.
  typedef Kokkos::View< size_t*, shmem_space, Kokkos::MemoryUnmanaged > shared_int_array_type;

  KOKKOS_INLINE_FUNCTION
  void operator()( const typename policy_type::member_type & ind, value_type & update ) const
  {
    const shared_int_array_type scratch_ptr( ind.team_scratch( 1 ), 3 * ind.team_size() );
    const shared_int_array_type scratch_A( ind.team_scratch( 1 ), SHARED_TEAM_COUNT );
    const shared_int_array_type scratch_B( ind.thread_scratch( 1 ), SHARED_THREAD_COUNT );

    if ( ( scratch_ptr.data() == nullptr ) ||
         ( scratch_A.  data() == nullptr && SHARED_TEAM_COUNT > 0 ) ||
         ( scratch_B.  data() == nullptr && SHARED_THREAD_COUNT > 0 ) )
    {
      printf( "Failed to allocate shared memory of size %lu\n",
              static_cast<unsigned long>( SHARED_TEAM_COUNT ) );

      ++update; // Failure to allocate is an error.
    }
    else {
      Kokkos::parallel_for( Kokkos::TeamThreadRange( ind, 0, (int) SHARED_TEAM_COUNT ), [&] ( const int & i ) {
        scratch_A[i] = i + ind.league_rank();
      });

      for ( int i = 0; i < SHARED_THREAD_COUNT; i++ ) {
        scratch_B[i] = 10000 * ind.league_rank() + 100 * ind.team_rank() + i;
      }

      scratch_ptr[ind.team_rank()] = (size_t) scratch_A.data();
      scratch_ptr[ind.team_rank() + ind.team_size()] = (size_t) scratch_B.data();

      ind.team_barrier();

      for ( int i = 0; i < SHARED_TEAM_COUNT; i++ ) {
        if ( scratch_A[i] != size_t( i + ind.league_rank() ) ) ++update;
      }

      for ( int i = 0; i < ind.team_size(); i++ ) {
        if ( scratch_ptr[0] != scratch_ptr[i] ) ++update;
      }

      if ( scratch_ptr[1 + ind.team_size()] - scratch_ptr[0 + ind.team_size()] < SHARED_THREAD_COUNT * sizeof( size_t ) ) {
        ++update;
      }

      for ( int i = 1; i < ind.team_size(); i++ ) {
        if ( ( scratch_ptr[i + ind.team_size()] - scratch_ptr[i - 1 + ind.team_size()] ) !=
             ( scratch_ptr[1 + ind.team_size()] - scratch_ptr[0 + ind.team_size()] ) )
        {
          ++update;
        }
      }
    }
  }
};

} // namespace Test

namespace {

template< class ExecSpace, class ScheduleType >
struct TestScratchTeam {
  TestScratchTeam() { run(); }

  void run()
  {
    typedef Test::ScratchTeamFunctor<ExecSpace, ScheduleType> Functor;
    typedef Kokkos::View< typename Functor::value_type, Kokkos::HostSpace, Kokkos::MemoryUnmanaged >  result_type;
    typedef Kokkos::TeamPolicy< ScheduleType,  ExecSpace > p_type;


    typename Functor::value_type error_count = 0;

    int thread_scratch_size = Functor::shared_int_array_type::shmem_size( Functor::SHARED_THREAD_COUNT );

    p_type team_exec = p_type(8192,1).set_scratch_size( 1, Kokkos::PerTeam( Functor::shared_int_array_type::shmem_size( Functor::SHARED_TEAM_COUNT ) ),
                                                           Kokkos::PerThread( thread_scratch_size + 3*sizeof(int)));

    const size_t team_size = team_exec.team_size_max( Functor(), Kokkos::ParallelReduceTag() );

    int team_scratch_size   = Functor::shared_int_array_type::shmem_size( Functor::SHARED_TEAM_COUNT ) +
                              Functor::shared_int_array_type::shmem_size( 3 * team_size );

    team_exec = p_type(8192 / team_size, team_size );

    Kokkos::parallel_reduce( team_exec.set_scratch_size( 1, Kokkos::PerTeam( team_scratch_size ),
                                                         Kokkos::PerThread( thread_scratch_size ) ),
                             Functor(), result_type( & error_count ) );
    ASSERT_EQ( error_count, 0 );
  }
};

} // namespace

namespace Test {

template< class ExecSpace >
KOKKOS_INLINE_FUNCTION
int test_team_mulit_level_scratch_loop_body( const typename Kokkos::TeamPolicy<ExecSpace>::member_type& team ) {
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > a_team1( team.team_scratch( 0 ), 128 );
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > a_thread1( team.thread_scratch( 0 ), 16 );
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > a_team2( team.team_scratch( 0 ), 128 );
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > a_thread2( team.thread_scratch( 0 ), 16 );

  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > b_team1( team.team_scratch( 1 ), 128000 );
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > b_thread1( team.thread_scratch( 1 ), 16000 );
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > b_team2( team.team_scratch( 1 ), 128000 );
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > b_thread2( team.thread_scratch( 1 ), 16000 );

  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > a_team3( team.team_scratch( 0 ), 128 );
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > a_thread3( team.thread_scratch( 0 ), 16 );
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > b_team3( team.team_scratch( 1 ), 128000 );
  Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > b_thread3( team.thread_scratch( 1 ), 16000 );

  // The explicit types for 0 and 128 are here to test TeamThreadRange accepting different
  // types for begin and end.
  Kokkos::parallel_for( Kokkos::TeamThreadRange( team, int( 0 ), unsigned( 128 ) ), [&] ( const int & i )
  {
    a_team1( i ) = 1000000 + i + team.league_rank() * 100000;
    a_team2( i ) = 2000000 + i + team.league_rank() * 100000;
    a_team3( i ) = 3000000 + i + team.league_rank() * 100000;
  });
  team.team_barrier();

  Kokkos::parallel_for( Kokkos::ThreadVectorRange( team, 16 ), [&] ( const int & i )
  {
    a_thread1( i ) = 1000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000;
    a_thread2( i ) = 2000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000;
    a_thread3( i ) = 3000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000;
  });

  Kokkos::parallel_for( Kokkos::TeamThreadRange( team, 0, 128000 ), [&] ( const int & i )
  {
    b_team1( i ) = 1000000 + i + team.league_rank() * 100000;
    b_team2( i ) = 2000000 + i + team.league_rank() * 100000;
    b_team3( i ) = 3000000 + i + team.league_rank() * 100000;
  });
  team.team_barrier();

  Kokkos::parallel_for( Kokkos::ThreadVectorRange( team, 16000 ), [&] ( const int & i )
  {
    b_thread1( i ) = 1000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000;
    b_thread2( i ) = 2000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000;
    b_thread3( i ) = 3000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000;
  });

  team.team_barrier();

  int error = 0;
  Kokkos::parallel_for( Kokkos::TeamThreadRange( team, 0, 128 ), [&] ( const int & i )
  {
    if ( a_team1( i ) != 1000000 + i + team.league_rank() * 100000 ) error++;
    if ( a_team2( i ) != 2000000 + i + team.league_rank() * 100000 ) error++;
    if ( a_team3( i ) != 3000000 + i + team.league_rank() * 100000 ) error++;
  });
  team.team_barrier();

  Kokkos::parallel_for( Kokkos::ThreadVectorRange( team, 16 ), [&] ( const int & i )
  {
    if ( a_thread1( i ) != 1000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000 ) error++;
    if ( a_thread2( i ) != 2000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000 ) error++;
    if ( a_thread3( i ) != 3000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000 ) error++;
  });

  Kokkos::parallel_for( Kokkos::TeamThreadRange( team, 0, 128000 ), [&] ( const int & i )
  {
    if ( b_team1( i ) != 1000000 + i + team.league_rank() * 100000 ) error++;
    if ( b_team2( i ) != 2000000 + i + team.league_rank() * 100000 ) error++;
    if ( b_team3( i ) != 3000000 + i + team.league_rank() * 100000 ) error++;
  });
  team.team_barrier();

  Kokkos::parallel_for( Kokkos::ThreadVectorRange( team, 16000 ), [&] ( const int & i )
  {
    if ( b_thread1( i ) != 1000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000 ) error++;
    if ( b_thread2( i ) != 2000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000 ) error++;
    if ( b_thread3( i ) != 3000000 + 100000 * team.team_rank() + 16 - i + team.league_rank() * 100000 ) error++;
  });

  return error;
}

struct TagReduce {};
struct TagFor {};

template< class ExecSpace, class ScheduleType >
struct ClassNoShmemSizeFunction {
  typedef typename Kokkos::TeamPolicy< ExecSpace, ScheduleType >::member_type member_type;

  Kokkos::View< int, ExecSpace, Kokkos::MemoryTraits<Kokkos::Atomic> > errors;

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagFor &, const member_type & team ) const {
    int error = test_team_mulit_level_scratch_loop_body< ExecSpace >( team );
    errors() += error;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const TagReduce &, const member_type & team, int & error ) const {
    error += test_team_mulit_level_scratch_loop_body< ExecSpace >( team );
  }

  void run() {
    Kokkos::View< int, ExecSpace > d_errors = Kokkos::View< int, ExecSpace >( "Errors" );
    errors = d_errors;

    const int per_team0 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 128 );
    const int per_thread0 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 16 );

    const int per_team1 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 128000 );
    const int per_thread1 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 16000 );

    int team_size = 8;
    if(team_size > ExecSpace::concurrency())
      team_size = ExecSpace::concurrency();
    {
      Kokkos::TeamPolicy< TagFor, ExecSpace, ScheduleType > policy( 10, team_size, 16 );

      Kokkos::parallel_for( policy.set_scratch_size( 0, Kokkos::PerTeam( per_team0 ), Kokkos::PerThread( per_thread0 ) ).set_scratch_size( 1, Kokkos::PerTeam( per_team1 ), Kokkos::PerThread( per_thread1 ) ), *this );
      Kokkos::fence();

      typename Kokkos::View< int, ExecSpace >::HostMirror h_errors = Kokkos::create_mirror_view( d_errors );
      Kokkos::deep_copy( h_errors, d_errors );
      ASSERT_EQ( h_errors(), 0 );
    }

    {
      int error = 0;
      Kokkos::TeamPolicy< TagReduce, ExecSpace, ScheduleType > policy( 10, team_size, 16 );

      Kokkos::parallel_reduce( policy.set_scratch_size( 0, Kokkos::PerTeam( per_team0 ), Kokkos::PerThread( per_thread0 ) ).set_scratch_size( 1, Kokkos::PerTeam( per_team1 ), Kokkos::PerThread( per_thread1 ) ), *this, error );
      Kokkos::fence();

      ASSERT_EQ( error, 0 );
    }
  };
};

template< class ExecSpace, class ScheduleType >
struct ClassWithShmemSizeFunction {
  typedef typename Kokkos::TeamPolicy< ExecSpace, ScheduleType >::member_type member_type;

  Kokkos::View< int, ExecSpace, Kokkos::MemoryTraits<Kokkos::Atomic> > errors;

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagFor &, const member_type & team ) const {
    int error = test_team_mulit_level_scratch_loop_body< ExecSpace >( team );
    errors() += error;
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const TagReduce &, const member_type & team, int & error ) const {
    error += test_team_mulit_level_scratch_loop_body< ExecSpace >( team );
  }

  void run() {
    Kokkos::View< int, ExecSpace > d_errors = Kokkos::View< int, ExecSpace >( "Errors" );
    errors = d_errors;

    const int per_team1 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 128000 );
    const int per_thread1 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 16000 );

    int team_size = 8;
    if(team_size > ExecSpace::concurrency())
      team_size = ExecSpace::concurrency();

    {
      Kokkos::TeamPolicy< TagFor, ExecSpace, ScheduleType > policy( 10, team_size, 16 );

      Kokkos::parallel_for( policy.set_scratch_size( 1, Kokkos::PerTeam( per_team1 ),
                                                     Kokkos::PerThread( per_thread1 ) ),
                            *this );
      Kokkos::fence();

      typename Kokkos::View< int, ExecSpace >::HostMirror h_errors = Kokkos::create_mirror_view( d_errors );
      Kokkos::deep_copy( h_errors, d_errors );
      ASSERT_EQ( h_errors(), 0 );
    }

    {
      int error = 0;
      Kokkos::TeamPolicy< TagReduce, ExecSpace, ScheduleType > policy( 10, team_size, 16 );

      Kokkos::parallel_reduce( policy.set_scratch_size( 1, Kokkos::PerTeam( per_team1 ),
                                                        Kokkos::PerThread( per_thread1 ) ),
                               *this, error );
      Kokkos::fence();

      ASSERT_EQ( error, 0 );
    }
  };

  unsigned team_shmem_size( int team_size ) const {
    const int per_team0 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 128 );
    const int per_thread0 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 16 );
    return per_team0 + team_size * per_thread0;
  }
};

template< class ExecSpace, class ScheduleType >
void test_team_mulit_level_scratch_test_lambda() {
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
#if !defined(KOKKOS_ENABLE_CUDA) || ( 8000 <= CUDA_VERSION )
  Kokkos::View< int, ExecSpace, Kokkos::MemoryTraits<Kokkos::Atomic> > errors;
  Kokkos::View< int, ExecSpace > d_errors( "Errors" );
  errors = d_errors;

  const int per_team0 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 128 );
  const int per_thread0 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 16 );

  const int per_team1 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 128000 );
  const int per_thread1 = 3 * Kokkos::View< double*, ExecSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >::shmem_size( 16000 );

  int team_size = 8;
  if(team_size > ExecSpace::concurrency())
    team_size = ExecSpace::concurrency();

  Kokkos::TeamPolicy< ExecSpace, ScheduleType > policy( 10, team_size, 16 );

  Kokkos::parallel_for( policy.set_scratch_size( 0, Kokkos::PerTeam( per_team0 ), Kokkos::PerThread( per_thread0 ) ).set_scratch_size( 1, Kokkos::PerTeam( per_team1 ), Kokkos::PerThread( per_thread1 ) ),
                        KOKKOS_LAMBDA ( const typename Kokkos::TeamPolicy< ExecSpace >::member_type & team )
  {
    int error = test_team_mulit_level_scratch_loop_body< ExecSpace >( team );
    errors() += error;
  });
  Kokkos::fence();

  typename Kokkos::View< int, ExecSpace >::HostMirror h_errors = Kokkos::create_mirror_view( errors );
  Kokkos::deep_copy( h_errors, d_errors );
  ASSERT_EQ( h_errors(), 0 );

  int error = 0;
  Kokkos::parallel_reduce( policy.set_scratch_size( 0, Kokkos::PerTeam( per_team0 ), Kokkos::PerThread( per_thread0 ) ).set_scratch_size( 1, Kokkos::PerTeam( per_team1 ), Kokkos::PerThread( per_thread1 ) ),
                           KOKKOS_LAMBDA ( const typename Kokkos::TeamPolicy< ExecSpace >::member_type & team, int & count )
  {
    count += test_team_mulit_level_scratch_loop_body< ExecSpace >( team );
  }, error );
  ASSERT_EQ( error, 0 );
  Kokkos::fence();
#endif
#endif
}

} // namespace Test

namespace {

template< class ExecSpace, class ScheduleType >
struct TestMultiLevelScratchTeam {
  TestMultiLevelScratchTeam() { run(); }

  void run()
  {
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
    Test::test_team_mulit_level_scratch_test_lambda< ExecSpace, ScheduleType >();
#endif
    Test::ClassNoShmemSizeFunction< ExecSpace, ScheduleType > c1;
    c1.run();

    Test::ClassWithShmemSizeFunction< ExecSpace, ScheduleType > c2;
    c2.run();
  }
};

} // namespace

namespace Test {

template< class ExecSpace >
struct TestShmemSize {
  TestShmemSize() { run(); }

  void run()
  {
    typedef Kokkos::View< long***, ExecSpace > view_type;

    size_t d1 = 5;
    size_t d2 = 6;
    size_t d3 = 7;

    size_t size = view_type::shmem_size( d1, d2, d3 );

    ASSERT_EQ( size, (d1 * d2 * d3 + 1)* sizeof( long ) );

    test_layout_stride();
  }

  void test_layout_stride()
  {
    int rank = 3;
    int order[3] = {2, 0, 1};
    int extents[3] = {100, 10, 3};
    auto s1 = Kokkos::View<double***, Kokkos::LayoutStride, ExecSpace>::shmem_size(Kokkos::LayoutStride::order_dimensions(rank, order, extents));
    auto s2 = Kokkos::View<double***, Kokkos::LayoutRight, ExecSpace>::shmem_size(extents[0], extents[1], extents[2]);
    ASSERT_EQ(s1, s2);
  }
};

} // namespace Test

/*--------------------------------------------------------------------------*/

namespace Test {

namespace {

template< class ExecSpace, class ScheduleType >
struct TestTeamBroadcast {
  typedef typename Kokkos::TeamPolicy< ScheduleType,  ExecSpace >::member_type team_member;

  TestTeamBroadcast( const size_t league_size ) {}

  struct BroadcastTag {};

  typedef long value_type;

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member &teamMember, value_type &update ) const
  {
    int lid = teamMember.league_rank();
    int tid = teamMember.team_rank();
    int ts  = teamMember.team_size();

    value_type parUpdate = 0;
    value_type value     = tid * 3 + 1;
	
    teamMember.team_broadcast(value, lid%ts); 

    Kokkos::parallel_reduce( Kokkos::TeamThreadRange( teamMember, ts ), [&] ( const int j, value_type &teamUpdate ) {
      teamUpdate += value;
    }, parUpdate );

    if ( teamMember.team_rank() == 0 ) update += parUpdate;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const BroadcastTag &, const team_member &teamMember, value_type &update ) const
  {
    int lid = teamMember.league_rank();
    int tid = teamMember.team_rank();
    int ts  = teamMember.team_size();

    value_type parUpdate = 0;
    value_type value     = tid * 3 + 1;

    teamMember.team_broadcast([&] (value_type & var) { var*=2; }, value, lid%ts);
    
    Kokkos::parallel_reduce( Kokkos::TeamThreadRange( teamMember, ts ), [&] ( const int j, value_type &teamUpdate ) {
      teamUpdate += value;
    }, parUpdate );

    if ( teamMember.team_rank() == 0 ) update += parUpdate;
  }

  static void test_teambroadcast( const size_t league_size )
  {
    TestTeamBroadcast functor( league_size );

    typedef Kokkos::TeamPolicy< ScheduleType, ExecSpace > policy_type;
    typedef Kokkos::TeamPolicy< ScheduleType, ExecSpace, BroadcastTag > policy_type_f;

    const int team_size = policy_type_f(league_size,1).team_size_max( functor, Kokkos::ParallelReduceTag() ); //printf("team_size=%d\n",team_size);

    //team_broadcast with value
    long total = 0;

    Kokkos::parallel_reduce( policy_type( league_size, team_size ), functor, total );
    
    value_type expected_result = 0;
    for (unsigned int i=0; i<league_size; i++){
      value_type val  = ((i%team_size)*3+1)*team_size;
      expected_result+= val;
    }
    ASSERT_EQ( size_t( expected_result ), size_t( total ) ); //printf("team_broadcast with value -- expected_result=%d, total=%d\n",expected_result, total);

    //team_broadcast with funtion object
    total = 0;

    Kokkos::parallel_reduce( policy_type_f( league_size, team_size ), functor, total );

    expected_result = 0;
    for (unsigned int i=0; i<league_size; i++){
      value_type val  = ((i%team_size)*3+1)*2*team_size;
      expected_result+= val;
    }
    ASSERT_EQ( size_t( expected_result ), size_t( total ) ); //printf("team_broadcast with funtion object -- expected_result=%d, total=%d\n",expected_result, total);
  }
};

template<class ExecSpace>
struct TestScratchAlignment {
  struct TestScalar {
    double x,y,z;
  };
  TestScratchAlignment() {
    test(true);
    test(false);
  }
  typedef Kokkos::View<TestScalar*,typename ExecSpace::scratch_memory_space> ScratchView;
  typedef Kokkos::View<int*,typename ExecSpace::scratch_memory_space> ScratchViewInt;
  void test(bool allocate_small) {
    int shmem_size = ScratchView::shmem_size(11);
    if(allocate_small) shmem_size += ScratchViewInt::shmem_size(1);
    Kokkos::parallel_for(Kokkos::TeamPolicy<ExecSpace>(1,1).set_scratch_size(0,Kokkos::PerTeam(shmem_size)),
     KOKKOS_LAMBDA (const typename Kokkos::TeamPolicy<ExecSpace>::member_type& team) {
     if(allocate_small) ScratchViewInt p(team.team_scratch(0),1);
     ScratchView a(team.team_scratch(0),11);
     if(ptrdiff_t(a.data())%sizeof(TestScalar)!=0)
       Kokkos::abort("Error: invalid scratch view alignment\n");
    });
    Kokkos::fence();
  }
};

} // namespace

} // namespace Test

/*--------------------------------------------------------------------------*/
