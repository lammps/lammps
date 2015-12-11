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

#include <stdio.h>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

/*--------------------------------------------------------------------------*/

namespace Test {
namespace {

template< class ExecSpace >
struct TestTeamPolicy {

  typedef typename Kokkos::TeamPolicy< ExecSpace >::member_type team_member ;
  typedef Kokkos::View<int**,ExecSpace> view_type ;

  view_type m_flags ;

  TestTeamPolicy( const size_t league_size )
    : m_flags( Kokkos::ViewAllocateWithoutInitializing("flags")
             , Kokkos::TeamPolicy< ExecSpace >::team_size_max( *this )
             , league_size )
    {}

  struct VerifyInitTag {};

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & member ) const
    {
      const int tid = member.team_rank() + member.team_size() * member.league_rank();

      m_flags( member.team_rank() , member.league_rank() ) = tid ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( const VerifyInitTag & , const team_member & member ) const
    {
      const int tid = member.team_rank() + member.team_size() * member.league_rank();

      if ( tid != m_flags( member.team_rank() , member.league_rank() ) ) {
        printf("TestTeamPolicy member(%d,%d) error %d != %d\n"
              , member.league_rank() , member.team_rank()
              , tid , m_flags( member.team_rank() , member.league_rank() ) );
      }
    }

  static void test_for( const size_t league_size )
    {
      TestTeamPolicy functor( league_size );

      const int team_size = Kokkos::TeamPolicy< ExecSpace >::team_size_max( functor );

      Kokkos::parallel_for( Kokkos::TeamPolicy< ExecSpace >( league_size , team_size ) , functor );
      Kokkos::parallel_for( Kokkos::TeamPolicy< ExecSpace , VerifyInitTag >( league_size , team_size ) , functor );
    }

  struct ReduceTag {};

  typedef long value_type ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & member , value_type & update ) const
    {
      update += member.team_rank() + member.team_size() * member.league_rank();
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( const ReduceTag & , const team_member & member , value_type & update ) const
    {
      update += 1 + member.team_rank() + member.team_size() * member.league_rank();
    }

  static void test_reduce( const size_t league_size )
    {
      TestTeamPolicy functor( league_size );

      const int team_size = Kokkos::TeamPolicy< ExecSpace >::team_size_max( functor );
      const long N = team_size * league_size ;

      long total = 0 ;

      Kokkos::parallel_reduce( Kokkos::TeamPolicy< ExecSpace >( league_size , team_size ) , functor , total );
      ASSERT_EQ( size_t((N-1)*(N))/2 , size_t(total) );

      Kokkos::parallel_reduce( Kokkos::TeamPolicy< ExecSpace , ReduceTag >( league_size , team_size ) , functor , total );
      ASSERT_EQ( (size_t(N)*size_t(N+1))/2 , size_t(total) );
    }
};

}
}

/*--------------------------------------------------------------------------*/

namespace Test {

template< typename ScalarType , class DeviceType >
class ReduceTeamFunctor
{
public:
  typedef DeviceType execution_space ;
  typedef Kokkos::TeamPolicy< execution_space >  policy_type ;
  typedef typename execution_space::size_type        size_type ;

  struct value_type {
    ScalarType value[3] ;
  };

  const size_type nwork ;

  ReduceTeamFunctor( const size_type & arg_nwork ) : nwork( arg_nwork ) {}

  ReduceTeamFunctor( const ReduceTeamFunctor & rhs )
    : nwork( rhs.nwork ) {}

  KOKKOS_INLINE_FUNCTION
  void init( value_type & dst ) const
  {
    dst.value[0] = 0 ;
    dst.value[1] = 0 ;
    dst.value[2] = 0 ;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dst ,
             const volatile value_type & src ) const
  {
    dst.value[0] += src.value[0] ;
    dst.value[1] += src.value[1] ;
    dst.value[2] += src.value[2] ;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const typename policy_type::member_type ind , value_type & dst ) const
  {
    const int thread_rank = ind.team_rank() + ind.team_size() * ind.league_rank();
    const int thread_size = ind.team_size() * ind.league_size();
    const int chunk = ( nwork + thread_size - 1 ) / thread_size ;

    size_type iwork = chunk * thread_rank ;
    const size_type iwork_end = iwork + chunk < nwork ? iwork + chunk : nwork ;

    for ( ; iwork < iwork_end ; ++iwork ) {
      dst.value[0] += 1 ;
      dst.value[1] += iwork + 1 ;
      dst.value[2] += nwork - iwork ;
    }
  }
};

} // namespace Test

namespace {

template< typename ScalarType , class DeviceType >
class TestReduceTeam
{
public:
  typedef DeviceType    execution_space ;
  typedef Kokkos::TeamPolicy< execution_space >  policy_type ;
  typedef typename execution_space::size_type    size_type ;

  //------------------------------------

  TestReduceTeam( const size_type & nwork )
  {
    run_test(nwork);
  }

  void run_test( const size_type & nwork )
  {
    typedef Test::ReduceTeamFunctor< ScalarType , execution_space > functor_type ;
    typedef typename functor_type::value_type value_type ;
    typedef Kokkos::View< value_type, Kokkos::HostSpace, Kokkos::MemoryUnmanaged > result_type ;

    enum { Count = 3 };
    enum { Repeat = 100 };

    value_type result[ Repeat ];

    const unsigned long nw   = nwork ;
    const unsigned long nsum = nw % 2 ? nw * (( nw + 1 )/2 )
                                      : (nw/2) * ( nw + 1 );

    const unsigned team_size   = policy_type::team_size_recommended( functor_type(nwork) );
    const unsigned league_size = ( nwork + team_size - 1 ) / team_size ;

    policy_type team_exec( league_size , team_size );

    for ( unsigned i = 0 ; i < Repeat ; ++i ) {
      result_type tmp( & result[i] );
      Kokkos::parallel_reduce( team_exec , functor_type(nwork) , tmp );
    }

    execution_space::fence();

    for ( unsigned i = 0 ; i < Repeat ; ++i ) {
      for ( unsigned j = 0 ; j < Count ; ++j ) {
        const unsigned long correct = 0 == j % 3 ? nw : nsum ;
        ASSERT_EQ( (ScalarType) correct , result[i].value[j] );
      }
    }
  }
};

}

/*--------------------------------------------------------------------------*/

namespace Test {

template< class DeviceType >
class ScanTeamFunctor
{
public:
  typedef DeviceType  execution_space ;
  typedef Kokkos::TeamPolicy< execution_space >  policy_type ;

  typedef long int    value_type ;
  Kokkos::View< value_type , execution_space > accum ;
  Kokkos::View< value_type , execution_space > total ;

  ScanTeamFunctor() : accum("accum"), total("total") {}

  KOKKOS_INLINE_FUNCTION
  void init( value_type & error ) const { error = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void join( value_type volatile & error ,
             value_type volatile const & input ) const
    { if ( input ) error = 1 ; }

  struct JoinMax {
    typedef long int value_type ;
    KOKKOS_INLINE_FUNCTION
    void join( value_type volatile & dst
             , value_type volatile const & input ) const
      { if ( dst < input ) dst = input ; }
  };

  KOKKOS_INLINE_FUNCTION
  void operator()( const typename policy_type::member_type ind , value_type & error ) const
  {
    if ( 0 == ind.league_rank() && 0 == ind.team_rank() ) {
      const long int thread_count = ind.league_size() * ind.team_size();
      total() = ( thread_count * ( thread_count + 1 ) ) / 2 ;
    }

    // Team max:
    const int long m = ind.team_reduce( (long int) ( ind.league_rank() + ind.team_rank() ) , JoinMax() );

    if ( m != ind.league_rank() + ( ind.team_size() - 1 ) ) {
      printf("ScanTeamFunctor[%d.%d of %d.%d] reduce_max_answer(%ld) != reduce_max(%ld)\n"
            , ind.league_rank(), ind.team_rank()
            , ind.league_size(), ind.team_size()
            , (long int)(ind.league_rank() + ( ind.team_size() - 1 )) , m );
    }

    // Scan:
    const long int answer =
      ( ind.league_rank() + 1 ) * ind.team_rank() +
      ( ind.team_rank() * ( ind.team_rank() + 1 ) ) / 2 ;

    const long int result =
      ind.team_scan( ind.league_rank() + 1 + ind.team_rank() + 1 );

    const long int result2 =
      ind.team_scan( ind.league_rank() + 1 + ind.team_rank() + 1 );

    if ( answer != result || answer != result2 ) {
      printf("ScanTeamFunctor[%d.%d of %d.%d] answer(%ld) != scan_first(%ld) or scan_second(%ld)\n",
             ind.league_rank(), ind.team_rank(),
             ind.league_size(), ind.team_size(),
             answer,result,result2);
      error = 1 ;
    }

    const long int thread_rank = ind.team_rank() +
                                 ind.team_size() * ind.league_rank();
    ind.team_scan( 1 + thread_rank , accum.ptr_on_device() );
  }
};

template< class DeviceType >
class TestScanTeam
{
public:
  typedef DeviceType  execution_space ;
  typedef long int    value_type ;

  typedef Kokkos::TeamPolicy< execution_space > policy_type ;
  typedef Test::ScanTeamFunctor<DeviceType> functor_type ;

  //------------------------------------

  TestScanTeam( const size_t nteam )
  {
    run_test(nteam);
  }

  void run_test( const size_t nteam )
  {
    typedef Kokkos::View< long int , Kokkos::HostSpace , Kokkos::MemoryUnmanaged >  result_type ;

    const unsigned REPEAT = 100000 ;
    const unsigned Repeat = ( REPEAT + nteam - 1 ) / nteam ;

    functor_type functor ;

    policy_type team_exec( nteam , policy_type::team_size_max( functor ) );

    for ( unsigned i = 0 ; i < Repeat ; ++i ) {
      long int accum = 0 ;
      long int total = 0 ;
      long int error = 0 ;
      Kokkos::deep_copy( functor.accum , total );
      Kokkos::parallel_reduce( team_exec , functor , result_type( & error ) );
      DeviceType::fence();
      Kokkos::deep_copy( accum , functor.accum );
      Kokkos::deep_copy( total , functor.total );

      ASSERT_EQ( error , 0 );
      ASSERT_EQ( total , accum );
    }

    execution_space::fence();
  }
};

} // namespace Test

/*--------------------------------------------------------------------------*/

namespace Test {

template< class ExecSpace >
struct SharedTeamFunctor {

  typedef ExecSpace  execution_space ;
  typedef int        value_type ;
  typedef Kokkos::TeamPolicy< execution_space >  policy_type ;

  enum { SHARED_COUNT = 1000 };

  typedef typename ExecSpace::scratch_memory_space shmem_space ;

  // tbd: MemoryUnmanaged should be the default for shared memory space
  typedef Kokkos::View<int*,shmem_space,Kokkos::MemoryUnmanaged> shared_int_array_type ;

  // Tell how much shared memory will be required by this functor:
  inline
  unsigned team_shmem_size( int /* team_size */ ) const
  {
    return shared_int_array_type::shmem_size( SHARED_COUNT ) +
           shared_int_array_type::shmem_size( SHARED_COUNT );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const typename policy_type::member_type & ind , value_type & update ) const
  {
    const shared_int_array_type shared_A( ind.team_shmem() , SHARED_COUNT );
    const shared_int_array_type shared_B( ind.team_shmem() , SHARED_COUNT );

    if ((shared_A.ptr_on_device () == NULL && SHARED_COUNT > 0) ||
        (shared_B.ptr_on_device () == NULL && SHARED_COUNT > 0)) {
      printf ("Failed to allocate shared memory of size %lu\n",
              static_cast<unsigned long> (SHARED_COUNT));
      ++update; // failure to allocate is an error
    }
    else {
      for ( int i = ind.team_rank() ; i < SHARED_COUNT ; i += ind.team_size() ) {
        shared_A[i] = i + ind.league_rank();
        shared_B[i] = 2 * i + ind.league_rank();
      }

      ind.team_barrier();

      if ( ind.team_rank() + 1 == ind.team_size() ) {
        for ( int i = 0 ; i < SHARED_COUNT ; ++i ) {
          if ( shared_A[i] != i + ind.league_rank() ) {
            ++update ;
          }
          if ( shared_B[i] != 2 * i + ind.league_rank() ) {
            ++update ;
          }
        }
      }
    }
  }
};

}

namespace {

template< class ExecSpace >
struct TestSharedTeam {

  TestSharedTeam()
  { run(); }

  void run()
  {
    typedef Test::SharedTeamFunctor<ExecSpace> Functor ;
    typedef Kokkos::View< typename Functor::value_type , Kokkos::HostSpace , Kokkos::MemoryUnmanaged >  result_type ;

    const size_t team_size = Kokkos::TeamPolicy< ExecSpace >::team_size_max( Functor() );

    Kokkos::TeamPolicy< ExecSpace > team_exec( 8192 / team_size , team_size );

    typename Functor::value_type error_count = 0 ;

    Kokkos::parallel_reduce( team_exec , Functor() , result_type( & error_count ) );

    ASSERT_EQ( error_count , 0 );
  }
};

#if defined (KOKKOS_HAVE_CXX11_DISPATCH_LAMBDA)

template< class ExecSpace >
struct TestLambdaSharedTeam {

  TestLambdaSharedTeam()
  { run(); }

  void run()
  {
    typedef Test::SharedTeamFunctor<ExecSpace> Functor ;
    typedef Kokkos::View< typename Functor::value_type , Kokkos::HostSpace , Kokkos::MemoryUnmanaged >  result_type ;
    typedef typename ExecSpace::scratch_memory_space shmem_space ;

    // tbd: MemoryUnmanaged should be the default for shared memory space
    typedef Kokkos::View<int*,shmem_space,Kokkos::MemoryUnmanaged> shared_int_array_type ;

    const int SHARED_COUNT = 1000;
    int team_size = 1;
#ifdef KOKKOS_HAVE_CUDA
    if(std::is_same<ExecSpace,Kokkos::Cuda>::value)
      team_size = 128;
#endif
    Kokkos::TeamPolicy< ExecSpace > team_exec( 8192 / team_size , team_size ,
        Kokkos::Experimental::TeamScratchRequest<shmem_space>(SHARED_COUNT*2*sizeof(int)));

    typename Functor::value_type error_count = 0 ;

    Kokkos::parallel_reduce( team_exec , KOKKOS_LAMBDA
        ( const typename Kokkos::TeamPolicy< ExecSpace >::member_type & ind , int & update ) {

      const shared_int_array_type shared_A( ind.team_shmem() , SHARED_COUNT );
      const shared_int_array_type shared_B( ind.team_shmem() , SHARED_COUNT );

      if ((shared_A.ptr_on_device () == NULL && SHARED_COUNT > 0) ||
          (shared_B.ptr_on_device () == NULL && SHARED_COUNT > 0)) {
        printf ("Failed to allocate shared memory of size %lu\n",
                static_cast<unsigned long> (SHARED_COUNT));
        ++update; // failure to allocate is an error
      } else {
        for ( int i = ind.team_rank() ; i < SHARED_COUNT ; i += ind.team_size() ) {
          shared_A[i] = i + ind.league_rank();
          shared_B[i] = 2 * i + ind.league_rank();
        }

        ind.team_barrier();

        if ( ind.team_rank() + 1 == ind.team_size() ) {
          for ( int i = 0 ; i < SHARED_COUNT ; ++i ) {
            if ( shared_A[i] != i + ind.league_rank() ) {
              ++update ;
            }
            if ( shared_B[i] != 2 * i + ind.league_rank() ) {
              ++update ;
            }
          }
        }
      }
    }, result_type( & error_count ) );

    ASSERT_EQ( error_count , 0 );
  }
};

#endif
}

/*--------------------------------------------------------------------------*/
