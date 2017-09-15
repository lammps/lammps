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

#include <Kokkos_Core.hpp>

namespace TestCXX11 {

template< class DeviceType >
struct FunctorAddTest {
  typedef Kokkos::View< double**, DeviceType > view_type;
  typedef DeviceType execution_space;
  typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member;

  view_type a_, b_;

  FunctorAddTest( view_type & a, view_type & b ) : a_( a ), b_( b ) {}

  KOKKOS_INLINE_FUNCTION
  void operator() ( const int& i ) const {
    b_( i, 0 ) = a_( i, 1 ) + a_( i, 2 );
    b_( i, 1 ) = a_( i, 0 ) - a_( i, 3 );
    b_( i, 2 ) = a_( i, 4 ) + a_( i, 0 );
    b_( i, 3 ) = a_( i, 2 ) - a_( i, 1 );
    b_( i, 4 ) = a_( i, 3 ) + a_( i, 4 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const team_member & dev ) const {
    const int begin = dev.league_rank() * 4;
    const int end   = begin + 4;
    for ( int i = begin + dev.team_rank(); i < end; i += dev.team_size() ) {
      b_( i, 0 ) = a_( i, 1 ) + a_( i, 2 );
      b_( i, 1 ) = a_( i, 0 ) - a_( i, 3 );
      b_( i, 2 ) = a_( i, 4 ) + a_( i, 0 );
      b_( i, 3 ) = a_( i, 2 ) - a_( i, 1 );
      b_( i, 4 ) = a_( i, 3 ) + a_( i, 4 );
    }
  }
};

template< class DeviceType, bool PWRTest >
double AddTestFunctor() {
  typedef Kokkos::TeamPolicy< DeviceType > policy_type;

  Kokkos::View< double**, DeviceType > a( "A", 100, 5 );
  Kokkos::View< double**, DeviceType > b( "B", 100, 5 );
  typename Kokkos::View< double**, DeviceType >::HostMirror h_a = Kokkos::create_mirror_view( a );
  typename Kokkos::View< double**, DeviceType >::HostMirror h_b = Kokkos::create_mirror_view( b );

  for ( int i = 0; i < 100; i++ ) {
    for  ( int j = 0; j < 5; j++ ) {
       h_a( i, j ) = 0.1 * i / ( 1.1 * j + 1.0 ) + 0.5 * j;
    }
  }
  Kokkos::deep_copy( a, h_a );

  if ( PWRTest == false ) {
    Kokkos::parallel_for( 100, FunctorAddTest< DeviceType >( a, b ) );
  }
  else {
    Kokkos::parallel_for( policy_type( 25, Kokkos::AUTO ), FunctorAddTest< DeviceType >( a, b ) );
  }
  Kokkos::deep_copy( h_b, b );

  double result = 0;
  for ( int i = 0; i < 100; i++ ) {
    for ( int j = 0; j < 5; j++ ) {
      result += h_b( i, j );
    }
  }

  return result;
}

#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
template< class DeviceType, bool PWRTest >
double AddTestLambda() {
  Kokkos::View< double**, DeviceType > a( "A", 100, 5 );
  Kokkos::View< double**, DeviceType > b( "B", 100, 5 );
  typename Kokkos::View< double**, DeviceType >::HostMirror h_a = Kokkos::create_mirror_view( a );
  typename Kokkos::View< double**, DeviceType >::HostMirror h_b = Kokkos::create_mirror_view( b );

  for ( int i = 0; i < 100; i++ ) {
    for ( int j = 0; j < 5; j++ ) {
       h_a( i, j ) = 0.1 * i / ( 1.1 * j + 1.0 ) + 0.5 * j;
    }
  }
  Kokkos::deep_copy( a, h_a );

  if ( PWRTest == false ) {
    Kokkos::parallel_for( 100, KOKKOS_LAMBDA( const int & i ) {
      b( i, 0 ) = a( i, 1 ) + a( i, 2 );
      b( i, 1 ) = a( i, 0 ) - a( i, 3 );
      b( i, 2 ) = a( i, 4 ) + a( i, 0 );
      b( i, 3 ) = a( i, 2 ) - a( i, 1 );
      b( i, 4 ) = a( i, 3 ) + a( i, 4 );
    });
  }
  else {
    typedef Kokkos::TeamPolicy< DeviceType > policy_type;
    typedef typename policy_type::member_type team_member;

    policy_type policy( 25, Kokkos::AUTO );

    Kokkos::parallel_for( policy, KOKKOS_LAMBDA( const team_member & dev ) {
      const int begin = dev.league_rank() * 4;
      const int end   = begin + 4;
      for ( int i = begin + dev.team_rank(); i < end; i += dev.team_size() ) {
        b( i, 0 ) = a( i, 1 ) + a( i, 2 );
        b( i, 1 ) = a( i, 0 ) - a( i, 3 );
        b( i, 2 ) = a( i, 4 ) + a( i, 0 );
        b( i, 3 ) = a( i, 2 ) - a( i, 1 );
        b( i, 4 ) = a( i, 3 ) + a( i, 4 );
      }
    });
  }
  Kokkos::deep_copy( h_b, b );

  double result = 0;
  for ( int i = 0; i < 100; i++ ) {
    for ( int j = 0; j < 5; j++ ) {
      result += h_b( i, j );
    }
  }

  return result;
}
#else
template< class DeviceType, bool PWRTest >
double AddTestLambda() {
  return AddTestFunctor< DeviceType, PWRTest >();
}
#endif

template< class DeviceType >
struct FunctorReduceTest {
  typedef Kokkos::View< double**, DeviceType > view_type;
  typedef DeviceType execution_space;
  typedef double value_type;
  typedef typename Kokkos::TeamPolicy< execution_space >::member_type team_member;

  view_type a_;

  FunctorReduceTest( view_type & a ) : a_( a ) {}

  KOKKOS_INLINE_FUNCTION
  void operator() ( const int & i, value_type & sum ) const {
    sum += a_( i, 1 ) + a_( i, 2 );
    sum += a_( i, 0 ) - a_( i, 3 );
    sum += a_( i, 4 ) + a_( i, 0 );
    sum += a_( i, 2 ) - a_( i, 1 );
    sum += a_( i, 3 ) + a_( i, 4 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const team_member & dev, value_type & sum ) const {
    const int begin = dev.league_rank() * 4;
    const int end   = begin + 4;
    for ( int i = begin + dev.team_rank(); i < end; i += dev.team_size() ) {
      sum += a_( i, 1 ) + a_( i, 2 );
      sum += a_( i, 0 ) - a_( i, 3 );
      sum += a_( i, 4 ) + a_( i, 0 );
      sum += a_( i, 2 ) - a_( i, 1 );
      sum += a_( i, 3 ) + a_( i, 4 );
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type & update ) const { update = 0.0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & update, volatile value_type const & input ) const { update += input; }
};

template< class DeviceType, bool PWRTest >
double ReduceTestFunctor() {
  typedef Kokkos::TeamPolicy< DeviceType > policy_type;
  typedef Kokkos::View< double**, DeviceType > view_type;
  typedef Kokkos::View< double, typename view_type::host_mirror_space, Kokkos::MemoryUnmanaged > unmanaged_result;

  view_type a( "A", 100, 5 );
  typename view_type::HostMirror h_a = Kokkos::create_mirror_view( a );

  for ( int i = 0; i < 100; i++ ) {
    for ( int j = 0; j < 5; j++ ) {
       h_a( i, j ) = 0.1 * i / ( 1.1 * j + 1.0 ) + 0.5 * j;
    }
  }
  Kokkos::deep_copy( a, h_a );

  double result = 0.0;
  if ( PWRTest == false ) {
    Kokkos::parallel_reduce( 100, FunctorReduceTest< DeviceType >( a ), unmanaged_result( & result ) );
  }
  else {
    Kokkos::parallel_reduce( policy_type( 25, Kokkos::AUTO ), FunctorReduceTest< DeviceType >( a ), unmanaged_result( & result ) );
  }

  return result;
}

#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
template< class DeviceType, bool PWRTest >
double ReduceTestLambda() {
  typedef Kokkos::TeamPolicy< DeviceType > policy_type;
  typedef Kokkos::View< double**, DeviceType > view_type;
  typedef Kokkos::View< double, typename view_type::host_mirror_space, Kokkos::MemoryUnmanaged > unmanaged_result;

  view_type a( "A", 100, 5 );
  typename view_type::HostMirror h_a = Kokkos::create_mirror_view( a );

  for ( int i = 0; i < 100; i++ ) {
    for ( int j = 0; j < 5; j++ ) {
       h_a( i, j ) = 0.1 * i / ( 1.1 * j + 1.0 ) + 0.5 * j;
    }
  }
  Kokkos::deep_copy( a, h_a );

  double result = 0.0;

  if ( PWRTest == false ) {
    Kokkos::parallel_reduce( 100, KOKKOS_LAMBDA( const int & i, double & sum ) {
      sum += a( i, 1 ) + a( i, 2 );
      sum += a( i, 0 ) - a( i, 3 );
      sum += a( i, 4 ) + a( i, 0 );
      sum += a( i, 2 ) - a( i, 1 );
      sum += a( i, 3 ) + a( i, 4 );
    }, unmanaged_result( & result ) );
  }
  else {
    typedef typename policy_type::member_type team_member;
    Kokkos::parallel_reduce( policy_type( 25, Kokkos::AUTO ), KOKKOS_LAMBDA( const team_member & dev, double & sum ) {
      const int begin = dev.league_rank() * 4;
      const int end   = begin + 4;
      for ( int i = begin + dev.team_rank(); i < end; i += dev.team_size() ) {
        sum += a( i, 1 ) + a( i, 2 );
        sum += a( i, 0 ) - a( i, 3 );
        sum += a( i, 4 ) + a( i, 0 );
        sum += a( i, 2 ) - a( i, 1 );
        sum += a( i, 3 ) + a( i, 4 );
      }
    }, unmanaged_result( & result ) );
  }

  return result;
}
#else
template< class DeviceType, bool PWRTest >
double ReduceTestLambda() {
  return ReduceTestFunctor< DeviceType, PWRTest >();
}
#endif

template< class DeviceType >
double TestVariantLambda( int test ) {
  switch ( test ) {
    case 1: return AddTestLambda< DeviceType, false >();
    case 2: return AddTestLambda< DeviceType, true >();
    case 3: return ReduceTestLambda< DeviceType, false >();
    case 4: return ReduceTestLambda< DeviceType, true >();
  }

  return 0;
}

template< class DeviceType >
double TestVariantFunctor( int test ) {
  switch ( test ) {
    case 1: return AddTestFunctor< DeviceType, false >();
    case 2: return AddTestFunctor< DeviceType, true >();
    case 3: return ReduceTestFunctor< DeviceType, false >();
    case 4: return ReduceTestFunctor< DeviceType, true >();
  }

  return 0;
}

template< class DeviceType >
bool Test( int test ) {
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
  double res_functor = TestVariantFunctor< DeviceType >( test );
  double res_lambda = TestVariantLambda< DeviceType >( test );

  char testnames[5][256] = { " "
                           , "AddTest", "AddTest TeamPolicy"
                           , "ReduceTest", "ReduceTest TeamPolicy"
                           };
  bool passed = true;

  if ( res_functor != res_lambda ) {
    passed = false;

    std::cout << "CXX11 ( test = '"
              << testnames[test] << "' FAILED : "
              << res_functor << " != " << res_lambda
              << std::endl;
  }

  return passed;
#else
  return true;
#endif
}

} // namespace TestCXX11

namespace Test {
TEST_F( TEST_CATEGORY, cxx11 )
{
  if ( std::is_same< Kokkos::DefaultExecutionSpace, TEST_EXECSPACE >::value ) {
    ASSERT_TRUE( ( TestCXX11::Test< TEST_EXECSPACE >( 1 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< TEST_EXECSPACE >( 2 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< TEST_EXECSPACE >( 3 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< TEST_EXECSPACE >( 4 ) ) );
  }
}

}

