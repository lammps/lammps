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

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <limits>

#include <Kokkos_Core.hpp>

namespace Test {

namespace ReduceCombinatorical {

template< class Scalar, class Space = Kokkos::HostSpace >
struct AddPlus {
public:
  // Required.
  typedef AddPlus reducer;
  typedef Scalar value_type;

  typedef Kokkos::View< value_type, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_type;

private:
  result_view_type result;

public:
  AddPlus( value_type & result_ ) : result( &result_ ) {}

  // Required.
  KOKKOS_INLINE_FUNCTION
  void join( value_type & dest, const value_type & src ) const {
    dest += src + 1;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & dest, const volatile value_type & src ) const {
    dest += src + 1;
  }

  // Optional.
  KOKKOS_INLINE_FUNCTION
  void init( value_type & val )  const {
    val = value_type();
  }

  KOKKOS_INLINE_FUNCTION
  value_type& reference() const {
    return result();
  }

  KOKKOS_INLINE_FUNCTION
  result_view_type view() const {
    return result;
  }
};

template< int ISTEAM >
struct FunctorScalar;

template<>
struct FunctorScalar< 0 > {
  Kokkos::View< double > result;

  FunctorScalar( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }
};

template<>
struct FunctorScalar< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalar( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }
};

template< int ISTEAM >
struct FunctorScalarInit;

template<>
struct FunctorScalarInit< 0 > {
  Kokkos::View< double > result;

  FunctorScalarInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template<>
struct FunctorScalarInit< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template< int ISTEAM >
struct FunctorScalarFinal;

template<>
struct FunctorScalarFinal< 0 > {
  Kokkos::View<double> result;

  FunctorScalarFinal( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }
};

template<>
struct FunctorScalarFinal< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarFinal( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }
};

template< int ISTEAM >
struct FunctorScalarJoin;

template<>
struct FunctorScalarJoin< 0 > {
  Kokkos::View<double> result;

  FunctorScalarJoin( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }
};

template<>
struct FunctorScalarJoin< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarJoin( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }
};

template< int ISTEAM >
struct FunctorScalarJoinFinal;

template<>
struct FunctorScalarJoinFinal< 0 > {
  Kokkos::View< double > result;

  FunctorScalarJoinFinal( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }
};

template<>
struct FunctorScalarJoinFinal< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarJoinFinal( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }
};

template< int ISTEAM >
struct FunctorScalarJoinInit;

template<>
struct FunctorScalarJoinInit< 0 > {
  Kokkos::View< double > result;

  FunctorScalarJoinInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template<>
struct FunctorScalarJoinInit< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarJoinInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template< int ISTEAM >
struct FunctorScalarJoinFinalInit;

template<>
struct FunctorScalarJoinFinalInit< 0 > {
  Kokkos::View<double> result;

  FunctorScalarJoinFinalInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

template<>
struct FunctorScalarJoinFinalInit< 1 > {
  typedef Kokkos::TeamPolicy<>::member_type team_type;

  Kokkos::View< double > result;

  FunctorScalarJoinFinalInit( Kokkos::View< double > r ) : result( r ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_type & team, double & update ) const {
    update += 1.0 / team.team_size() * team.league_rank();
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double & dst, const volatile double & update ) const {
    dst += update;
  }

  KOKKOS_INLINE_FUNCTION
  void final( double & update ) const {
    result() = update;
  }

  KOKKOS_INLINE_FUNCTION
  void init( double & update ) const {
    update = 0.0;
  }
};

struct Functor1 {
  KOKKOS_INLINE_FUNCTION
  void operator()( const int & i, double & update ) const {
    update += i;
  }
};

struct Functor2 {
  typedef double value_type[];

  const unsigned value_count;

  Functor2( unsigned n ) : value_count( n ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned & i, double update[] ) const {
    for ( unsigned j = 0; j < value_count; j++ ) {
      update[j] += i;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( double dst[] ) const
  {
    for ( unsigned i = 0; i < value_count; ++i ) dst[i] = 0;
  }

  KOKKOS_INLINE_FUNCTION
  void join( volatile double dst[],
             const volatile double src[] ) const
  {
    for ( unsigned i = 0; i < value_count; ++i ) dst[i] += src[i];
  }
};

} // namespace ReduceCombinatorical

template< class ExecSpace = Kokkos::DefaultExecutionSpace >
struct TestReduceCombinatoricalInstantiation {
  template< class ... Args >
  static void CallParallelReduce( Args... args ) {
    Kokkos::parallel_reduce( args... );
  }

  template< class ... Args >
  static void AddReturnArgument( Args... args ) {
    Kokkos::View< double, Kokkos::HostSpace > result_view( "ResultView" );
    double expected_result = 1000.0 * 999.0 / 2.0;

    double value = 0;
    Kokkos::parallel_reduce( args..., value );
    ASSERT_EQ( expected_result, value );

    result_view() = 0;
    CallParallelReduce( args..., result_view );
    Kokkos::fence();
    ASSERT_EQ( expected_result, result_view() );

    value = 0;
    CallParallelReduce( args..., Kokkos::View< double, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> >( &value ) );
    Kokkos::fence();
    ASSERT_EQ( expected_result, value );

    result_view() = 0;
    const Kokkos::View< double, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > result_view_const_um = result_view;
    CallParallelReduce( args..., result_view_const_um );
    Kokkos::fence();
    ASSERT_EQ( expected_result, result_view_const_um() );

    value = 0;
    CallParallelReduce( args..., Test::ReduceCombinatorical::AddPlus< double >( value ) );
    if ( ( Kokkos::DefaultExecutionSpace::concurrency() > 1 ) && ( ExecSpace::concurrency() > 1 ) ) {
      ASSERT_TRUE( expected_result < value );
    }
    else if ( ( Kokkos::DefaultExecutionSpace::concurrency() > 1 ) || ( ExecSpace::concurrency() > 1 ) ) {
      ASSERT_TRUE( expected_result <= value );
    }
    else {
      ASSERT_EQ( expected_result, value );
    }

    value = 0;
    Test::ReduceCombinatorical::AddPlus< double > add( value );
    CallParallelReduce( args..., add );
    if ( ( Kokkos::DefaultExecutionSpace::concurrency() > 1 ) && ( ExecSpace::concurrency() > 1 ) ) {
      ASSERT_TRUE( expected_result < value );
    }
    else if ( ( Kokkos::DefaultExecutionSpace::concurrency() > 1 ) || ( ExecSpace::concurrency() > 1 ) ) {
      ASSERT_TRUE( expected_result <= value );
    }
    else {
      ASSERT_EQ( expected_result, value );
    }
  }

  template< class ... Args >
  static void AddLambdaRange( void*, Args... args ) {
    AddReturnArgument( args..., KOKKOS_LAMBDA ( const int & i, double & lsum ) {
      lsum += i;
    });
  }

  template< class ... Args >
  static void AddLambdaTeam( void*, Args... args ) {
    AddReturnArgument( args..., KOKKOS_LAMBDA ( const Kokkos::TeamPolicy<>::member_type & team, double & update ) {
      update += 1.0 / team.team_size() * team.league_rank();
    });
  }

  template< class ... Args >
  static void AddLambdaRange( Kokkos::InvalidType, Args... args ) {}

  template< class ... Args >
  static void AddLambdaTeam( Kokkos::InvalidType, Args... args ) {}

  template< int ISTEAM, class ... Args >
  static void AddFunctor( Args... args ) {
    Kokkos::View< double > result_view( "FunctorView" );
    auto h_r = Kokkos::create_mirror_view( result_view );
    Test::ReduceCombinatorical::FunctorScalar< ISTEAM > functor( result_view );
    double expected_result = 1000.0 * 999.0 / 2.0;

    AddReturnArgument( args..., functor );
    AddReturnArgument( args..., Test::ReduceCombinatorical::FunctorScalar< ISTEAM >( result_view ) );
    AddReturnArgument( args..., Test::ReduceCombinatorical::FunctorScalarInit< ISTEAM >( result_view ) );
    AddReturnArgument( args..., Test::ReduceCombinatorical::FunctorScalarJoin< ISTEAM >( result_view ) );
    AddReturnArgument( args..., Test::ReduceCombinatorical::FunctorScalarJoinInit< ISTEAM >( result_view ) );

    h_r() = 0;
    Kokkos::deep_copy( result_view, h_r );
    CallParallelReduce( args..., Test::ReduceCombinatorical::FunctorScalarFinal< ISTEAM >( result_view ) );
    Kokkos::fence();
    Kokkos::deep_copy( h_r, result_view );
    ASSERT_EQ( expected_result, h_r() );

    h_r() = 0;
    Kokkos::deep_copy( result_view, h_r );
    CallParallelReduce( args..., Test::ReduceCombinatorical::FunctorScalarJoinFinal< ISTEAM >( result_view ) );
    Kokkos::fence();
    Kokkos::deep_copy( h_r, result_view );
    ASSERT_EQ( expected_result, h_r() );

    h_r() = 0;
    Kokkos::deep_copy( result_view, h_r );
    CallParallelReduce( args..., Test::ReduceCombinatorical::FunctorScalarJoinFinalInit< ISTEAM >( result_view ) );
    Kokkos::fence();
    Kokkos::deep_copy( h_r, result_view );
    ASSERT_EQ( expected_result, h_r() );
  }

  template< class ... Args >
  static void AddFunctorLambdaRange( Args... args ) {
    AddFunctor< 0, Args... >( args... );
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
    AddLambdaRange( typename std::conditional< std::is_same<ExecSpace, Kokkos::DefaultExecutionSpace>::value, void*, Kokkos::InvalidType >::type(), args... );
#endif
  }

  template< class ... Args >
  static void AddFunctorLambdaTeam( Args... args ) {
    AddFunctor< 1, Args... >( args... );
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
    AddLambdaTeam( typename std::conditional< std::is_same<ExecSpace, Kokkos::DefaultExecutionSpace>::value, void*, Kokkos::InvalidType >::type(), args... );
#endif
  }

  template< class ... Args >
  static void AddPolicy_1( Args... args ) {
    int N = 1000;
    Kokkos::RangePolicy< ExecSpace > policy( 0, N );

    AddFunctorLambdaRange( args..., 1000 );
    AddFunctorLambdaRange( args..., N );
    AddFunctorLambdaRange( args..., policy );
  }

  template< class ... Args >
  static void AddPolicy_2( Args... args ) {
    int N = 1000;
    Kokkos::RangePolicy< ExecSpace > policy( 0, N );

    AddFunctorLambdaRange( args..., Kokkos::RangePolicy< ExecSpace >( 0, N ) );
    AddFunctorLambdaRange( args..., Kokkos::RangePolicy< ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >( 0, N ) );
    AddFunctorLambdaRange( args..., Kokkos::RangePolicy< ExecSpace, Kokkos::Schedule<Kokkos::Static> >( 0, N ).set_chunk_size( 10 ) );
    AddFunctorLambdaRange( args..., Kokkos::RangePolicy< ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >( 0, N ).set_chunk_size( 10 ) );

  }

  template< class ... Args >
  static void AddPolicy_3( Args... args ) {
    int N = 1000;
    Kokkos::RangePolicy< ExecSpace > policy( 0, N );

    AddFunctorLambdaTeam( args..., Kokkos::TeamPolicy< ExecSpace >( N, Kokkos::AUTO ) );
    AddFunctorLambdaTeam( args..., Kokkos::TeamPolicy< ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >( N, Kokkos::AUTO ) );
    AddFunctorLambdaTeam( args..., Kokkos::TeamPolicy< ExecSpace, Kokkos::Schedule<Kokkos::Static> >( N, Kokkos::AUTO ).set_chunk_size( 10 ) );
    AddFunctorLambdaTeam( args..., Kokkos::TeamPolicy< ExecSpace, Kokkos::Schedule<Kokkos::Dynamic> >( N, Kokkos::AUTO ).set_chunk_size( 10 ) );
  }

  static void execute_a1() {
    AddPolicy_1();
  }

  static void execute_b1() {
    std::string s( "Std::String" );
    AddPolicy_1( s.c_str() );
    AddPolicy_1( "Char Constant" );
  }

  static void execute_c1() {
    std::string s( "Std::String" );
    AddPolicy_1( s );
  }

  static void execute_a2() {
    AddPolicy_2();
  }

  static void execute_b2() {
    std::string s( "Std::String" );
    AddPolicy_2( s.c_str() );
    AddPolicy_2( "Char Constant" );
  }

  static void execute_c2() {
    std::string s( "Std::String" );
    AddPolicy_2( s );
  }

  static void execute_a3() {
    AddPolicy_1();
  }

  static void execute_b3() {
    std::string s( "Std::String" );
    AddPolicy_1( s.c_str() );
    AddPolicy_1( "Char Constant" );
  }

  static void execute_c3() {
    std::string s( "Std::String" );
    AddPolicy_1( s );
  }

};

} // namespace Test

