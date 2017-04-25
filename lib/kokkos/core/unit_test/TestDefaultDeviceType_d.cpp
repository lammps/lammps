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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#if !defined( KOKKOS_ENABLE_CUDA ) || defined( __CUDACC__ )

#include <TestAtomic.hpp>
#include <TestViewAPI.hpp>
#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestTeam.hpp>
#include <TestAggregate.hpp>
#include <TestCompilerMacros.hpp>
#include <TestCXX11.hpp>
#include <TestTeamVector.hpp>
#include <TestUtilities.hpp>

namespace Test {

class defaultdevicetype : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    Kokkos::initialize();
  }

  static void TearDownTestCase()
  {
    Kokkos::finalize();
  }
};

TEST_F( defaultdevicetype, test_utilities )
{
  test_utilities();
}

TEST_F( defaultdevicetype, long_reduce )
{
  TestReduce< long, Kokkos::DefaultExecutionSpace >( 100000 );
}

TEST_F( defaultdevicetype, double_reduce )
{
  TestReduce< double, Kokkos::DefaultExecutionSpace >( 100000 );
}

TEST_F( defaultdevicetype, long_reduce_dynamic )
{
  TestReduceDynamic< long, Kokkos::DefaultExecutionSpace >( 100000 );
}

TEST_F( defaultdevicetype, double_reduce_dynamic )
{
  TestReduceDynamic< double, Kokkos::DefaultExecutionSpace >( 100000 );
}

TEST_F( defaultdevicetype, long_reduce_dynamic_view )
{
  TestReduceDynamicView< long, Kokkos::DefaultExecutionSpace >( 100000 );
}

TEST_F( defaultdevicetype, atomics )
{
  const int loop_count = 1e4;

  ASSERT_TRUE( ( TestAtomic::Loop< int, Kokkos::DefaultExecutionSpace >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< int, Kokkos::DefaultExecutionSpace >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< int, Kokkos::DefaultExecutionSpace >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< unsigned int, Kokkos::DefaultExecutionSpace >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< unsigned int, Kokkos::DefaultExecutionSpace >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< unsigned int, Kokkos::DefaultExecutionSpace >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< long int, Kokkos::DefaultExecutionSpace >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< long int, Kokkos::DefaultExecutionSpace >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< long int, Kokkos::DefaultExecutionSpace >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< unsigned long int, Kokkos::DefaultExecutionSpace >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< unsigned long int, Kokkos::DefaultExecutionSpace >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< unsigned long int, Kokkos::DefaultExecutionSpace >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< long long int, Kokkos::DefaultExecutionSpace >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< long long int, Kokkos::DefaultExecutionSpace >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< long long int, Kokkos::DefaultExecutionSpace >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< double, Kokkos::DefaultExecutionSpace >( loop_count, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< double, Kokkos::DefaultExecutionSpace >( loop_count, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< double, Kokkos::DefaultExecutionSpace >( loop_count, 3 ) ) );

  ASSERT_TRUE( ( TestAtomic::Loop< float, Kokkos::DefaultExecutionSpace >( 100, 1 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< float, Kokkos::DefaultExecutionSpace >( 100, 2 ) ) );
  ASSERT_TRUE( ( TestAtomic::Loop< float, Kokkos::DefaultExecutionSpace >( 100, 3 ) ) );
}

/*TEST_F( defaultdevicetype, view_remap )
{
  enum { N0 = 3, N1 = 2, N2 = 8, N3 = 9 };

  typedef Kokkos::View< double*[N1][N2][N3],
                        Kokkos::LayoutRight,
                        Kokkos::DefaultExecutionSpace > output_type;

  typedef Kokkos::View< int**[N2][N3],
                        Kokkos::LayoutLeft,
                        Kokkos::DefaultExecutionSpace > input_type;

  typedef Kokkos::View< int*[N0][N2][N3],
                        Kokkos::LayoutLeft,
                        Kokkos::DefaultExecutionSpace > diff_type;

  output_type output( "output", N0 );
  input_type  input ( "input", N0, N1 );
  diff_type   diff  ( "diff", N0 );

  int value = 0;
  for ( size_t i3 = 0; i3 < N3; ++i3 ) {
    for ( size_t i2 = 0; i2 < N2; ++i2 ) {
      for ( size_t i1 = 0; i1 < N1; ++i1 ) {
        for ( size_t i0 = 0; i0 < N0; ++i0 ) {
          input( i0, i1, i2, i3 ) = ++value;
        }
      }
    }
  }

  // Kokkos::deep_copy( diff, input ); // Throw with incompatible shape.
  Kokkos::deep_copy( output, input );

  value = 0;
  for ( size_t i3 = 0; i3 < N3; ++i3 ) {
    for ( size_t i2 = 0; i2 < N2; ++i2 ) {
      for ( size_t i1 = 0; i1 < N1; ++i1 ) {
        for ( size_t i0 = 0; i0 < N0; ++i0 ) {
          ++value;
          ASSERT_EQ( value, ( (int) output( i0, i1, i2, i3 ) ) );
        }
      }
    }
  }
}*/

TEST_F( defaultdevicetype, view_aggregate )
{
  TestViewAggregate< Kokkos::DefaultExecutionSpace >();
}

TEST_F( defaultdevicetype, scan )
{
  TestScan< Kokkos::DefaultExecutionSpace >::test_range( 1, 1000 );
  TestScan< Kokkos::DefaultExecutionSpace >( 1000000 );
  TestScan< Kokkos::DefaultExecutionSpace >( 10000000 );
  Kokkos::DefaultExecutionSpace::fence();
}

TEST_F( defaultdevicetype, compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::DefaultExecutionSpace >() ) );
}

TEST_F( defaultdevicetype, cxx11 )
{
  ASSERT_TRUE( ( TestCXX11::Test< Kokkos::DefaultExecutionSpace >( 1 ) ) );
  ASSERT_TRUE( ( TestCXX11::Test< Kokkos::DefaultExecutionSpace >( 2 ) ) );
  ASSERT_TRUE( ( TestCXX11::Test< Kokkos::DefaultExecutionSpace >( 3 ) ) );
  ASSERT_TRUE( ( TestCXX11::Test< Kokkos::DefaultExecutionSpace >( 4 ) ) );
}

#if !defined(KOKKOS_CUDA_CLANG_WORKAROUND) && !defined(KOKKOS_ARCH_PASCAL)
TEST_F( defaultdevicetype, team_vector )
{
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::DefaultExecutionSpace >( 0 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::DefaultExecutionSpace >( 1 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::DefaultExecutionSpace >( 2 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::DefaultExecutionSpace >( 3 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::DefaultExecutionSpace >( 4 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::DefaultExecutionSpace >( 5 ) ) );
}
#endif

TEST_F( defaultdevicetype, malloc )
{
  int* data = (int*) Kokkos::kokkos_malloc( 100 * sizeof( int ) );
  ASSERT_NO_THROW( data = (int*) Kokkos::kokkos_realloc( data, 120 * sizeof( int ) ) );
  Kokkos::kokkos_free( data );

  int* data2 = (int*) Kokkos::kokkos_malloc( 0 );
  ASSERT_TRUE( data2 == NULL );
  Kokkos::kokkos_free( data2 );
}

} // namespace Test

#endif
