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

#include <qthreads/TestQthreads.hpp>

namespace Test {

TEST_F( qthreads, init )
{
  ;
}

TEST_F( qthreads, md_range )
{
#if 0
  TestMDRange_2D< Kokkos::Qthreads >::test_for2( 100, 100 );
  TestMDRange_3D< Kokkos::Qthreads >::test_for3( 100, 100, 100 );
#endif
}

TEST_F( qthreads, policy_construction )
{
#if 0
  TestRangePolicyConstruction< Kokkos::Qthreads >();
  TestTeamPolicyConstruction< Kokkos::Qthreads >();
#endif
}

TEST_F( qthreads, range_tag )
{
#if 0
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_for( 0 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 0 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_scan( 0 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 0 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 0 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_scan( 0 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_dynamic_policy( 0 );

  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_for( 2 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 2 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_scan( 2 );

  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 3 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 3 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_scan( 3 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_dynamic_policy( 3 );

  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_for( 1000 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 1000 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_scan( 1000 );

  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 1001 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 1001 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_scan( 1001 );
  TestRange< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_dynamic_policy( 1000 );
#endif
}

//----------------------------------------------------------------------------

TEST_F( qthreads, compiler_macros )
{
#if 0
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Qthreads >() ) );
#endif
}

//----------------------------------------------------------------------------

TEST_F( qthreads, memory_pool )
{
#if 0

#endif
}

//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_TASKDAG )

TEST_F( qthreads, task_fib )
{
#if 0
  const int N = 24 ; // 25 triggers tbd bug on Cuda/Pascal
  for ( int i = 0; i < N; ++i ) {
    TestTaskScheduler::TestFib< Kokkos::Qthreads >::run( i, ( i + 1 ) * ( i + 1 ) * 10000 );
  }
#endif
}

TEST_F( qthreads, task_depend )
{
#if 0
  for ( int i = 0; i < 25; ++i ) {
    TestTaskScheduler::TestTaskDependence< Kokkos::Qthreads >::run( i );
  }
#endif
}

TEST_F( qthreads, task_team )
{
#if 0
  TestTaskScheduler::TestTaskTeam< Kokkos::Qthreads >::run( 1000 );
  //TestTaskScheduler::TestTaskTeamValue< Kokkos::Qthreads >::run( 1000 ); // Put back after testing.
#endif
}

#endif // #if defined( KOKKOS_ENABLE_TASKDAG )

//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_QTHREADS )

TEST_F( qthreads, cxx11 )
{
#if 0
  if ( std::is_same< Kokkos::DefaultExecutionSpace, Kokkos::Qthreads >::value ) {
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Qthreads >( 1 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Qthreads >( 2 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Qthreads >( 3 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Qthreads >( 4 ) ) );
  }
#endif
}

#endif

TEST_F( qthreads, tile_layout )
{
#if 0
  TestTile::test< Kokkos::Qthreads, 1, 1 >( 1, 1 );
  TestTile::test< Kokkos::Qthreads, 1, 1 >( 2, 3 );
  TestTile::test< Kokkos::Qthreads, 1, 1 >( 9, 10 );

  TestTile::test< Kokkos::Qthreads, 2, 2 >( 1, 1 );
  TestTile::test< Kokkos::Qthreads, 2, 2 >( 2, 3 );
  TestTile::test< Kokkos::Qthreads, 2, 2 >( 4, 4 );
  TestTile::test< Kokkos::Qthreads, 2, 2 >( 9, 9 );

  TestTile::test< Kokkos::Qthreads, 2, 4 >( 9, 9 );
  TestTile::test< Kokkos::Qthreads, 4, 2 >( 9, 9 );

  TestTile::test< Kokkos::Qthreads, 4, 4 >( 1, 1 );
  TestTile::test< Kokkos::Qthreads, 4, 4 >( 4, 4 );
  TestTile::test< Kokkos::Qthreads, 4, 4 >( 9, 9 );
  TestTile::test< Kokkos::Qthreads, 4, 4 >( 9, 11 );

  TestTile::test< Kokkos::Qthreads, 8, 8 >( 1, 1 );
  TestTile::test< Kokkos::Qthreads, 8, 8 >( 4, 4 );
  TestTile::test< Kokkos::Qthreads, 8, 8 >( 9, 9 );
  TestTile::test< Kokkos::Qthreads, 8, 8 >( 9, 11 );
#endif
}

TEST_F( qthreads, dispatch )
{
#if 0
  const int repeat = 100;
  for ( int i = 0; i < repeat; ++i ) {
    for ( int j = 0; j < repeat; ++j ) {
      Kokkos::parallel_for( Kokkos::RangePolicy< Kokkos::Qthreads >( 0, j )
                          , KOKKOS_LAMBDA( int ) {} );
    }
  }
#endif
}

} // namespace Test
