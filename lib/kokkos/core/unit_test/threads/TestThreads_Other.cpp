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

#include <threads/TestThreads.hpp>

namespace Test {

TEST_F( threads, init )
{
  ;
}

TEST_F( threads , mdrange_for ) {
  TestMDRange_2D< Kokkos::Threads >::test_for2( 100, 100 );
  TestMDRange_3D< Kokkos::Threads >::test_for3( 100, 10, 100 );
  TestMDRange_4D< Kokkos::Threads >::test_for4( 100, 10, 10, 10 );
  TestMDRange_5D< Kokkos::Threads >::test_for5( 100, 10, 10, 10, 5 );
  TestMDRange_6D< Kokkos::Threads >::test_for6( 10, 10, 10, 10, 5, 5 );
}

TEST_F( threads , mdrange_reduce ) {
  TestMDRange_2D< Kokkos::Threads >::test_reduce2( 100, 100 );
  TestMDRange_3D< Kokkos::Threads >::test_reduce3( 100, 10, 100 );
}

TEST_F( threads, policy_construction )
{
  TestRangePolicyConstruction< Kokkos::Threads >();
  TestTeamPolicyConstruction< Kokkos::Threads >();
}

TEST_F( threads, range_tag )
{
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_for( 0 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 0 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_scan( 0 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 0 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 0 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_scan( 0 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_dynamic_policy( 0 );

  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_for( 2 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 2 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_scan( 2 );

  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 3 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 3 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_scan( 3 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_dynamic_policy( 3 );

  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_for( 1000 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 1000 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_scan( 1000 );

  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 1001 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 1001 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_scan( 1001 );
  TestRange< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_dynamic_policy( 1000 );
}

//----------------------------------------------------------------------------

TEST_F( threads, compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Threads >() ) );
}

//----------------------------------------------------------------------------

TEST_F( threads, memory_pool )
{
  bool val = TestMemoryPool::test_mempool< Kokkos::Threads >( 128, 128000000 );
  ASSERT_TRUE( val );

  TestMemoryPool::test_mempool2< Kokkos::Threads >( 64, 4, 1000000, 2000000 );

  TestMemoryPool::test_memory_exhaustion< Kokkos::Threads >();
}

//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_TASKDAG )
/*
TEST_F( threads, task_fib )
{
  for ( int i = 0; i < 25; ++i ) {
    TestTaskScheduler::TestFib< Kokkos::Threads >::run( i );
  }
}

TEST_F( threads, task_depend )
{
  for ( int i = 0; i < 25; ++i ) {
    TestTaskScheduler::TestTaskDependence< Kokkos::Threads >::run( i );
  }
}

TEST_F( threads, task_team )
{
  TestTaskScheduler::TestTaskTeam< Kokkos::Threads >::run( 1000 );
  //TestTaskScheduler::TestTaskTeamValue< Kokkos::Threads >::run( 1000 ); // Put back after testing.
}
*/
#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */

//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS )
TEST_F( threads, cxx11 )
{
  if ( std::is_same< Kokkos::DefaultExecutionSpace, Kokkos::Threads >::value ) {
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >( 1 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >( 2 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >( 3 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >( 4 ) ) );
  }
}
#endif

TEST_F( threads, tile_layout )
{
  TestTile::test< Kokkos::Threads, 1, 1 >( 1, 1 );
  TestTile::test< Kokkos::Threads, 1, 1 >( 2, 3 );
  TestTile::test< Kokkos::Threads, 1, 1 >( 9, 10 );

  TestTile::test< Kokkos::Threads, 2, 2 >( 1, 1 );
  TestTile::test< Kokkos::Threads, 2, 2 >( 2, 3 );
  TestTile::test< Kokkos::Threads, 2, 2 >( 4, 4 );
  TestTile::test< Kokkos::Threads, 2, 2 >( 9, 9 );

  TestTile::test< Kokkos::Threads, 2, 4 >( 9, 9 );
  TestTile::test< Kokkos::Threads, 4, 2 >( 9, 9 );

  TestTile::test< Kokkos::Threads, 4, 4 >( 1, 1 );
  TestTile::test< Kokkos::Threads, 4, 4 >( 4, 4 );
  TestTile::test< Kokkos::Threads, 4, 4 >( 9, 9 );
  TestTile::test< Kokkos::Threads, 4, 4 >( 9, 11 );

  TestTile::test< Kokkos::Threads, 8, 8 >( 1, 1 );
  TestTile::test< Kokkos::Threads, 8, 8 >( 4, 4 );
  TestTile::test< Kokkos::Threads, 8, 8 >( 9, 9 );
  TestTile::test< Kokkos::Threads, 8, 8 >( 9, 11 );
}

TEST_F( threads, dispatch )
{
  const int repeat = 100;
  for ( int i = 0; i < repeat; ++i ) {
    for ( int j = 0; j < repeat; ++j ) {
      Kokkos::parallel_for( Kokkos::RangePolicy< Kokkos::Threads >( 0, j )
                          , KOKKOS_LAMBDA( int ) {} );
    }
  }
}

} // namespace Test
