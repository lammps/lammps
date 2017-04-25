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

#include <serial/TestSerial.hpp>

namespace Test {

TEST_F( serial , mdrange_for )
{
  TestMDRange_2D< Kokkos::Serial >::test_for2( 100, 100 );
  TestMDRange_3D< Kokkos::Serial >::test_for3( 100, 10, 100 );
  TestMDRange_4D< Kokkos::Serial >::test_for4( 100, 10, 10, 10 );
  TestMDRange_5D< Kokkos::Serial >::test_for5( 100, 10, 10, 10, 5 );
  TestMDRange_6D< Kokkos::Serial >::test_for6( 10, 10, 10, 10, 5, 5 );
}

TEST_F( serial , mdrange_reduce )
{
  TestMDRange_2D< Kokkos::Serial >::test_reduce2( 100, 100 );
  TestMDRange_3D< Kokkos::Serial >::test_reduce3( 100, 10, 100 );
}

TEST_F( serial, policy_construction )
{
  TestRangePolicyConstruction< Kokkos::Serial >();
  TestTeamPolicyConstruction< Kokkos::Serial >();
}

TEST_F( serial, range_tag )
{
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Static> >::test_for( 0 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 0 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Static> >::test_scan( 0 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 0 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 0 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Dynamic> >::test_scan( 0 );

  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Static> >::test_for( 1000 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 1000 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Static> >::test_scan( 1000 );

  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 1001 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 1001 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Dynamic> >::test_scan( 1001 );
  TestRange< Kokkos::Serial, Kokkos::Schedule<Kokkos::Dynamic> >::test_dynamic_policy( 1000 );
}

//----------------------------------------------------------------------------

TEST_F( serial, compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Serial >() ) );
}

//----------------------------------------------------------------------------

TEST_F( serial, memory_pool )
{
  bool val = TestMemoryPool::test_mempool< Kokkos::Serial >( 128, 128000000 );
  ASSERT_TRUE( val );

  TestMemoryPool::test_mempool2< Kokkos::Serial >( 64, 4, 1000000, 2000000 );

  TestMemoryPool::test_memory_exhaustion< Kokkos::Serial >();
}

//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_TASKDAG )

TEST_F( serial, task_fib )
{
  for ( int i = 0; i < 25; ++i ) {
    TestTaskScheduler::TestFib< Kokkos::Serial >::run( i );
  }
}

TEST_F( serial, task_depend )
{
  for ( int i = 0; i < 25; ++i ) {
    TestTaskScheduler::TestTaskDependence< Kokkos::Serial >::run( i );
  }
}

TEST_F( serial, task_team )
{
  TestTaskScheduler::TestTaskTeam< Kokkos::Serial >::run( 1000 );
  //TestTaskScheduler::TestTaskTeamValue< Kokkos::Serial >::run( 1000 ); // Put back after testing.
}

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */

//----------------------------------------------------------------------------

#if defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL )
TEST_F( serial, cxx11 )
{
  if ( std::is_same< Kokkos::DefaultExecutionSpace, Kokkos::Serial >::value ) {
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >( 1 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >( 2 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >( 3 ) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >( 4 ) ) );
  }
}
#endif

TEST_F( serial, tile_layout )
{
  TestTile::test< Kokkos::Serial, 1, 1 >( 1, 1 );
  TestTile::test< Kokkos::Serial, 1, 1 >( 2, 3 );
  TestTile::test< Kokkos::Serial, 1, 1 >( 9, 10 );

  TestTile::test< Kokkos::Serial, 2, 2 >( 1, 1 );
  TestTile::test< Kokkos::Serial, 2, 2 >( 2, 3 );
  TestTile::test< Kokkos::Serial, 2, 2 >( 4, 4 );
  TestTile::test< Kokkos::Serial, 2, 2 >( 9, 9 );

  TestTile::test< Kokkos::Serial, 2, 4 >( 9, 9 );
  TestTile::test< Kokkos::Serial, 4, 2 >( 9, 9 );

  TestTile::test< Kokkos::Serial, 4, 4 >( 1, 1 );
  TestTile::test< Kokkos::Serial, 4, 4 >( 4, 4 );
  TestTile::test< Kokkos::Serial, 4, 4 >( 9, 9 );
  TestTile::test< Kokkos::Serial, 4, 4 >( 9, 11 );

  TestTile::test< Kokkos::Serial, 8, 8 >( 1, 1 );
  TestTile::test< Kokkos::Serial, 8, 8 >( 4, 4 );
  TestTile::test< Kokkos::Serial, 8, 8 >( 9, 9 );
  TestTile::test< Kokkos::Serial, 8, 8 >( 9, 11 );
}

} // namespace Test
