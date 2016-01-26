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

#include <impl/Kokkos_ViewTileLeft.hpp>
#include <TestTile.hpp>

#include <impl/Kokkos_Serial_TaskPolicy.hpp>

//----------------------------------------------------------------------------

#include <TestSharedAlloc.hpp>
#include <TestViewMapping.hpp>

#include <TestViewImpl.hpp>

#include <TestViewAPI.hpp>
#include <TestViewOfClass.hpp>
#include <TestViewSubview.hpp>
#include <TestAtomic.hpp>
#include <TestRange.hpp>
#include <TestTeam.hpp>
#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestAggregate.hpp>
#include <TestAggregateReduction.hpp>
#include <TestCompilerMacros.hpp>
#include <TestTaskPolicy.hpp>
#include <TestCXX11.hpp>
#include <TestCXX11Deduction.hpp>
#include <TestTeamVector.hpp>
#include <TestMemorySpaceTracking.hpp>
#include <TestTemplateMetaFunctions.hpp>

namespace Test {

class serial : public ::testing::Test {
protected:
  static void SetUpTestCase()
    {
      Kokkos::HostSpace::execution_space::initialize();
    }
  static void TearDownTestCase()
    {
      Kokkos::HostSpace::execution_space::finalize();
    }
};

TEST_F( serial , impl_shared_alloc ) {
  test_shared_alloc< Kokkos::HostSpace , Kokkos::Serial >();
}

TEST_F( serial , impl_view_mapping ) {
  test_view_mapping< Kokkos::Serial >();
  test_view_mapping_subview< Kokkos::Serial >();
  test_view_mapping_operator< Kokkos::Serial >();
  TestViewMappingAtomic< Kokkos::Serial >::run();
}

TEST_F( serial, view_impl) {
  test_view_impl< Kokkos::Serial >();
}

TEST_F( serial, view_api) {
  TestViewAPI< double , Kokkos::Serial >();
}

TEST_F( serial , view_nested_view )
{
  ::Test::view_nested_view< Kokkos::Serial >();
}

TEST_F( serial, view_subview_auto_1d_left ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutLeft,Kokkos::Serial >();
}

TEST_F( serial, view_subview_auto_1d_right ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutRight,Kokkos::Serial >();
}

TEST_F( serial, view_subview_auto_1d_stride ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutStride,Kokkos::Serial >();
}

TEST_F( serial, view_subview_assign_strided ) {
  TestViewSubview::test_1d_strided_assignment< Kokkos::Serial >();
}

TEST_F( serial, view_subview_left_0 ) {
  TestViewSubview::test_left_0< Kokkos::Serial >();
}

TEST_F( serial, view_subview_left_1 ) {
  TestViewSubview::test_left_1< Kokkos::Serial >();
}

TEST_F( serial, view_subview_left_2 ) {
  TestViewSubview::test_left_2< Kokkos::Serial >();
}

TEST_F( serial, view_subview_left_3 ) {
  TestViewSubview::test_left_3< Kokkos::Serial >();
}

TEST_F( serial, view_subview_right_0 ) {
  TestViewSubview::test_right_0< Kokkos::Serial >();
}

TEST_F( serial, view_subview_right_1 ) {
  TestViewSubview::test_right_1< Kokkos::Serial >();
}

TEST_F( serial, view_subview_right_3 ) {
  TestViewSubview::test_right_3< Kokkos::Serial >();
}

TEST_F( serial , range_tag )
{
  TestRange< Kokkos::Serial >::test_for(1000);
  TestRange< Kokkos::Serial >::test_reduce(1000);
  TestRange< Kokkos::Serial >::test_scan(1000);
}

TEST_F( serial , team_tag )
{
  TestTeamPolicy< Kokkos::Serial >::test_for( 1000 );
  TestTeamPolicy< Kokkos::Serial >::test_reduce( 1000 );
}

TEST_F( serial, long_reduce) {
  TestReduce< long ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial, double_reduce) {
  TestReduce< double ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial, long_reduce_dynamic ) {
  TestReduceDynamic< long ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial, double_reduce_dynamic ) {
  TestReduceDynamic< double ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial, long_reduce_dynamic_view ) {
  TestReduceDynamicView< long ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial , scan )
{
  TestScan< Kokkos::Serial >::test_range( 1 , 1000 );
  TestScan< Kokkos::Serial >( 10 );
  TestScan< Kokkos::Serial >( 10000 );
}

TEST_F( serial , team_long_reduce) {
  TestReduceTeam< long ,   Kokkos::Serial >( 100000 );
}

TEST_F( serial , team_double_reduce) {
  TestReduceTeam< double ,   Kokkos::Serial >( 100000 );
}

TEST_F( serial , team_shared_request) {
  TestSharedTeam< Kokkos::Serial >();
}

#if defined(KOKKOS_HAVE_CXX11_DISPATCH_LAMBDA) && !defined(KOKKOS_HAVE_CUDA)
TEST_F( serial , team_lambda_shared_request) {
  TestLambdaSharedTeam< Kokkos::Serial >();
}
#endif

TEST_F( serial  , team_scan )
{
  TestScanTeam< Kokkos::Serial >( 10 );
  TestScanTeam< Kokkos::Serial >( 10000 );
}


TEST_F( serial , view_remap )
{
  enum { N0 = 3 , N1 = 2 , N2 = 8 , N3 = 9 };

  typedef Kokkos::View< double*[N1][N2][N3] ,
                             Kokkos::LayoutRight ,
                             Kokkos::Serial > output_type ;

  typedef Kokkos::View< int**[N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Serial > input_type ;

  typedef Kokkos::View< int*[N0][N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Serial > diff_type ;

  output_type output( "output" , N0 );
  input_type  input ( "input" , N0 , N1 );
  diff_type   diff  ( "diff" , N0 );

  int value = 0 ;
  for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
  for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
  for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
  for ( size_t i0 = 0 ; i0 < N0 ; ++i0 ) {
    input(i0,i1,i2,i3) = ++value ;
  }}}}

  // Kokkos::deep_copy( diff , input ); // throw with incompatible shape
  Kokkos::deep_copy( output , input );

  value = 0 ;
  for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
  for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
  for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
  for ( size_t i0 = 0 ; i0 < N0 ; ++i0 ) {
    ++value ;
    ASSERT_EQ( value , ((int) output(i0,i1,i2,i3) ) );
  }}}}
}

//----------------------------------------------------------------------------

TEST_F( serial , view_aggregate )
{
  TestViewAggregate< Kokkos::Serial >();
  TestViewAggregateReduction< Kokkos::Serial >();
}

//----------------------------------------------------------------------------

TEST_F( serial , atomics )
{
  const int loop_count = 1e6 ;

  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Serial>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Serial>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Serial>(100,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Serial>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Serial>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Serial>(100,3) ) );
}

//----------------------------------------------------------------------------

TEST_F( serial, tile_layout )
{
  TestTile::test< Kokkos::Serial , 1 , 1 >( 1 , 1 );
  TestTile::test< Kokkos::Serial , 1 , 1 >( 2 , 3 );
  TestTile::test< Kokkos::Serial , 1 , 1 >( 9 , 10 );

  TestTile::test< Kokkos::Serial , 2 , 2 >( 1 , 1 );
  TestTile::test< Kokkos::Serial , 2 , 2 >( 2 , 3 );
  TestTile::test< Kokkos::Serial , 2 , 2 >( 4 , 4 );
  TestTile::test< Kokkos::Serial , 2 , 2 >( 9 , 9 );

  TestTile::test< Kokkos::Serial , 2 , 4 >( 9 , 9 );
  TestTile::test< Kokkos::Serial , 4 , 2 >( 9 , 9 );

  TestTile::test< Kokkos::Serial , 4 , 4 >( 1 , 1 );
  TestTile::test< Kokkos::Serial , 4 , 4 >( 4 , 4 );
  TestTile::test< Kokkos::Serial , 4 , 4 >( 9 , 9 );
  TestTile::test< Kokkos::Serial , 4 , 4 >( 9 , 11 );

  TestTile::test< Kokkos::Serial , 8 , 8 >( 1 , 1 );
  TestTile::test< Kokkos::Serial , 8 , 8 >( 4 , 4 );
  TestTile::test< Kokkos::Serial , 8 , 8 >( 9 , 9 );
  TestTile::test< Kokkos::Serial , 8 , 8 >( 9 , 11 );
}

//----------------------------------------------------------------------------

TEST_F( serial , compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Serial >() ) );
}

//----------------------------------------------------------------------------

TEST_F( serial , memory_space )
{
  TestMemorySpace< Kokkos::Serial >();
}

//----------------------------------------------------------------------------

TEST_F( serial , task_policy )
{
  TestTaskPolicy::test_task_dep< Kokkos::Serial >( 10 );
  // TestTaskPolicy::test_norm2< Kokkos::Serial >( 1000 );
  // for ( long i = 0 ; i < 30 ; ++i ) TestTaskPolicy::test_fib< Kokkos::Serial >(i);
  // for ( long i = 0 ; i < 40 ; ++i ) TestTaskPolicy::test_fib2< Kokkos::Serial >(i);
  for ( long i = 0 ; i < 20 ; ++i ) TestTaskPolicy::test_fib< Kokkos::Serial >(i);
  for ( long i = 0 ; i < 25 ; ++i ) TestTaskPolicy::test_fib2< Kokkos::Serial >(i);
}

TEST_F( serial , task_team )
{
  TestTaskPolicy::test_task_team< Kokkos::Serial >(1000);
}

//----------------------------------------------------------------------------

TEST_F( serial , template_meta_functions )
{
  TestTemplateMetaFunctions<int, Kokkos::Serial >();
}

//----------------------------------------------------------------------------

#if defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
TEST_F( serial , cxx11 )
{
  if ( Kokkos::Impl::is_same< Kokkos::DefaultExecutionSpace , Kokkos::Serial >::value ) {
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >(1) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >(2) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >(3) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >(4) ) );
  }
}
#endif

TEST_F( serial , reduction_deduction )
{
  TestCXX11::test_reduction_deduction< Kokkos::Serial >();
}

TEST_F( serial , team_vector )
{
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(0) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(1) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(2) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(3) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(4) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(5) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(6) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(7) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(8) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(9) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(10) ) );
}

} // namespace test

