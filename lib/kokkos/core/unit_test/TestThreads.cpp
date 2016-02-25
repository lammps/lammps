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

#include <Kokkos_Macros.hpp>

#if defined( KOKKOS_HAVE_PTHREAD )

#include <Kokkos_Core.hpp>

#include <Threads/Kokkos_Threads_TaskPolicy.hpp>

//----------------------------------------------------------------------------

#include <TestSharedAlloc.hpp>
#include <TestViewMapping.hpp>

#include <TestViewImpl.hpp>

#include <TestViewAPI.hpp>
#include <TestViewSubview.hpp>
#include <TestViewOfClass.hpp>
#include <TestAtomic.hpp>

#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestRange.hpp>
#include <TestTeam.hpp>
#include <TestAggregate.hpp>
#include <TestAggregateReduction.hpp>
#include <TestCompilerMacros.hpp>
#include <TestCXX11.hpp>
#include <TestCXX11Deduction.hpp>
#include <TestTeamVector.hpp>
#include <TestMemorySpaceTracking.hpp>
#include <TestTemplateMetaFunctions.hpp>

#include <TestTaskPolicy.hpp>

namespace Test {

class threads : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    // Finalize without initialize is a no-op:
    Kokkos::Threads::finalize();

    const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
    const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
    const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

    unsigned threads_count = 0 ;

    // Initialize and finalize with no threads:
    Kokkos::Threads::initialize( 1u );
    Kokkos::Threads::finalize();

    threads_count = std::max( 1u , numa_count )
                  * std::max( 2u , cores_per_numa * threads_per_core );

    Kokkos::Threads::initialize( threads_count );
    Kokkos::Threads::finalize();

    
    threads_count = std::max( 1u , numa_count * 2 )
                  * std::max( 2u , ( cores_per_numa * threads_per_core ) / 2 );

    Kokkos::Threads::initialize( threads_count );
    Kokkos::Threads::finalize();

    // Quick attempt to verify thread start/terminate don't have race condition:
    threads_count = std::max( 1u , numa_count )
                  * std::max( 2u , ( cores_per_numa * threads_per_core ) / 2 );
    for ( unsigned i = 0 ; i < 10 ; ++i ) {
      Kokkos::Threads::initialize( threads_count );
      Kokkos::Threads::sleep();
      Kokkos::Threads::wake();
      Kokkos::Threads::finalize();
    }

    Kokkos::Threads::initialize( threads_count );
    Kokkos::Threads::print_configuration( std::cout , true /* detailed */ );
  }

  static void TearDownTestCase()
  {
    Kokkos::Threads::finalize();
  }
};

TEST_F( threads , init ) {
  ;
}

TEST_F( threads , dispatch )
{
  const int repeat = 100 ;
  for ( int i = 0 ; i < repeat ; ++i ) {
  for ( int j = 0 ; j < repeat ; ++j ) {
    Kokkos::parallel_for( Kokkos::RangePolicy< Kokkos::Threads >(0,j)
                        , KOKKOS_LAMBDA( int ) {} );
  }}
}

TEST_F( threads , impl_shared_alloc ) {
  test_shared_alloc< Kokkos::HostSpace , Kokkos::Threads >();
}

TEST_F( threads , impl_view_mapping ) {
  test_view_mapping< Kokkos::Threads >();
  test_view_mapping_subview< Kokkos::Threads >();
  test_view_mapping_operator< Kokkos::Threads >();
  TestViewMappingAtomic< Kokkos::Threads >::run();
}


TEST_F( threads, view_impl) {
  test_view_impl< Kokkos::Threads >();
}

TEST_F( threads, view_api) {
  TestViewAPI< double , Kokkos::Threads >();
}

TEST_F( threads , view_nested_view )
{
  ::Test::view_nested_view< Kokkos::Threads >();
}

TEST_F( threads, view_subview_auto_1d_left ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutLeft,Kokkos::Threads >();
}

TEST_F( threads, view_subview_auto_1d_right ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutRight,Kokkos::Threads >();
}

TEST_F( threads, view_subview_auto_1d_stride ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutStride,Kokkos::Threads >();
}

TEST_F( threads, view_subview_assign_strided ) {
  TestViewSubview::test_1d_strided_assignment< Kokkos::Threads >();
}

TEST_F( threads, view_subview_left_0 ) {
  TestViewSubview::test_left_0< Kokkos::Threads >();
}

TEST_F( threads, view_subview_left_1 ) {
  TestViewSubview::test_left_1< Kokkos::Threads >();
}

TEST_F( threads, view_subview_left_2 ) {
  TestViewSubview::test_left_2< Kokkos::Threads >();
}

TEST_F( threads, view_subview_left_3 ) {
  TestViewSubview::test_left_3< Kokkos::Threads >();
}

TEST_F( threads, view_subview_right_0 ) {
  TestViewSubview::test_right_0< Kokkos::Threads >();
}

TEST_F( threads, view_subview_right_1 ) {
  TestViewSubview::test_right_1< Kokkos::Threads >();
}

TEST_F( threads, view_subview_right_3 ) {
  TestViewSubview::test_right_3< Kokkos::Threads >();
}


TEST_F( threads, view_aggregate ) {
  TestViewAggregate< Kokkos::Threads >();
  TestViewAggregateReduction< Kokkos::Threads >();
}

TEST_F( threads , range_tag )
{
  TestRange< Kokkos::Threads >::test_for(1000);
  TestRange< Kokkos::Threads >::test_reduce(1000);
  TestRange< Kokkos::Threads >::test_scan(1000);
}

TEST_F( threads , team_tag )
{
  TestTeamPolicy< Kokkos::Threads >::test_for(1000);
  TestTeamPolicy< Kokkos::Threads >::test_reduce(1000);
}

TEST_F( threads, long_reduce) {
  TestReduce< long ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, double_reduce) {
  TestReduce< double ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, team_long_reduce) {
  TestReduceTeam< long ,   Kokkos::Threads >( 100000 );
}

TEST_F( threads, team_double_reduce) {
  TestReduceTeam< double ,   Kokkos::Threads >( 100000 );
}

TEST_F( threads, long_reduce_dynamic ) {
  TestReduceDynamic< long ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, double_reduce_dynamic ) {
  TestReduceDynamic< double ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, long_reduce_dynamic_view ) {
  TestReduceDynamicView< long ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, team_shared_request) {
  TestSharedTeam< Kokkos::Threads >();
}

#if defined(KOKKOS_HAVE_CXX11_DISPATCH_LAMBDA) && !defined(KOKKOS_HAVE_CUDA)
TEST_F( threads, team_lambda_shared_request) {
  TestLambdaSharedTeam< Kokkos::Threads >();
}
#endif

TEST_F( threads , view_remap )
{
  enum { N0 = 3 , N1 = 2 , N2 = 8 , N3 = 9 };

  typedef Kokkos::View< double*[N1][N2][N3] ,
                             Kokkos::LayoutRight ,
                             Kokkos::Threads > output_type ;

  typedef Kokkos::View< int**[N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Threads > input_type ;

  typedef Kokkos::View< int*[N0][N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Threads > diff_type ;

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

TEST_F( threads , atomics )
{
  const int loop_count = 1e6 ;

  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Threads>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Threads>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Threads>(100,3) ) );

#if defined( KOKKOS_ENABLE_ASM )
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Threads>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Threads>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Threads>(100,3) ) );
#endif

  ASSERT_TRUE( ( TestAtomic::Loop<TestAtomic::SuperScalar<3>, Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<TestAtomic::SuperScalar<3>, Kokkos::Threads>(loop_count,2) ) );
}

//----------------------------------------------------------------------------

#if 0
TEST_F( threads , scan_small )
{
  typedef TestScan< Kokkos::Threads , Kokkos::Impl::ThreadsExecUseScanSmall > TestScanFunctor ;
  for ( int i = 0 ; i < 1000 ; ++i ) {
    TestScanFunctor( 10 );
    TestScanFunctor( 10000 );
  }
  TestScanFunctor( 1000000 );
  TestScanFunctor( 10000000 );

  Kokkos::Threads::fence();
}
#endif

TEST_F( threads , scan )
{
  TestScan< Kokkos::Threads >::test_range( 1 , 1000 );
  TestScan< Kokkos::Threads >( 1000000 );
  TestScan< Kokkos::Threads >( 10000000 );
  Kokkos::Threads::fence();
}

//----------------------------------------------------------------------------

TEST_F( threads , team_scan )
{
  TestScanTeam< Kokkos::Threads >( 10 );
  TestScanTeam< Kokkos::Threads >( 10000 );
}

//----------------------------------------------------------------------------

TEST_F( threads , compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Threads >() ) );
}

TEST_F( threads , memory_space )
{
  TestMemorySpace< Kokkos::Threads >();
}

//----------------------------------------------------------------------------

TEST_F( threads , template_meta_functions )
{
  TestTemplateMetaFunctions<int, Kokkos::Threads >();
}

//----------------------------------------------------------------------------

#if defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS )
TEST_F( threads , cxx11 )
{
  if ( Kokkos::Impl::is_same< Kokkos::DefaultExecutionSpace , Kokkos::Threads >::value ) {
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >(1) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >(2) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >(3) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >(4) ) );
  }
}

TEST_F( threads , reduction_deduction )
{
  TestCXX11::test_reduction_deduction< Kokkos::Threads >();
}
#endif /* #if defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) */

TEST_F( threads , team_vector )
{
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(0) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(1) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(2) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(3) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(4) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(5) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(6) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(7) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(8) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(9) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >(10) ) );
}

TEST_F( threads , task_policy )
{
  TestTaskPolicy::test_task_dep< Kokkos::Threads >( 10 );
  for ( long i = 0 ; i < 25 ; ++i ) TestTaskPolicy::test_fib< Kokkos::Threads >(i);
  for ( long i = 0 ; i < 35 ; ++i ) TestTaskPolicy::test_fib2< Kokkos::Threads >(i);
}

TEST_F( threads , task_team )
{
  TestTaskPolicy::test_task_team< Kokkos::Threads >(1000);
}

} // namespace Test

#endif /* #if defined( KOKKOS_HAVE_PTHREAD ) */
