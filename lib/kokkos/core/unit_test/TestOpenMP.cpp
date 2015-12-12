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

//----------------------------------------------------------------------------

#include <TestViewImpl.hpp>
#include <TestAtomic.hpp>

#include <TestViewAPI.hpp>
#include <TestViewSubview.hpp>
#include <TestViewOfClass.hpp>

#include <TestSharedAlloc.hpp>
#include <TestViewMapping.hpp>

#include <TestRange.hpp>
#include <TestTeam.hpp>
#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestAggregate.hpp>
#include <TestAggregateReduction.hpp>
#include <TestCompilerMacros.hpp>
#include <TestCXX11.hpp>
#include <TestCXX11Deduction.hpp>
#include <TestTeamVector.hpp>
#include <TestMemorySpaceTracking.hpp>
#include <TestTemplateMetaFunctions.hpp>

namespace Test {

class openmp : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
    const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
    const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

    const unsigned threads_count = std::max( 1u , numa_count ) *
                                   std::max( 2u , ( cores_per_numa * threads_per_core ) / 2 );

    Kokkos::OpenMP::initialize( threads_count );
    Kokkos::OpenMP::print_configuration( std::cout , true );
  }

  static void TearDownTestCase()
  {
    Kokkos::OpenMP::finalize();

    omp_set_num_threads(1);

    ASSERT_EQ( 1 , omp_get_max_threads() );
  }
};


TEST_F( openmp , impl_shared_alloc ) {
  test_shared_alloc< Kokkos::HostSpace , Kokkos::OpenMP >();
}

TEST_F( openmp , impl_view_mapping ) {
  test_view_mapping< Kokkos::OpenMP >();
  test_view_mapping_subview< Kokkos::OpenMP >();
  test_view_mapping_operator< Kokkos::OpenMP >();
  TestViewMappingAtomic< Kokkos::OpenMP >::run();
}

TEST_F( openmp, view_impl) {
  test_view_impl< Kokkos::OpenMP >();
}

TEST_F( openmp, view_api) {
  TestViewAPI< double , Kokkos::OpenMP >();
}

TEST_F( openmp , view_nested_view )
{
  ::Test::view_nested_view< Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_auto_1d_left ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutLeft,Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_auto_1d_right ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutRight,Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_auto_1d_stride ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutStride,Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_assign_strided ) {
  TestViewSubview::test_1d_strided_assignment< Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_left_0 ) {
  TestViewSubview::test_left_0< Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_left_1 ) {
  TestViewSubview::test_left_1< Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_left_2 ) {
  TestViewSubview::test_left_2< Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_left_3 ) {
  TestViewSubview::test_left_3< Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_right_0 ) {
  TestViewSubview::test_right_0< Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_right_1 ) {
  TestViewSubview::test_right_1< Kokkos::OpenMP >();
}

TEST_F( openmp, view_subview_right_3 ) {
  TestViewSubview::test_right_3< Kokkos::OpenMP >();
}



TEST_F( openmp , range_tag )
{
  TestRange< Kokkos::OpenMP >::test_for(1000);
  TestRange< Kokkos::OpenMP >::test_reduce(1000);
  TestRange< Kokkos::OpenMP >::test_scan(1000);
}

TEST_F( openmp , team_tag )
{
  TestTeamPolicy< Kokkos::OpenMP >::test_for(1000);
  TestTeamPolicy< Kokkos::OpenMP >::test_reduce(1000);
}

TEST_F( openmp, long_reduce) {
  TestReduce< long ,   Kokkos::OpenMP >( 1000000 );
}

TEST_F( openmp, double_reduce) {
  TestReduce< double ,   Kokkos::OpenMP >( 1000000 );
}

TEST_F( openmp, long_reduce_dynamic ) {
  TestReduceDynamic< long ,   Kokkos::OpenMP >( 1000000 );
}

TEST_F( openmp, double_reduce_dynamic ) {
  TestReduceDynamic< double ,   Kokkos::OpenMP >( 1000000 );
}

TEST_F( openmp, long_reduce_dynamic_view ) {
  TestReduceDynamicView< long ,   Kokkos::OpenMP >( 1000000 );
}

TEST_F( openmp, team_long_reduce) {
  TestReduceTeam< long ,   Kokkos::OpenMP >( 100000 );
}

TEST_F( openmp, team_double_reduce) {
  TestReduceTeam< double ,   Kokkos::OpenMP >( 100000 );
}

TEST_F( openmp, team_shared_request) {
  TestSharedTeam< Kokkos::OpenMP >();
}

#if defined(KOKKOS_HAVE_CXX11_DISPATCH_LAMBDA) && !defined(KOKKOS_HAVE_CUDA)
TEST_F( openmp, team_lambda_shared_request) {
  TestLambdaSharedTeam< Kokkos::OpenMP >();
}
#endif

TEST_F( openmp , atomics )
{
  const int loop_count = 1e4 ;

  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::OpenMP>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::OpenMP>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::OpenMP>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::OpenMP>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::OpenMP>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::OpenMP>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::OpenMP>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::OpenMP>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::OpenMP>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::OpenMP>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::OpenMP>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::OpenMP>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::OpenMP>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::OpenMP>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::OpenMP>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::OpenMP>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::OpenMP>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::OpenMP>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::OpenMP>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::OpenMP>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::OpenMP>(100,3) ) );

#if defined( KOKKOS_ENABLE_ASM )
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::OpenMP>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::OpenMP>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::OpenMP>(100,3) ) );
#endif
}

TEST_F( openmp , view_remap )
{
  enum { N0 = 3 , N1 = 2 , N2 = 8 , N3 = 9 };

  typedef Kokkos::View< double*[N1][N2][N3] ,
                             Kokkos::LayoutRight ,
                             Kokkos::OpenMP > output_type ;

  typedef Kokkos::View< int**[N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::OpenMP > input_type ;

  typedef Kokkos::View< int*[N0][N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::OpenMP > diff_type ;

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


TEST_F( openmp , view_aggregate )
{
  TestViewAggregate< Kokkos::OpenMP >();
  TestViewAggregateReduction< Kokkos::OpenMP >();
}

//----------------------------------------------------------------------------

TEST_F( openmp , scan )
{
  TestScan< Kokkos::OpenMP >::test_range( 1 , 1000 );
  TestScan< Kokkos::OpenMP >( 1000000 );
  TestScan< Kokkos::OpenMP >( 10000000 );
  Kokkos::OpenMP::fence();
}


TEST_F( openmp , team_scan )
{
  TestScanTeam< Kokkos::OpenMP >( 10000 );
  TestScanTeam< Kokkos::OpenMP >( 10000 );
}

//----------------------------------------------------------------------------

TEST_F( openmp , compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::OpenMP >() ) );
}

//----------------------------------------------------------------------------

TEST_F( openmp , memory_space )
{
  TestMemorySpace< Kokkos::OpenMP >();
}

//----------------------------------------------------------------------------

TEST_F( openmp , template_meta_functions )
{
  TestTemplateMetaFunctions<int, Kokkos::OpenMP >();
}

//----------------------------------------------------------------------------

#if defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP )
TEST_F( openmp , cxx11 )
{
  if ( Kokkos::Impl::is_same< Kokkos::DefaultExecutionSpace , Kokkos::OpenMP >::value ) {
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::OpenMP >(1) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::OpenMP >(2) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::OpenMP >(3) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::OpenMP >(4) ) );
  }
}
#endif

TEST_F( openmp , reduction_deduction )
{
  TestCXX11::test_reduction_deduction< Kokkos::OpenMP >();
}

TEST_F( openmp , team_vector )
{
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(0) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(1) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(2) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(3) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(4) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(5) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(6) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(7) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(8) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(9) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::OpenMP >(10) ) );
}
} // namespace test

