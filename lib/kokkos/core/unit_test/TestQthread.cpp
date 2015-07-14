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
#include <Kokkos_Qthread.hpp>

#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp>

//----------------------------------------------------------------------------

#include <TestViewImpl.hpp>
#include <TestAtomic.hpp>

#include <TestViewAPI.hpp>

#include <TestTeam.hpp>
#include <TestRange.hpp>
#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestAggregate.hpp>
#include <TestCompilerMacros.hpp>
#include <TestTaskPolicy.hpp>
// #include <TestTeamVector.hpp>

namespace Test {

class qthread : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
    const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
    const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

    int threads_count = std::max( 1u , numa_count )
                      * std::max( 2u , ( cores_per_numa * threads_per_core ) / 2 );
    Kokkos::Qthread::initialize( threads_count );
    Kokkos::Qthread::print_configuration( std::cout , true );
  }

  static void TearDownTestCase()
  {
    Kokkos::Qthread::finalize();
  }
};

TEST_F( qthread , compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Qthread >() ) );
}

TEST_F( qthread, view_impl) {
  test_view_impl< Kokkos::Qthread >();
}

TEST_F( qthread, view_api) {
  TestViewAPI< double , Kokkos::Qthread >();
}

TEST_F( qthread , range_tag )
{
  TestRange< Kokkos::Qthread >::test_for(1000);
  TestRange< Kokkos::Qthread >::test_reduce(1000);
  TestRange< Kokkos::Qthread >::test_scan(1000);
}

TEST_F( qthread , team_tag )
{
  TestTeamPolicy< Kokkos::Qthread >::test_for( 1000 );
  TestTeamPolicy< Kokkos::Qthread >::test_reduce( 1000 );
}

TEST_F( qthread, long_reduce) {
  TestReduce< long ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, double_reduce) {
  TestReduce< double ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, long_reduce_dynamic ) {
  TestReduceDynamic< long ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, double_reduce_dynamic ) {
  TestReduceDynamic< double ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, long_reduce_dynamic_view ) {
  TestReduceDynamicView< long ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, team_long_reduce) {
  TestReduceTeam< long ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, team_double_reduce) {
  TestReduceTeam< double ,   Kokkos::Qthread >( 1000000 );
}


TEST_F( qthread , atomics )
{
  const int loop_count = 1e4 ;

  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Qthread>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Qthread>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Qthread>(100,3) ) );

#if defined( KOKKOS_ENABLE_ASM )
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Qthread>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Qthread>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<Kokkos::complex<double> ,Kokkos::Qthread>(100,3) ) );
#endif

}

TEST_F( qthread , view_remap )
{
  enum { N0 = 3 , N1 = 2 , N2 = 8 , N3 = 9 };

  typedef Kokkos::View< double*[N1][N2][N3] ,
                             Kokkos::LayoutRight ,
                             Kokkos::Qthread > output_type ;

  typedef Kokkos::View< int**[N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Qthread > input_type ;

  typedef Kokkos::View< int*[N0][N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Qthread > diff_type ;

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

TEST_F( qthread , view_aggregate )
{
  TestViewAggregate< Kokkos::Qthread >();
}

//----------------------------------------------------------------------------

TEST_F( qthread , scan )
{
  TestScan< Kokkos::Qthread >::test_range( 1 , 1000 );
  TestScan< Kokkos::Qthread >( 1000000 );
  TestScan< Kokkos::Qthread >( 10000000 );
  Kokkos::Qthread::fence();
}

TEST_F( qthread, team_shared ) {
  TestSharedTeam< Kokkos::Qthread >();
}

TEST_F( qthread , team_scan )
{
  TestScanTeam< Kokkos::Qthread >( 10 );
  TestScanTeam< Kokkos::Qthread >( 10000 );
}

#if defined (KOKKOS_HAVE_CXX11) && 0 /* disable */
TEST_F( qthread , team_vector )
{
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthread >(0) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthread >(1) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthread >(2) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthread >(3) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthread >(4) ) );
}
#endif

//----------------------------------------------------------------------------

TEST_F( qthread , task_policy )
{
  TestTaskPolicy::test_task_dep< Kokkos::Qthread >( 10 );
  for ( long i = 0 ; i < 25 ; ++i ) TestTaskPolicy::test_fib< Kokkos::Qthread >(i);
  for ( long i = 0 ; i < 35 ; ++i ) TestTaskPolicy::test_fib2< Kokkos::Qthread >(i);
}

#if defined( KOKKOS_HAVE_CXX11 )
TEST_F( qthread , task_team )
{
  std::cout << "qthread.task_team test disabled due to unresolved error causing the test to hang." << std::endl ;
  // TestTaskPolicy::test_task_team< Kokkos::Qthread >(1000);
}
#endif

//----------------------------------------------------------------------------

} // namespace test

