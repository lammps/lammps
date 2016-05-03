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
#ifdef KOKKOS_LAMBDA
#undef KOKKOS_LAMBDA
#endif
#define KOKKOS_LAMBDA [=]

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
#include <TestMemoryPool.hpp>


#include <TestCXX11.hpp>
#include <TestCXX11Deduction.hpp>
#include <TestTeamVector.hpp>
#include <TestMemorySpaceTracking.hpp>
#include <TestTemplateMetaFunctions.hpp>

#include <TestPolicyConstruction.hpp>


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

TEST_F( openmp , range_tag )
{
  TestRange< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >::test_for(1000);
  TestRange< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >::test_reduce(1000);
  TestRange< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >::test_scan(1000);
  TestRange< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >::test_for(1001);
  TestRange< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(1001);
  TestRange< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >::test_scan(1001);
  TestRange< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >::test_dynamic_policy(1000);
}

TEST_F( openmp , team_tag )
{
  TestTeamPolicy< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >::test_for(1000);
  TestTeamPolicy< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >::test_reduce(1000);
  TestTeamPolicy< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >::test_for(1000);
  TestTeamPolicy< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(1000);
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
  TestReduceTeam< long ,   Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >( 3 );
  TestReduceTeam< long ,   Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >( 3 );
  TestReduceTeam< long ,   Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >( 100000 );
  TestReduceTeam< long ,   Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >( 100000 );
}

TEST_F( openmp, team_double_reduce) {
  TestReduceTeam< double ,   Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >( 3 );
  TestReduceTeam< double ,   Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >( 3 );
  TestReduceTeam< double ,   Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >( 100000 );
  TestReduceTeam< double ,   Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >( 100000 );
}

TEST_F( openmp, team_shared_request) {
  TestSharedTeam< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >();
  TestSharedTeam< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >();
}

TEST_F( openmp, team_scratch_request) {
  TestScratchTeam< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >();
  TestScratchTeam< Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >();
}

#if defined(KOKKOS_HAVE_CXX11_DISPATCH_LAMBDA) 
TEST_F( openmp, team_lambda_shared_request) {
  TestLambdaSharedTeam< Kokkos::HostSpace, Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Static> >();
  TestLambdaSharedTeam< Kokkos::HostSpace, Kokkos::OpenMP , Kokkos::Schedule<Kokkos::Dynamic> >();
}
#endif

} // namespace test

