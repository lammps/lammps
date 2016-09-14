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

#include <iostream>

#include <Kokkos_Core.hpp>

//----------------------------------------------------------------------------

#include <Cuda/Kokkos_Cuda_TaskPolicy.hpp>
#include <impl/Kokkos_ViewTileLeft.hpp>
#include <TestTile.hpp>

//----------------------------------------------------------------------------

#include <TestSharedAlloc.hpp>
#include <TestViewMapping.hpp>

#include <TestViewImpl.hpp>
#include <TestAtomic.hpp>

#include <TestViewAPI.hpp>
#include <TestViewSubview.hpp>
#include <TestViewOfClass.hpp>

#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestRange.hpp>
#include <TestTeam.hpp>
#include <TestAggregate.hpp>
#include <TestAggregateReduction.hpp>
#include <TestCompilerMacros.hpp>
#include <TestMemorySpaceTracking.hpp>
#include <TestMemoryPool.hpp>
#include <TestTeamVector.hpp>
#include <TestTemplateMetaFunctions.hpp>
#include <TestCXX11Deduction.hpp>

#include <TestTaskPolicy.hpp>
#include <TestPolicyConstruction.hpp>

//----------------------------------------------------------------------------

class cuda : public ::testing::Test {
protected:
  static void SetUpTestCase();
  static void TearDownTestCase();
};

//----------------------------------------------------------------------------

namespace Test {

TEST_F( cuda, range_tag )
{
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_for(3);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_reduce(3);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_scan(3);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_for(3);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(3);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_scan(3);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_for(1000);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_reduce(1000);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_scan(1000);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_for(1001);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(1001);
  TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_scan(1001);
  //TestRange< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_dynamic_policy(1000);
}

TEST_F( cuda, team_tag )
{
  TestTeamPolicy< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_for(3);
  TestTeamPolicy< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_reduce(3);
  TestTeamPolicy< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_for(3);
  TestTeamPolicy< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(3);
  TestTeamPolicy< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_for(1000);
  TestTeamPolicy< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >::test_reduce(1000);
  TestTeamPolicy< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_for(1000);
  TestTeamPolicy< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(1000);
}

TEST_F( cuda, reduce )
{
  TestReduce< long ,   Kokkos::Cuda >( 10000000 );
  TestReduce< double , Kokkos::Cuda >( 1000000 );
  TestReduce< int , Kokkos::Cuda >( 0 );
}

TEST_F( cuda , reducers )
{
  TestReducers<int, Kokkos::Cuda>::execute_integer();
  TestReducers<size_t, Kokkos::Cuda>::execute_integer();
  TestReducers<double, Kokkos::Cuda>::execute_float();
  TestReducers<Kokkos::complex<double>, Kokkos::Cuda>::execute_basic();
}

TEST_F( cuda, reduce_team )
{
  TestReduceTeam< long ,   Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >( 3 );
  TestReduceTeam< long ,   Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >( 3 );
  TestReduceTeam< long ,   Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >( 100000 );
  TestReduceTeam< long ,   Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >( 100000 );
  TestReduceTeam< double ,   Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >( 3 );
  TestReduceTeam< double ,   Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >( 3 );
  TestReduceTeam< double ,   Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >( 100000 );
  TestReduceTeam< double ,   Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >( 100000 );
}

TEST_F( cuda, shared_team )
{
  TestSharedTeam< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >();
  TestSharedTeam< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >();
}

#if defined (KOKKOS_HAVE_CXX11_DISPATCH_LAMBDA)
TEST_F( cuda, lambda_shared_team )
{
  TestLambdaSharedTeam< Kokkos::CudaSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >();
  TestLambdaSharedTeam< Kokkos::CudaUVMSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >();
  TestLambdaSharedTeam< Kokkos::CudaHostPinnedSpace, Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static>  >();
  TestLambdaSharedTeam< Kokkos::CudaSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >();
  TestLambdaSharedTeam< Kokkos::CudaUVMSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >();
  TestLambdaSharedTeam< Kokkos::CudaHostPinnedSpace, Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic>  >();
}
#endif

TEST_F( cuda, shmem_size) {
  TestShmemSize< Kokkos::Cuda >();
}

TEST_F( cuda, multi_level_scratch) {
  TestMultiLevelScratchTeam< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Static> >();
  TestMultiLevelScratchTeam< Kokkos::Cuda , Kokkos::Schedule<Kokkos::Dynamic> >();
}

TEST_F( cuda, reduce_dynamic )
{
  TestReduceDynamic< long ,   Kokkos::Cuda >( 10000000 );
  TestReduceDynamic< double , Kokkos::Cuda >( 1000000 );
}

TEST_F( cuda, reduce_dynamic_view )
{
  TestReduceDynamicView< long ,   Kokkos::Cuda >( 10000000 );
  TestReduceDynamicView< double , Kokkos::Cuda >( 1000000 );
}

}
