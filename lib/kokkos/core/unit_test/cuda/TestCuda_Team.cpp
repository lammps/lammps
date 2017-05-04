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

#include <cuda/TestCuda.hpp>

namespace Test {

TEST_F( cuda, team_tag )
{
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >::test_for( 0 );
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 0 );
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 0 );
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 0 );

  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >::test_for( 2 );
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 2 );
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 2 );
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 2 );

  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >::test_for( 1000 );
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 1000 );
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 1000 );
  TestTeamPolicy< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 1000 );
}

TEST_F( cuda, team_shared_request )
{
  TestSharedTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >();
  TestSharedTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >();
}

// This tests request to much L0 scratch.
//TEST_F( cuda, team_scratch_request )
//{
//  TestScratchTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >();
//  TestScratchTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >();
//}

#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
TEST_F( cuda, team_lambda_shared_request )
{
  TestLambdaSharedTeam< Kokkos::CudaSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >();
  TestLambdaSharedTeam< Kokkos::CudaUVMSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >();
  TestLambdaSharedTeam< Kokkos::CudaHostPinnedSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >();
  TestLambdaSharedTeam< Kokkos::CudaSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >();
  TestLambdaSharedTeam< Kokkos::CudaUVMSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >();
  TestLambdaSharedTeam< Kokkos::CudaHostPinnedSpace, Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >();
}
#endif

TEST_F( cuda, shmem_size )
{
  TestShmemSize< Kokkos::Cuda >();
}

TEST_F( cuda, multi_level_scratch )
{
  TestMultiLevelScratchTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Static> >();
  TestMultiLevelScratchTeam< Kokkos::Cuda, Kokkos::Schedule<Kokkos::Dynamic> >();
}

#if !defined(KOKKOS_CUDA_CLANG_WORKAROUND) && !defined(KOKKOS_ARCH_PASCAL)
TEST_F( cuda, team_vector )
{
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 0 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 1 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 2 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 3 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 4 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 5 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 6 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 7 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 8 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 9 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Cuda >( 10 ) ) );
}
#endif

TEST_F( cuda, triple_nested_parallelism )
{
  TestTripleNestedReduce< double, Kokkos::Cuda >( 8192, 2048, 32, 32 );
  TestTripleNestedReduce< double, Kokkos::Cuda >( 8192, 2048, 32, 16 );
  TestTripleNestedReduce< double, Kokkos::Cuda >( 8192, 2048, 16, 16 );
}

} // namespace Test
