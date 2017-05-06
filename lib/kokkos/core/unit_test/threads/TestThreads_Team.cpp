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

TEST_F( threads, team_tag )
{
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_for( 0 );
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 0 );
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 0 );
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 0 );

  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_for( 2 );
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 2 );
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 2 );
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 2 );

  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_for( 1000 );
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 1000 );
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 1000 );
  TestTeamPolicy< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 1000 );
}

TEST_F( threads, team_shared_request )
{
  TestSharedTeam< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >();
  TestSharedTeam< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >();
}

TEST_F( threads, team_scratch_request )
{
  TestScratchTeam< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >();
  TestScratchTeam< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >();
}

#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
TEST_F( threads, team_lambda_shared_request )
{
  TestLambdaSharedTeam< Kokkos::HostSpace, Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >();
  TestLambdaSharedTeam< Kokkos::HostSpace, Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >();
}
#endif

TEST_F( threads, shmem_size )
{
  TestShmemSize< Kokkos::Threads >();
}

TEST_F( threads, multi_level_scratch )
{
  TestMultiLevelScratchTeam< Kokkos::Threads, Kokkos::Schedule<Kokkos::Static> >();
  TestMultiLevelScratchTeam< Kokkos::Threads, Kokkos::Schedule<Kokkos::Dynamic> >();
}

TEST_F( threads, team_vector )
{
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 0 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 1 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 2 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 3 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 4 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 5 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 6 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 7 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 8 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 9 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Threads >( 10 ) ) );
}

#ifdef KOKKOS_COMPILER_GNU
#if ( KOKKOS_COMPILER_GNU == 472 )
#define SKIP_TEST
#endif
#endif

#ifndef SKIP_TEST
TEST_F( threads, triple_nested_parallelism )
{
  TestTripleNestedReduce< double, Kokkos::Threads >( 8192, 2048, 32, 32 );
  TestTripleNestedReduce< double, Kokkos::Threads >( 8192, 2048, 32, 16 );
  TestTripleNestedReduce< double, Kokkos::Threads >( 8192, 2048, 16, 16 );
}
#endif

} // namespace Test
