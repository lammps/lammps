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

#include <qthreads/TestQthreads.hpp>

namespace Test {

TEST_F( qthreads, team_tag )
{
#if 0
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_for( 0 );
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 0 );
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 0 );
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 0 );

  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_for( 2 );
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 2 );
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 2 );
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 2 );

  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_for( 1000 );
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 1000 );
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 1000 );
  TestTeamPolicy< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 1000 );
#endif
}

TEST_F( qthreads, team_shared_request )
{
#if 0
  TestSharedTeam< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >();
  TestSharedTeam< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >();
#endif
}

TEST_F( qthreads, team_scratch_request )
{
#if 0
  TestScratchTeam< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >();
  TestScratchTeam< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >();
#endif
}

#if defined( KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA )
TEST_F( qthreads, team_lambda_shared_request )
{
#if 0
  TestLambdaSharedTeam< Kokkos::HostSpace, Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >();
  TestLambdaSharedTeam< Kokkos::HostSpace, Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >();
#endif
}
#endif

TEST_F( qthreads, shmem_size )
{
#if 0
  TestShmemSize< Kokkos::Qthreads >();
#endif
}

TEST_F( qthreads, multi_level_scratch )
{
#if 0
  TestMultiLevelScratchTeam< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Static> >();
  TestMultiLevelScratchTeam< Kokkos::Qthreads, Kokkos::Schedule<Kokkos::Dynamic> >();
#endif
}

TEST_F( qthreads, team_vector )
{
#if 0
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 0 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 1 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 2 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 3 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 4 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 5 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 6 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 7 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 8 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 9 ) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Qthreads >( 10 ) ) );
#endif
}

#ifdef KOKKOS_COMPILER_GNU
#if ( KOKKOS_COMPILER_GNU == 472 )
#define SKIP_TEST
#endif
#endif

#ifndef SKIP_TEST
TEST_F( qthreads, triple_nested_parallelism )
{
#if 0
  TestTripleNestedReduce< double, Kokkos::Qthreads >( 8192, 2048, 32, 32 );
  TestTripleNestedReduce< double, Kokkos::Qthreads >( 8192, 2048, 32, 16 );
  TestTripleNestedReduce< double, Kokkos::Qthreads >( 8192, 2048, 16, 16 );
#endif
}
#endif

} // namespace Test
