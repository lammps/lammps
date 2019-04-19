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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <cuda/TestCuda_Category.hpp>
#include <TestTeam.hpp>

namespace Test {

TEST_F( TEST_CATEGORY, team_for )
{
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_for( 0 );
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 0 );

  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_for( 2 );
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 2 );

  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_for( 1000 );
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_for( 1000 );
}


TEST_F( TEST_CATEGORY, team_reduce )
{
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 0 );
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 0 );
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 2 );
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 2 );
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_reduce( 1000 );
  TestTeamPolicy< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce( 1000 );
}

TEST_F( TEST_CATEGORY, team_broadcast )
{
  TestTeamBroadcast< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_teambroadcast( 0 );
  TestTeamBroadcast< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_teambroadcast( 0 );

  TestTeamBroadcast< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_teambroadcast( 2 );
  TestTeamBroadcast< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_teambroadcast( 2 );

  TestTeamBroadcast< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_teambroadcast( 16 );
  TestTeamBroadcast< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_teambroadcast( 16 );

  TestTeamBroadcast< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_teambroadcast( 1000 );
  TestTeamBroadcast< TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_teambroadcast( 1000 );
}

}

#include <TestTeamVector.hpp>


