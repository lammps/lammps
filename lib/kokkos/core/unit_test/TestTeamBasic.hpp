/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef KOKKOS_TEST_TEAM_BASIC_HPP
#define KOKKOS_TEST_TEAM_BASIC_HPP
#include <TestTeam.hpp>

namespace Test {

TEST(TEST_CATEGORY, team_for) {
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_for(
      0);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_for(
      0);

  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_for(
      2);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_for(
      2);

  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static> >::test_for(
      1000);
  TestTeamPolicy<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic> >::test_for(
      1000);
}

// FIXME_OPENMPTARGET wrong results
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, team_reduce) {
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Static> >::test_reduce(0);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(0);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Static> >::test_reduce(2);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(2);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Static> >::test_reduce(1000);
  TestTeamPolicy<TEST_EXECSPACE,
                 Kokkos::Schedule<Kokkos::Dynamic> >::test_reduce(1000);
}
#endif

TEST(TEST_CATEGORY, team_broadcast_long) {
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long>::test_teambroadcast(0, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long>::test_teambroadcast(0, 1);

  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long>::test_teambroadcast(2, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long>::test_teambroadcast(2, 1);

  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long>::test_teambroadcast(16, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long>::test_teambroadcast(16, 1);

  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                    long>::test_teambroadcast(1000, 1);
  TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                    long>::test_teambroadcast(1000, 1);
}

TEST(TEST_CATEGORY, team_broadcast_char) {
  {
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      unsigned char>::test_teambroadcast(0, 1);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      unsigned char>::test_teambroadcast(0, 1);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      unsigned char>::test_teambroadcast(2, 1);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      unsigned char>::test_teambroadcast(2, 1);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      unsigned char>::test_teambroadcast(16, 1);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      unsigned char>::test_teambroadcast(16, 1);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      long>::test_teambroadcast(1000, 1);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      long>::test_teambroadcast(1000, 1);
  }
}

TEST(TEST_CATEGORY, team_broadcast_float) {
  {
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      float>::test_teambroadcast(0, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      float>::test_teambroadcast(0, 1.3);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      float>::test_teambroadcast(2, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      float>::test_teambroadcast(2, 1.3);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      float>::test_teambroadcast(16, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      float>::test_teambroadcast(16, 1.3);

    // FIXME_CUDA
#ifdef KOKKOS_ENABLE_CUDA
    if (!std::is_same<TEST_EXECSPACE, Kokkos::Cuda>::value)
#endif
    // FIXME_HIP
#ifdef KOKKOS_ENABLE_HIP
      if (!std::is_same<TEST_EXECSPACE, Kokkos::Experimental::HIP>::value)
#endif
      {
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                          float>::test_teambroadcast(1000, 1.3);
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                          float>::test_teambroadcast(1000, 1.3);
      }
  }
}

TEST(TEST_CATEGORY, team_broadcast_double) {
  {
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      double>::test_teambroadcast(0, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      double>::test_teambroadcast(0, 1.3);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      double>::test_teambroadcast(2, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      double>::test_teambroadcast(2, 1.3);

    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                      double>::test_teambroadcast(16, 1.3);
    TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,
                      double>::test_teambroadcast(16, 1.3);

    // FIXME_CUDA
#ifdef KOKKOS_ENABLE_CUDA
    if (!std::is_same<TEST_EXECSPACE, Kokkos::Cuda>::value)
#endif
    // FIXME_HIP
#ifdef KOKKOS_ENABLE_HIP
      if (!std::is_same<TEST_EXECSPACE, Kokkos::Experimental::HIP>::value)
#endif
      {
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Static>,
                          double>::test_teambroadcast(1000, 1.3);
        TestTeamBroadcast<TEST_EXECSPACE, Kokkos::Schedule<Kokkos::Dynamic>,

                          double>::test_teambroadcast(1000, 1.3);
      }
  }
}

TEST(TEST_CATEGORY, team_handle_by_value) {
  { TestTeamPolicyHandleByValue<TEST_EXECSPACE>(); }
}

}  // namespace Test

#ifndef KOKKOS_ENABLE_OPENMPTARGET
#include <TestTeamVector.hpp>
#endif
#endif
