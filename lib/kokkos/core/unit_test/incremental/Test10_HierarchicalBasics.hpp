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

// @Kokkos_Feature_Level_Required:10
// Unit test for hierarchical parallelism
// Create concurrent work hierarchically and verify if
// contributions of paticipating processing units corresponds to expected value

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace Test {

template <class ExecSpace>
struct HierarchicalBasics {
  using policy_t = Kokkos::TeamPolicy<ExecSpace>;
  using team_t   = typename policy_t::member_type;

  void run(const int nP, int nT) {
    if (nT > ExecSpace::concurrency()) nT = ExecSpace::concurrency();

    policy_t pol(nP, nT);

    ASSERT_EQ(pol.league_size(), nP);
    ASSERT_LE(pol.team_size(), nT);

    nT = pol.team_size();

    Kokkos::View<int **, ExecSpace> v("Array_A", nP, nT);
    Kokkos::parallel_for(
        "Teams", pol, KOKKOS_LAMBDA(const team_t &team) {
          const int tR = team.team_rank();
          const int tS = team.team_size();
          const int lR = team.league_rank();
          const int lS = team.league_size();
          if (lR < lS) {
            v(lR, tR) = lR * tS + tR;
          } else {
            v(lR, tR) = 100000;
          }
        });
    Kokkos::fence();
    auto h_v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);

    size_t check = 0;
    size_t ref   = nP * nT;
    for (int i = 0; i < nP; ++i)
      for (int j = 0; j < nT; ++j) check += h_v(i, j);

    ASSERT_EQ(check, ref * (ref - 1) / 2);
  }
};

TEST(TEST_CATEGORY, IncrTest_10_Hierarchical_Basics) {
  HierarchicalBasics<TEST_EXECSPACE> test;

  // OpenMPTarget backend only accepts >= 32 threads per team
#if defined(KOKKOS_ENABLE_OPENMPTARGET)
  test.run(1, 32);
  test.run(8, 64);
  test.run(11, 128);
#else
  test.run(1, 4);
  test.run(8, 16);
  test.run(11, 13);
#endif
}

}  // namespace Test
