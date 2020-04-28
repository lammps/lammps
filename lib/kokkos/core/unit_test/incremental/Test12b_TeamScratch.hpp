/*
//@HEADER
// ************************************************************************
//
// Kokkos v. 3.0
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

// @Kokkos_Feature_Level_Required:12
// Unit test for hierarchical parallelism
// Create concurrent work hierarchically and verify if
// contributions of paticipating processing units corresponds to expected value
// Use a scratch pad memory for each team
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace Test {

template <class ExecSpace>
struct TeamScratch {
  void run(const int pN, const int sX, const int sY) {
    using policy_t = Kokkos::TeamPolicy<ExecSpace>;
    using team_t   = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
    using data_t   = Kokkos::View<size_t **, ExecSpace>;
    data_t v("Matrix", pN, sX);

    using scratch_t = Kokkos::View<size_t **, ExecSpace,
                                   Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
    int scratchSize = scratch_t::shmem_size(sX, sY);

    Kokkos::parallel_for(
        "Team",
        policy_t(pN, Kokkos::AUTO)
            .set_scratch_size(1, Kokkos::PerTeam(scratchSize)),
        KOKKOS_LAMBDA(const team_t &team) {
          // Allocate and use scratch pad memory
          scratch_t v_S(team.team_scratch(1), sX, sY);
          int n = team.league_rank();

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, sX), [&](const int m) {
                Kokkos::parallel_for(
                    Kokkos::ThreadVectorRange(team, sY), [&](const int k) {
                      v_S(m, k) = v_S.extent(0) * v_S.extent(1) * n +
                                  v_S.extent(1) * m + k;
                    });
              });

          team.team_barrier();

          // Sum up contributions and reduce by one dimension
          Kokkos::parallel_for(Kokkos::TeamThreadRange(team, sX),
                               [&](const int m) {
                                 for (int i = 0; i < sY; ++i)
                                   v(n, m) += v_S(m, i);
                               });
        });

    Kokkos::fence();
    auto v_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);

    size_t check   = 0;
    const size_t s = pN * sX * sY;
    for (int n = 0; n < pN; ++n)
      for (int m = 0; m < sX; ++m) check += v_H(n, m);
    ASSERT_EQ(check, s * (s - 1) / 2);
  }
};

TEST(TEST_CATEGORY, IncrTest_12b_TeamScratch) {
  TeamScratch<TEST_EXECSPACE> test;
  test.run(1, 4, 4);
  test.run(4, 7, 10);
  test.run(14, 277, 321);
}

}  // namespace Test
