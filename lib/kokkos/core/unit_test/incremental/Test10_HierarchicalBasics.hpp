// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// @Kokkos_Feature_Level_Required:10
// Unit test for hierarchical parallelism
// Create concurrent work hierarchically and verify if
// contributions of paticipating processing units corresponds to expected value

#include <gtest/gtest.h>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

namespace Test {

template <class ExecSpace>
struct HierarchicalBasics {
  using policy_t = Kokkos::TeamPolicy<ExecSpace>;
  using team_t   = typename policy_t::member_type;

  void run(const int nP, int nT) {
    int const concurrency = ExecSpace().concurrency();
    if (nT > concurrency) nT = concurrency;

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

    int check = 0;
    int ref   = nP * nT;
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
