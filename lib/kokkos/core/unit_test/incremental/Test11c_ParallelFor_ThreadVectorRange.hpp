//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

// @Kokkos_Feature_Level_Required:11
// Unit test for hierarchical parallelism
// Create concurrent work hierarchically and verify if
// contributions of paticipating processing units corresponds to expected value

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace Test {

template <class ExecSpace>
struct Hierarchical_ForLoop_C {
  void run(const int pN, const int sX, const int sY, const int sZ) {
    using team_policy = Kokkos::TeamPolicy<ExecSpace>;
    using member_type = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

    using viewDataType = Kokkos::View<size_t ***, ExecSpace>;
    viewDataType v("Matrix", sX, sY, sZ);

    Kokkos::parallel_for(
        "Team", team_policy(pN, Kokkos::AUTO),
        KOKKOS_LAMBDA(const member_type &team) {
          int n  = team.league_rank();
          int ls = team.league_size();

          int startDim1 = n * (int)(sX / ls);
          int modDim1   = n == ls - 1 ? sX % ls : 0;

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, v.extent(1)), [&](const int m) {
                Kokkos::parallel_for(
                    Kokkos::ThreadVectorRange(team, v.extent(2)),
                    [&](const int k) {
                      for (int i = startDim1;
                           i < (startDim1 + (int)(sX / ls) + modDim1); ++i)
                        v(i, m, k) =
                            i * v.extent(1) * v.extent(2) + m * v.extent(2) + k;
                    });
              });
        });

    Kokkos::fence();
    auto v_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);

    size_t check   = 0;
    const size_t s = sX * sY * sZ;
    for (int i = 0; i < sX; ++i)
      for (int j = 0; j < sY; ++j)
        for (int k = 0; k < sZ; ++k) check += v_H(i, j, k);
    ASSERT_EQ(check, s * (s - 1) / 2);
  }
};

TEST(TEST_CATEGORY, IncrTest_11c_Hierarchical_ForLoop) {
  Hierarchical_ForLoop_C<TEST_EXECSPACE> test;
  test.run(4, 16, 16, 16);
  test.run(8, 12, 333, 16);
  test.run(12, 277, 321, 345);
}

}  // namespace Test
