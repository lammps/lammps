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

// @Kokkos_Feature_Level_Required:13
// Unit test for hierarchical parallelism
// Create concurrent work hierarchically and verify if
// sum of created processing units corresponds to expected value

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

using SCALAR_TYPE = int;

namespace Test {

template <class ExecSpace>
struct Hierarchical_Red_C {
  void run(const int pN, const int sX, const int sY) {
    using team_policy = Kokkos::TeamPolicy<ExecSpace>;
    using member_type = typename Kokkos::TeamPolicy<ExecSpace>::member_type;

    using viewDataType = Kokkos::View<SCALAR_TYPE *, ExecSpace>;
    viewDataType v("Vector", pN);

    Kokkos::parallel_for(
        "Team", team_policy(pN, Kokkos::AUTO),
        KOKKOS_LAMBDA(const member_type &team) {
          int n           = team.league_rank();
          SCALAR_TYPE out = 0;

          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(team, sX),
              [=](const int i, SCALAR_TYPE &tmp) {
                SCALAR_TYPE out_inner = 0;
                Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange(team, sY),
                    [=](const int k, int &tmp_inner) {
                      tmp_inner += n * sX * v.extent(0) + sX * i + k;
                    },
                    out_inner);

                Kokkos::single(Kokkos::PerThread(team),
                               [&]() { tmp += out_inner; });
              },
              out);

          Kokkos::single(Kokkos::PerTeam(team), [&]() { v(n) += out; });
        });

    Kokkos::fence();
    auto v_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);

    SCALAR_TYPE check = 0;
    SCALAR_TYPE ref   = 0;
    for (int i = 0; i < pN; ++i) {
      check += v_H(i);
      for (int j = 0; j < sX; ++j)
        for (int k = 0; k < sY; ++k) ref += i * sX * pN + sX * j + k;
    }
    ASSERT_EQ(check, ref);
  }
};

TEST(TEST_CATEGORY, IncrTest_13c_Hierarchical_Red) {
  Hierarchical_Red_C<TEST_EXECSPACE> test;
  test.run(1, 4, 8);
  test.run(2, 39, 12);
  test.run(39, 3, 235);
}

}  // namespace Test
