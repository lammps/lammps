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

    const int scratch_level = 1;

    Kokkos::parallel_for(
        "Team",
        policy_t(pN, Kokkos::AUTO, 1)
            .set_scratch_size(scratch_level, Kokkos::PerTeam(scratchSize)),
        KOKKOS_LAMBDA(const team_t &team) {
          // Allocate and use scratch pad memory
          scratch_t v_S(team.team_scratch(scratch_level), sX, sY);
          int n = team.league_rank();

          Kokkos::parallel_for(
              Kokkos::TeamThreadRange(team, sX), [&](const int m) {
      // FIXME_SYCL This deadlocks in the subgroup_barrier
      // when running on CUDA devices.
#ifdef KOKKOS_ENABLE_SYCL
                for (int k = 0; k < sY; ++k) {
                  v_S(m, k) =
                      v_S.extent(0) * v_S.extent(1) * n + v_S.extent(1) * m + k;
                }
#else
                Kokkos::parallel_for(
                    Kokkos::ThreadVectorRange(team, sY), [&](const int k) {
                      v_S(m, k) = v_S.extent(0) * v_S.extent(1) * n +
                                  v_S.extent(1) * m + k;
                    });
#endif
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
#ifdef KOKKOS_ENABLE_OPENACC  // FIXME_OPENACC
  GTEST_SKIP() << "skipping since scratch memory is not yet implemented in the "
                  "OpenACC backend";
#endif
  // FIXME_OPENMPTARGET - team_size has to be a multiple of 32 for the tests to
  // pass in the Release and RelWithDebInfo builds. Does not need the team_size
  // to be a multiple of 32 for the Debug builds.
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  test.run(1, 32, 4);
  test.run(4, 64, 10);
  test.run(14, 128, 20);
#else
  test.run(1, 4, 4);
  test.run(4, 7, 10);
  test.run(14, 277, 321);
#endif
}

}  // namespace Test
