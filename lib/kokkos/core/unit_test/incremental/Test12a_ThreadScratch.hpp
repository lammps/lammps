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
struct ThreadScratch {
  using policy_t = Kokkos::TeamPolicy<ExecSpace>;
  using team_t   = typename Kokkos::TeamPolicy<ExecSpace>::member_type;
  using data_t   = Kokkos::View<size_t **, ExecSpace>;

  using scratch_t = Kokkos::View<size_t *, ExecSpace,
                                 Kokkos::MemoryTraits<Kokkos::Unmanaged> >;

  int sX, sY;
  data_t v;

  const int scratch_level = 1;
  KOKKOS_FUNCTION
  void operator()(const team_t &team) const {
    // Allocate and use scratch pad memory
    scratch_t v_S(team.thread_scratch(scratch_level), sY);
    int n = team.league_rank();

    for (int i = 0; i < sY; ++i) v_S(i) = 0;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, sX), [&](const int m) {
    // FIXME_SYCL This deadlocks in the subgroup_barrier when running on CUDA
    // devices.
#ifdef KOKKOS_ENABLE_SYCL
      for (int k = 0; k < sY; ++k) v_S(k) += sX * sY * n + sY * m + k;
#else
      Kokkos::parallel_for(
          Kokkos::ThreadVectorRange(team, sY),
          [&](const int k) { v_S(k) += sX * sY * n + sY * m + k; });
#endif
    });

    team.team_barrier();

    for (int i = 0; i < sY; ++i) {
      v(n, team.team_rank()) += v_S(i);
    }
  }

  void run(const int pN, const int sX_, const int sY_) {
    sX = sX_;
    sY = sY_;

    int scratchSize = scratch_t::shmem_size(sY);
    // So this works with deprecated code enabled:
    policy_t policy =
        policy_t(pN, Kokkos::AUTO, 1)
            .set_scratch_size(scratch_level, Kokkos::PerThread(scratchSize));

    int max_team_size = policy.team_size_max(*this, Kokkos::ParallelForTag());
    v                 = data_t("Matrix", pN, max_team_size);

    Kokkos::parallel_for(
        "Test12a_ThreadScratch",
        policy_t(pN, max_team_size)
            .set_scratch_size(scratch_level, Kokkos::PerThread(scratchSize)),
        *this);

    Kokkos::fence();
    auto v_H = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);

    size_t check   = 0;
    const size_t s = pN * sX * sY;
    for (int n = 0; n < pN; ++n)
      for (int m = 0; m < max_team_size; ++m) {
        check += v_H(n, m);
      }
    ASSERT_EQ(s * (s - 1) / 2, check);
  }
};

TEST(TEST_CATEGORY, IncrTest_12a_ThreadScratch) {
  ThreadScratch<TEST_EXECSPACE> test;
#ifdef KOKKOS_ENABLE_OPENACC  // FIXME_OPENACC
  GTEST_SKIP() << "skipping since scratch memory is not yet implemented in the "
                  "OpenACC backend";
#endif
  // FIXME_OPENMPTARGET - team_size has to be a multiple of 32 for the tests to
  // pass in the Release and RelWithDebInfo builds. Does not need the team_size
  // to be a multiple of 32 for the Debug builds.
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  test.run(1, 32, 9);
  test.run(2, 64, 22);
  test.run(14, 128, 321);
#else
  test.run(1, 55, 9);
  test.run(2, 4, 22);
  test.run(14, 277, 321);
#endif
}

}  // namespace Test
