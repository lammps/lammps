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

#include <Test_InterOp_Streams.hpp>

namespace {

void test_policies(TEST_EXECSPACE exec0, Kokkos::View<int *, TEST_EXECSPACE> v0,
                   TEST_EXECSPACE exec, Kokkos::View<int *, TEST_EXECSPACE> v) {
  using MemorySpace = typename TEST_EXECSPACE::memory_space;

  exec.fence();
  exec0.fence();

  Kokkos::deep_copy(exec, v, 5);
  Kokkos::deep_copy(exec0, v0, 5);

  Kokkos::deep_copy(v, v0);

  int sum;
  int sum0;

  Kokkos::parallel_for("Test::Range_0",
                       Kokkos::RangePolicy<TEST_EXECSPACE>(exec0, 0, 100),
                       Test::FunctorRange<MemorySpace>(v0));
  Kokkos::parallel_for("Test::Range",
                       Kokkos::RangePolicy<TEST_EXECSPACE>(exec, 0, 100),
                       Test::FunctorRange<MemorySpace>(v));
  exec.fence();
  exec0.fence();
  Kokkos::parallel_reduce(
      "Test::RangeReduce_0",
      Kokkos::RangePolicy<TEST_EXECSPACE, Kokkos::LaunchBounds<128, 2>>(exec0,
                                                                        0, 100),
      Test::FunctorRangeReduce<MemorySpace>(v0), sum0);
  Kokkos::parallel_reduce(
      "Test::RangeReduce",
      Kokkos::RangePolicy<TEST_EXECSPACE, Kokkos::LaunchBounds<128, 2>>(exec, 0,
                                                                        100),
      Test::FunctorRangeReduce<MemorySpace>(v), sum);
  ASSERT_EQ(600, sum0);
  ASSERT_EQ(600, sum);

  Kokkos::parallel_for("Test::MDRange_0",
                       Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                           exec0, {0, 0}, {10, 10}),
                       Test::FunctorMDRange<MemorySpace>(v0));
  Kokkos::parallel_for("Test::MDRange",
                       Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                           exec, {0, 0}, {10, 10}),
                       Test::FunctorMDRange<MemorySpace>(v));
  Kokkos::parallel_reduce("Test::MDRangeReduce_0",
                          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>,
                                                Kokkos::LaunchBounds<128, 2>>(
                              exec0, {0, 0}, {10, 10}),
                          Test::FunctorMDRangeReduce<MemorySpace>(v0), sum0);
  Kokkos::parallel_reduce("Test::MDRangeReduce",
                          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>,
                                                Kokkos::LaunchBounds<128, 2>>(
                              exec, {0, 0}, {10, 10}),
                          Test::FunctorMDRangeReduce<MemorySpace>(v), sum);
  ASSERT_EQ(700, sum0);
  ASSERT_EQ(700, sum);

  Kokkos::parallel_for("Test::Team_0",
                       Kokkos::TeamPolicy<TEST_EXECSPACE>(exec0, 10, 10),
                       Test::FunctorTeam<MemorySpace, TEST_EXECSPACE>(v0));
  Kokkos::parallel_for("Test::Team",
                       Kokkos::TeamPolicy<TEST_EXECSPACE>(exec, 10, 10),
                       Test::FunctorTeam<MemorySpace, TEST_EXECSPACE>(v));
  Kokkos::parallel_reduce(
      "Test::Team_0",
      Kokkos::TeamPolicy<TEST_EXECSPACE, Kokkos::LaunchBounds<128, 2>>(exec0,
                                                                       10, 10),
      Test::FunctorTeamReduce<MemorySpace, TEST_EXECSPACE>(v0), sum0);
  Kokkos::parallel_reduce(
      "Test::Team",
      Kokkos::TeamPolicy<TEST_EXECSPACE, Kokkos::LaunchBounds<128, 2>>(exec, 10,
                                                                       10),
      Test::FunctorTeamReduce<MemorySpace, TEST_EXECSPACE>(v), sum);
  ASSERT_EQ(800, sum0);
  ASSERT_EQ(800, sum);
}

struct ScratchFunctor {
  int scratch_size;
  int R;

  ScratchFunctor(int scratch_size_, int R_)
      : scratch_size(scratch_size_), R(R_) {}

  KOKKOS_FUNCTION
  void operator()(const Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type &team,
                  int &error_accum) const {
    Kokkos::View<int *, TEST_EXECSPACE::scratch_memory_space> scratch_mem(
        team.team_scratch(1), scratch_size);

    // Initialize scratch memory
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, scratch_size),
                         [&](int i) { scratch_mem(i) = 0; });
    team.team_barrier();

    // Increment each entry in scratch memory R times
    for (int r = 0; r < R; ++r) {
      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, 0, scratch_size),
                           [&](int i) { scratch_mem(i) += 1; });
    }
    team.team_barrier();

    // Check that each scratch entry has been incremented exactly R times
    int team_error_accum;
    auto R_loc = R;  // avoid implicit capture of this
    Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, 0, scratch_size),
        [&](int i, int &tsum) {
          if (scratch_mem(i) != R_loc) {
            tsum += 1;
          }
        },
        team_error_accum);
    Kokkos::single(Kokkos::PerTeam(team),
                   [&]() { error_accum += team_error_accum; });
  }
};

void test_scratch(TEST_EXECSPACE exec0, TEST_EXECSPACE exec1) {
  constexpr int N            = 10;
  constexpr int R            = 1000;
  constexpr int scratch_size = 100;
  using ScratchType = Kokkos::View<int *, TEST_EXECSPACE::scratch_memory_space>;

  // Test allocating and using scratch space
  ScratchFunctor f(scratch_size, R);

  auto policy0 =
      Kokkos::TeamPolicy<TEST_EXECSPACE>(exec0, N, 10)
          .set_scratch_size(
              1, Kokkos::PerTeam(ScratchType::shmem_size(scratch_size)));
  auto policy1 =
      Kokkos::TeamPolicy<TEST_EXECSPACE>(exec1, N, 10)
          .set_scratch_size(
              1, Kokkos::PerTeam(ScratchType::shmem_size(scratch_size)));

  int error0, error1;

  Kokkos::parallel_reduce("test_scratch_device_0", policy0, f, error0);
  Kokkos::parallel_reduce("test_scratch_device_1", policy1, f, error1);
  ASSERT_EQ(error0, 0);
  ASSERT_EQ(error1, 0);

  // Request larger scratch size to trigger a realloc and test
  const auto new_scratch_size = scratch_size + 10;
  ScratchFunctor f_more_scratch(new_scratch_size, R);

  auto policy0_more_scratch =
      Kokkos::TeamPolicy<TEST_EXECSPACE>(exec0, N, 10)
          .set_scratch_size(
              1, Kokkos::PerTeam(ScratchType::shmem_size(new_scratch_size)));
  auto policy1_more_scratch =
      Kokkos::TeamPolicy<TEST_EXECSPACE>(exec1, N, 10)
          .set_scratch_size(
              1, Kokkos::PerTeam(ScratchType::shmem_size(new_scratch_size)));

  Kokkos::parallel_reduce("test_realloc_scratch_device_0", policy0_more_scratch,
                          f_more_scratch, error0);
  Kokkos::parallel_reduce("test_realloc_scratch_device_1", policy1_more_scratch,
                          f_more_scratch, error1);
  ASSERT_EQ(error0, 0);
  ASSERT_EQ(error1, 0);
}
}  // namespace
