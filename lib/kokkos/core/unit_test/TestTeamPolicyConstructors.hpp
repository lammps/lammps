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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace {

struct SomeTag {};

struct FunctorFor {
  KOKKOS_FUNCTION
  void operator()(
      Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type const&) const {}

  KOKKOS_FUNCTION
  void operator()(
      SomeTag, Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type const&) const {}
};

template <typename Policy>
void test_run_time_parameters() {
  int league_size = 131;

  using ExecutionSpace = typename Policy::execution_space;
  using ParallelTag    = Kokkos::ParallelForTag;
  int team_size =
      4 < ExecutionSpace().concurrency() ? 4 : ExecutionSpace().concurrency();
#ifdef KOKKOS_ENABLE_HPX
  team_size = 1;
#endif
#ifdef KOKKOS_ENABLE_OPENMPTARGET
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::OpenMPTarget>::value)
    team_size = 32;
#endif
  int chunk_size         = 4;
  int per_team_scratch   = 1024;
  int per_thread_scratch = 16;
  int scratch_size       = per_team_scratch + per_thread_scratch * team_size;

  Policy p1(league_size, team_size);
  ASSERT_EQ(p1.league_size(), league_size);
  ASSERT_EQ(p1.team_size(), team_size);
  ASSERT_GT(p1.chunk_size(), 0);
  ASSERT_EQ(p1.scratch_size(0), 0u);
  ASSERT_GT(p1.team_size_max(FunctorFor(), ParallelTag()), 0);
  ASSERT_GT(p1.team_size_recommended(FunctorFor(), ParallelTag()), 0);

  Policy p2 = p1.set_chunk_size(chunk_size);
  ASSERT_EQ(p1.league_size(), league_size);
  ASSERT_EQ(p1.team_size(), team_size);
  ASSERT_EQ(p1.chunk_size(), chunk_size);
  ASSERT_EQ(p1.scratch_size(0), 0u);

  ASSERT_EQ(p2.league_size(), league_size);
  ASSERT_EQ(p2.team_size(), team_size);
  ASSERT_EQ(p2.chunk_size(), chunk_size);
  ASSERT_EQ(p2.scratch_size(0), 0u);

  Policy p3 = p2.set_scratch_size(0, Kokkos::PerTeam(per_team_scratch));
  ASSERT_EQ(p2.league_size(), league_size);
  ASSERT_EQ(p2.team_size(), team_size);
  ASSERT_EQ(p2.chunk_size(), chunk_size);
  ASSERT_EQ(p2.scratch_size(0), size_t(per_team_scratch));
  ASSERT_EQ(p3.league_size(), league_size);
  ASSERT_EQ(p3.team_size(), team_size);
  ASSERT_EQ(p3.chunk_size(), chunk_size);
  ASSERT_EQ(p3.scratch_size(0), size_t(per_team_scratch));

  Policy p4 = p2.set_scratch_size(0, Kokkos::PerThread(per_thread_scratch));
  ASSERT_EQ(p2.league_size(), league_size);
  ASSERT_EQ(p2.team_size(), team_size);
  ASSERT_EQ(p2.chunk_size(), chunk_size);
  ASSERT_EQ(p2.scratch_size(0), size_t(scratch_size));
  ASSERT_EQ(p4.league_size(), league_size);
  ASSERT_EQ(p4.team_size(), team_size);
  ASSERT_EQ(p4.chunk_size(), chunk_size);
  ASSERT_EQ(p4.scratch_size(0), size_t(scratch_size));

  Policy p5 = p2.set_scratch_size(0, Kokkos::PerThread(per_thread_scratch),
                                  Kokkos::PerTeam(per_team_scratch));
  ASSERT_EQ(p2.league_size(), league_size);
  ASSERT_EQ(p2.team_size(), team_size);
  ASSERT_EQ(p2.chunk_size(), chunk_size);
  ASSERT_EQ(p2.scratch_size(0), size_t(scratch_size));
  ASSERT_EQ(p5.league_size(), league_size);
  ASSERT_EQ(p5.team_size(), team_size);
  ASSERT_EQ(p5.chunk_size(), chunk_size);
  ASSERT_EQ(p5.scratch_size(0), size_t(scratch_size));

  Policy p6 = p2.set_scratch_size(0, Kokkos::PerTeam(per_team_scratch),
                                  Kokkos::PerThread(per_thread_scratch));
  ASSERT_EQ(p2.league_size(), league_size);
  ASSERT_EQ(p2.team_size(), team_size);
  ASSERT_EQ(p2.chunk_size(), chunk_size);
  ASSERT_EQ(p2.scratch_size(0), size_t(scratch_size));
  ASSERT_EQ(p6.league_size(), league_size);
  ASSERT_EQ(p6.team_size(), team_size);
  ASSERT_EQ(p6.chunk_size(), chunk_size);
  ASSERT_EQ(p6.scratch_size(0), size_t(scratch_size));

  Policy p7 = p3.set_scratch_size(0, Kokkos::PerTeam(per_team_scratch),
                                  Kokkos::PerThread(per_thread_scratch));
  ASSERT_EQ(p3.league_size(), league_size);
  ASSERT_EQ(p3.team_size(), team_size);
  ASSERT_EQ(p3.chunk_size(), chunk_size);
  ASSERT_EQ(p3.scratch_size(0), size_t(scratch_size));
  ASSERT_EQ(p7.league_size(), league_size);
  ASSERT_EQ(p7.team_size(), team_size);
  ASSERT_EQ(p7.chunk_size(), chunk_size);
  ASSERT_EQ(p7.scratch_size(0), size_t(scratch_size));

  Policy p8;  // default constructed
  ASSERT_EQ(p8.league_size(), 0);
  ASSERT_EQ(p8.scratch_size(0), 0u);
  ASSERT_GT(p8.team_size_max(FunctorFor(), ParallelTag()), 0);
  ASSERT_GT(p8.team_size_recommended(FunctorFor(), ParallelTag()), 0);
  p8 = p3;  // call assignment operator
  ASSERT_EQ(p3.league_size(), league_size);
  ASSERT_EQ(p3.team_size(), team_size);
  ASSERT_EQ(p3.chunk_size(), chunk_size);
  ASSERT_EQ(p3.scratch_size(0), size_t(scratch_size));
  ASSERT_EQ(p8.league_size(), league_size);
  ASSERT_EQ(p8.team_size(), team_size);
  ASSERT_EQ(p8.chunk_size(), chunk_size);
  ASSERT_EQ(p8.scratch_size(0), size_t(scratch_size));

  Policy p9(league_size, Kokkos::AUTO);
  ASSERT_EQ(p9.league_size(), league_size);
  ASSERT_GT(p9.team_size_max(FunctorFor(), ParallelTag()), 0);
  ASSERT_GT(p9.team_size_recommended(FunctorFor(), ParallelTag()), 0);

  Policy p10(league_size, team_size, Kokkos::AUTO);
  ASSERT_EQ(p10.league_size(), league_size);
  ASSERT_EQ(p10.team_size(), team_size);
  ASSERT_GT(p10.team_size_max(FunctorFor(), ParallelTag()), 0);
  ASSERT_GT(p10.team_size_recommended(FunctorFor(), ParallelTag()), 0);

  Policy p11(league_size, Kokkos::AUTO, Kokkos::AUTO);
  ASSERT_EQ(p11.league_size(), league_size);
  ASSERT_GT(p11.team_size_max(FunctorFor(), ParallelTag()), 0);
  ASSERT_GT(p11.team_size_recommended(FunctorFor(), ParallelTag()), 0);
}

TEST(TEST_CATEGORY, team_policy_runtime_parameters) {
  using TestExecSpace   = TEST_EXECSPACE;
  using DynamicSchedule = Kokkos::Schedule<Kokkos::Dynamic>;
  using LongIndex       = Kokkos::IndexType<long>;

  // clang-format off
  test_run_time_parameters<Kokkos::TeamPolicy<TestExecSpace                                             >>();
  test_run_time_parameters<Kokkos::TeamPolicy<TestExecSpace,   DynamicSchedule, LongIndex               >>();
  test_run_time_parameters<Kokkos::TeamPolicy<LongIndex,       TestExecSpace,   DynamicSchedule         >>();
  test_run_time_parameters<Kokkos::TeamPolicy<DynamicSchedule, LongIndex,       TestExecSpace,   SomeTag>>();
  // clang-format on
}

}  // namespace
