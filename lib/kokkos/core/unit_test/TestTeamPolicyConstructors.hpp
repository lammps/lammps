// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif

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

// The execution space is defaulted if not given to the constructor.
TEST(TEST_CATEGORY, team_policy_default_space) {
  using policy_t = Kokkos::TeamPolicy<TEST_EXECSPACE>;

  policy_t defaulted(42, Kokkos::AUTO);

  ASSERT_EQ(defaulted.space(), TEST_EXECSPACE{});
}

// The execution space instance can be updated.
TEST(TEST_CATEGORY, team_policy_impl_set_space) {
  using policy_t = Kokkos::TeamPolicy<TEST_EXECSPACE>;

  const auto [exec_old, exec_new] =
      Kokkos::Experimental::partition_space(TEST_EXECSPACE{}, 1, 1);

  const policy_t policy_old(exec_old, 42, Kokkos::AUTO);
  ASSERT_EQ(policy_old.space(), exec_old);

  const policy_t policy_new(Kokkos::Impl::PolicyUpdate{}, policy_old, exec_new);

  ASSERT_EQ(policy_new.space(), exec_new);
  ASSERT_EQ(policy_new.league_size(), 42);
}

}  // namespace
