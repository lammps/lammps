#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>
#include <PerfTest_Category.hpp>

namespace Test {

namespace {

template <class ExecutionSpace>
bool is_overlapping(const ExecutionSpace&) {
  return false;
}

#ifndef KOKKOS_ENABLE_DEBUG
#ifdef KOKKOS_ENABLE_CUDA
template <>
bool is_overlapping<Kokkos::Cuda>(const Kokkos::Cuda&) {
  bool value          = true;
  auto local_rank_str = std::getenv("CUDA_LAUNCH_BLOCKING");
  if (local_rank_str) {
    value = (std::stoi(local_rank_str) == 0);
  }
  return value;
}
#endif

#ifdef KOKKOS_ENABLE_HIP
template <>
bool is_overlapping<Kokkos::Experimental::HIP>(
    const Kokkos::Experimental::HIP&) {
  // FIXME_HIP This doesn't pass yet in CI.
  return false;
  // bool value          = true;
  // auto local_rank_str = std::getenv("HIP_LAUNCH_BLOCKING");
  // if (local_rank_str) {
  //  value = (std::stoi(local_rank_str) == 0);
  //}
  // return value;
}
#endif

#ifdef KOKKOS_ENABLE_SYCL
template <>
bool is_overlapping<Kokkos::Experimental::SYCL>(
    const Kokkos::Experimental::SYCL&) {
  return true;
}
#endif
#endif

}  // namespace

struct FunctorRange {
  int M, R;
  Kokkos::View<double**, TEST_EXECSPACE> a;
  FunctorRange(int M_, int R_, Kokkos::View<double**, TEST_EXECSPACE> a_)
      : M(M_), R(R_), a(a_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    for (int r = 0; r < R; r++)
      for (int j = 0; j < M; j++) {
        a(i, j) += 1.0;
      }
  }
};

struct FunctorMDRange {
  int M, R;
  Kokkos::View<double**, TEST_EXECSPACE> a;
  FunctorMDRange(int M_, int R_, Kokkos::View<double**, TEST_EXECSPACE> a_)
      : M(M_), R(R_), a(a_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int) const {
    for (int j = 0; j < M; j++) a(i, j) += 1.0;
  }
};

struct FunctorTeam {
  int M, R;
  Kokkos::View<double**, Kokkos::LayoutRight, TEST_EXECSPACE> a;
  FunctorTeam(int M_, int R_,
              Kokkos::View<double**, Kokkos::LayoutRight, TEST_EXECSPACE> a_)
      : M(M_), R(R_), a(a_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(
      const Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type& team) const {
    int i = team.league_rank();
    for (int r = 0; r < R; r++) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, M),
                           [&](const int j) { a(i, j) += 1.0; });
    }
  }
};

struct FunctorRangeReduce {
  int M, R;
  Kokkos::View<double**, TEST_EXECSPACE> a;
  FunctorRangeReduce(int M_, int R_, Kokkos::View<double**, TEST_EXECSPACE> a_)
      : M(M_), R(R_), a(a_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, double& tmp) const {
    for (int r = 0; r < R; r++)
      for (int j = 0; j < M; j++) {
        tmp += a(i, j);
      }
  }
};

struct FunctorMDRangeReduce {
  int M, R;
  Kokkos::View<double**, TEST_EXECSPACE> a;
  FunctorMDRangeReduce(int M_, int R_,
                       Kokkos::View<double**, TEST_EXECSPACE> a_)
      : M(M_), R(R_), a(a_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int, double& tmp) const {
    for (int j = 0; j < M; j++) tmp += a(i, j);
  }
};

struct FunctorTeamReduce {
  int M, R;
  Kokkos::View<double**, Kokkos::LayoutRight, TEST_EXECSPACE> a;
  FunctorTeamReduce(
      int M_, int R_,
      Kokkos::View<double**, Kokkos::LayoutRight, TEST_EXECSPACE> a_)
      : M(M_), R(R_), a(a_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const Kokkos::TeamPolicy<TEST_EXECSPACE>::member_type& team,
                  double& tmp) const {
    int i = team.league_rank();
    for (int r = 0; r < R; r++) {
      double val;
      Kokkos::parallel_reduce(
          Kokkos::TeamThreadRange(team, M),
          [&](const int j, double& tmp2) { tmp2 += a(i, j); }, val);
      tmp += val;
    }
  }
};

TEST(default_exec, overlap_range_policy) {
  int N = 2000;
  int M = 10000;
  int R = 10;

  TEST_EXECSPACE space;
  std::vector<TEST_EXECSPACE> execution_space_instances =
      Kokkos::Experimental::partition_space(space, 1, 1);
  TEST_EXECSPACE space1 = execution_space_instances[0];
  TEST_EXECSPACE space2 = execution_space_instances[1];

  Kokkos::View<double**, TEST_EXECSPACE> a("A", N, M);
  FunctorRange f(M, R, a);
  FunctorRangeReduce fr(M, R, a);
  Kokkos::parallel_for("default_exec::overlap_range_policy::kernel0",
                       Kokkos::RangePolicy<TEST_EXECSPACE>(0, N),
                       FunctorRange(M, R, a));

  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel1",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space1, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel2",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space2, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::fence();

  Kokkos::Timer timer;
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel3",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel4",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel5",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space1, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorRange(M, R, a));
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel6",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space2, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorRange(M, R, a));
  Kokkos::fence();
  double time_overlap = timer.seconds();

  timer.reset();
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel7",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel8",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::fence();
  double time_end = timer.seconds();

  if (is_overlapping(space)) {
    ASSERT_GT(time_end, 1.5 * time_overlap);
  }
  printf("Time RangePolicy: NonOverlap: %lf Time Overlap: %lf\n", time_end,
         time_overlap);

  Kokkos::View<double, TEST_EXECSPACE> result("result");
  Kokkos::View<double, TEST_EXECSPACE> result1("result1");
  Kokkos::View<double, TEST_EXECSPACE> result2("result2");
  Kokkos::View<double, Kokkos::HostSpace> h_result("h_result");
  Kokkos::View<double, Kokkos::HostSpace> h_result1("h_result1");
  Kokkos::View<double, Kokkos::HostSpace> h_result2("h_result2");

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_range_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  Kokkos::fence();
  double time_fenced = timer.seconds();
  Kokkos::deep_copy(h_result, result);

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_range_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  double time_not_fenced = timer.seconds();
  Kokkos::fence();
  if (is_overlapping(space)) {
    ASSERT_GT(time_fenced, 2.0 * time_not_fenced);
  }

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_range_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  Kokkos::parallel_reduce(
      "default_exec::overlap_range_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  Kokkos::fence();
  double time_no_overlapped_reduce = timer.seconds();

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_range_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space1, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result1);
  Kokkos::parallel_reduce(
      "default_exec::overlap_range_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<TEST_EXECSPACE>(space2, 0, N),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result2);
  Kokkos::fence();
  double time_overlapped_reduce = timer.seconds();

  Kokkos::deep_copy(h_result2, result2);
  Kokkos::deep_copy(h_result1, result1);

  ASSERT_EQ(h_result1(), h_result());
  ASSERT_EQ(h_result2(), h_result());

  if (is_overlapping(space)) {
    ASSERT_LT(time_overlapped_reduce, 1.5 * time_no_overlapped_reduce);
  }
  printf("Time RangePolicy Reduce: NonOverlap: %lf Time Overlap: %lf\n",
         time_no_overlapped_reduce, time_overlapped_reduce);
}

TEST(default_exec, overlap_mdrange_policy) {
  int N = 200;
  int M = 10000;
  int R = 10;

  TEST_EXECSPACE space;
  std::vector<TEST_EXECSPACE> execution_space_instances =
      Kokkos::Experimental::partition_space(space, 1, 1);
  TEST_EXECSPACE space1 = execution_space_instances[0];
  TEST_EXECSPACE space2 = execution_space_instances[1];

  Kokkos::View<double**, TEST_EXECSPACE> a("A", N, M);
  FunctorMDRange f(M, R, a);
  FunctorMDRangeReduce fr(M, R, a);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel0",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>({0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorMDRange(M, R, a));

  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel1",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space1, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel2",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space2, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::fence();

  Kokkos::Timer timer;
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel3",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel4",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel5",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space1, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorMDRange(M, R, a));
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel6",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space2, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorMDRange(M, R, a));
  Kokkos::fence();
  double time_overlap = timer.seconds();

  timer.reset();
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel7",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel8",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::fence();
  double time_end = timer.seconds();

  if (is_overlapping(space)) {
    ASSERT_GT(time_end, 1.5 * time_overlap);
  }
  printf("Time MDRangePolicy: NonOverlap: %lf Time Overlap: %lf\n", time_end,
         time_overlap);

  Kokkos::View<double, TEST_EXECSPACE> result("result");
  Kokkos::View<double, TEST_EXECSPACE> result1("result1");
  Kokkos::View<double, TEST_EXECSPACE> result2("result2");
  Kokkos::View<double, Kokkos::HostSpace> h_result("h_result");
  Kokkos::View<double, Kokkos::HostSpace> h_result1("h_result1");
  Kokkos::View<double, Kokkos::HostSpace> h_result2("h_result2");

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_mdrange_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  Kokkos::fence();
  double time_fenced = timer.seconds();
  Kokkos::deep_copy(h_result, result);

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_mdrange_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  double time_not_fenced = timer.seconds();
  Kokkos::fence();
  if (is_overlapping(space)) {
    ASSERT_GT(time_fenced, 2.0 * time_not_fenced);
  }

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_mdrange_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  Kokkos::parallel_reduce(
      "default_exec::overlap_mdrange_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  Kokkos::fence();
  double time_no_overlapped_reduce = timer.seconds();

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_mdrange_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space1, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result1);
  Kokkos::parallel_reduce(
      "default_exec::overlap_mdrange_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(space2, {0, 0},
                                                                 {N, R}),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result2);
  Kokkos::fence();
  double time_overlapped_reduce = timer.seconds();

  Kokkos::deep_copy(h_result2, result2);
  Kokkos::deep_copy(h_result1, result1);

  ASSERT_EQ(h_result1(), h_result());
  ASSERT_EQ(h_result2(), h_result());

  if (is_overlapping(space)) {
    ASSERT_LT(time_overlapped_reduce, 1.5 * time_no_overlapped_reduce);
  }
  printf("Time MDRangePolicy Reduce: NonOverlap: %lf Time Overlap: %lf\n",
         time_no_overlapped_reduce, time_overlapped_reduce);
}

TEST(default_exec, overlap_team_policy) {
  int N = 20;
  int M = 1000000;
  int R = 10;

  TEST_EXECSPACE space;
  std::vector<TEST_EXECSPACE> execution_space_instances =
      Kokkos::Experimental::partition_space(space, 1, 1);
  TEST_EXECSPACE space1 = execution_space_instances[0];
  TEST_EXECSPACE space2 = execution_space_instances[1];

  Kokkos::View<double**, Kokkos::LayoutRight, TEST_EXECSPACE> a("A", N, M);
  FunctorTeam f(M, R, a);
  FunctorTeamReduce fr(M, R, a);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel0",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorTeam(M, R, a));

  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel1",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space1, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel2",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space2, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::fence();

  Kokkos::Timer timer;
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel3",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel4",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::fence();

  timer.reset();
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel5",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space1, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorTeam(M, R, a));
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel6",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space2, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorTeam(M, R, a));
  Kokkos::fence();
  double time_overlap = timer.seconds();

  timer.reset();
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel7",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::parallel_for(
      "default_exec::overlap_range_policy::kernel8",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      f);
  Kokkos::fence();
  double time_end = timer.seconds();

  if (is_overlapping(space)) {
    ASSERT_GT(time_end, 1.5 * time_overlap);
  }
  printf("Time TeamPolicy: NonOverlap: %lf Time Overlap: %lf\n", time_end,
         time_overlap);

  Kokkos::View<double, TEST_EXECSPACE> result("result");
  Kokkos::View<double, TEST_EXECSPACE> result1("result1");
  Kokkos::View<double, TEST_EXECSPACE> result2("result2");
  Kokkos::View<double, Kokkos::HostSpace> h_result("h_result");
  Kokkos::View<double, Kokkos::HostSpace> h_result1("h_result1");
  Kokkos::View<double, Kokkos::HostSpace> h_result2("h_result2");

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_team_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  Kokkos::fence();
  double time_fenced = timer.seconds();
  Kokkos::deep_copy(h_result, result);

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_team_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  double time_not_fenced = timer.seconds();
  Kokkos::fence();
  if (is_overlapping(space)) {
    ASSERT_GT(time_fenced, 2.0 * time_not_fenced);
  }
  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_team_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  Kokkos::parallel_reduce(
      "default_exec::overlap_team_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result);
  Kokkos::fence();
  double time_no_overlapped_reduce = timer.seconds();

  timer.reset();
  Kokkos::parallel_reduce(
      "default_exec::overlap_team_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space1, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result1);
  Kokkos::parallel_reduce(
      "default_exec::overlap_team_policy::kernel_reduce",
      Kokkos::Experimental::require(
          Kokkos::TeamPolicy<TEST_EXECSPACE>(space2, N, Kokkos::AUTO),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      fr, result2);
  Kokkos::fence();
  double time_overlapped_reduce = timer.seconds();

  Kokkos::deep_copy(h_result2, result2);
  Kokkos::deep_copy(h_result1, result1);

  ASSERT_EQ(h_result1(), h_result());
  ASSERT_EQ(h_result2(), h_result());

  if (is_overlapping(space)) {
    ASSERT_LT(time_overlapped_reduce, 1.5 * time_no_overlapped_reduce);
  }
  printf("Time TeamPolicy Reduce: NonOverlap: %lf Time Overlap: %lf\n",
         time_no_overlapped_reduce, time_overlapped_reduce);
}
}  // namespace Test
