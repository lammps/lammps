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

#include <Kokkos_Core.hpp>
#include <benchmark/benchmark.h>
#include "PerfTest_Category.hpp"

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
bool is_overlapping<Kokkos::HIP>(const Kokkos::HIP&) {
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

#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_ARCH_INTEL_GPU)
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

static void OverlapRangePolicy(benchmark::State& state) {
  int N = state.range(0);
  int M = state.range(1);
  int R = state.range(2);

  TEST_EXECSPACE space;
  std::vector<TEST_EXECSPACE> execution_space_instances =
      Kokkos::Experimental::partition_space(space, 1, 1);
  TEST_EXECSPACE space1 = execution_space_instances[0];
  TEST_EXECSPACE space2 = execution_space_instances[1];

  for (auto _ : state) {
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
      KOKKOS_ASSERT(time_end > 1.5 * time_overlap);
    }
    state.counters["Time NonOverlap"] = benchmark::Counter(time_end);
    state.counters["Time Overlap"]    = benchmark::Counter(time_overlap);

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
      KOKKOS_ASSERT(time_fenced > 2.0 * time_not_fenced);
    }

    state.counters["Time fenced"]     = benchmark::Counter(time_fenced);
    state.counters["Time not fenced"] = benchmark::Counter(time_not_fenced);

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

    KOKKOS_ASSERT(h_result1() == h_result());
    KOKKOS_ASSERT(h_result2() == h_result());

    if (is_overlapping(space)) {
      KOKKOS_ASSERT(time_overlapped_reduce < 1.5 * time_no_overlapped_reduce);
    }

    state.counters["Time Reduce: NonOverlap"] =
        benchmark::Counter(time_no_overlapped_reduce);
    state.counters["Time Reduce: Overlap"] =
        benchmark::Counter(time_overlapped_reduce);
  }
}

BENCHMARK(OverlapRangePolicy)
    ->ArgNames({"N", "M", "R"})
    ->Args({2'000, 10'000, 10});

static void OverlapMDRangePolicy(benchmark::State& state) {
  int N = state.range(0);
  int M = state.range(1);
  int R = state.range(2);

  TEST_EXECSPACE space;
  std::vector<TEST_EXECSPACE> execution_space_instances =
      Kokkos::Experimental::partition_space(space, 1, 1);
  TEST_EXECSPACE space1 = execution_space_instances[0];
  TEST_EXECSPACE space2 = execution_space_instances[1];

  for (auto _ : state) {
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
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space1, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        f);
    Kokkos::parallel_for(
        "default_exec::overlap_range_policy::kernel2",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space2, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        f);
    Kokkos::fence();

    Kokkos::Timer timer;
    Kokkos::parallel_for(
        "default_exec::overlap_range_policy::kernel3",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        f);
    Kokkos::parallel_for(
        "default_exec::overlap_range_policy::kernel4",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        f);
    Kokkos::fence();

    timer.reset();
    Kokkos::parallel_for(
        "default_exec::overlap_range_policy::kernel5",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space1, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        FunctorMDRange(M, R, a));
    Kokkos::parallel_for(
        "default_exec::overlap_range_policy::kernel6",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space2, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        FunctorMDRange(M, R, a));
    Kokkos::fence();
    double time_overlap = timer.seconds();

    timer.reset();
    Kokkos::parallel_for(
        "default_exec::overlap_range_policy::kernel7",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        f);
    Kokkos::parallel_for(
        "default_exec::overlap_range_policy::kernel8",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        f);
    Kokkos::fence();
    double time_end = timer.seconds();

    if (is_overlapping(space)) {
      KOKKOS_ASSERT(time_end > 1.5 * time_overlap);
    }

    state.counters["Time NonOverlap"] = benchmark::Counter(time_end);
    state.counters["Time Overlap"]    = benchmark::Counter(time_overlap);

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
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        fr, result);
    Kokkos::fence();
    double time_fenced = timer.seconds();
    Kokkos::deep_copy(h_result, result);

    timer.reset();
    Kokkos::parallel_reduce(
        "default_exec::overlap_mdrange_policy::kernel_reduce",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        fr, result);
    double time_not_fenced = timer.seconds();
    Kokkos::fence();
    if (is_overlapping(space)) {
      KOKKOS_ASSERT(time_fenced > 2.0 * time_not_fenced);
    }

    state.counters["Time fenced"]     = benchmark::Counter(time_fenced);
    state.counters["Time not fenced"] = benchmark::Counter(time_not_fenced);

    timer.reset();
    Kokkos::parallel_reduce(
        "default_exec::overlap_mdrange_policy::kernel_reduce",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        fr, result);
    Kokkos::parallel_reduce(
        "default_exec::overlap_mdrange_policy::kernel_reduce",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        fr, result);
    Kokkos::fence();
    double time_no_overlapped_reduce = timer.seconds();

    timer.reset();
    Kokkos::parallel_reduce(
        "default_exec::overlap_mdrange_policy::kernel_reduce",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space1, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        fr, result1);
    Kokkos::parallel_reduce(
        "default_exec::overlap_mdrange_policy::kernel_reduce",
        Kokkos::Experimental::require(
            Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>(
                space2, {0, 0}, {N, R}),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        fr, result2);
    Kokkos::fence();
    double time_overlapped_reduce = timer.seconds();

    Kokkos::deep_copy(h_result2, result2);
    Kokkos::deep_copy(h_result1, result1);

    KOKKOS_ASSERT(h_result1() == h_result());
    KOKKOS_ASSERT(h_result2() == h_result());

    if (is_overlapping(space)) {
      KOKKOS_ASSERT(time_overlapped_reduce < 1.5 * time_no_overlapped_reduce);
    }

    state.counters["Time Reduce: NonOverlap"] =
        benchmark::Counter(time_no_overlapped_reduce);
    state.counters["Time Reduce: Time Overlap"] =
        benchmark::Counter(time_overlapped_reduce);
  }
}

BENCHMARK(OverlapMDRangePolicy)
    ->ArgNames({"N", "M", "R"})
    ->Args({200, 10'000, 10});

static void OverlapTeamPolicy(benchmark::State& state) {
  int N = state.range(0);
  int M = state.range(1);
  int R = state.range(2);

  TEST_EXECSPACE space;
  std::vector<TEST_EXECSPACE> execution_space_instances =
      Kokkos::Experimental::partition_space(space, 1, 1);
  TEST_EXECSPACE space1 = execution_space_instances[0];
  TEST_EXECSPACE space2 = execution_space_instances[1];

  for (auto _ : state) {
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
      KOKKOS_ASSERT(time_end > 1.5 * time_overlap);
    }

    state.counters["Time NonOverlap"] = benchmark::Counter(time_end);
    state.counters["Time Overlap"]    = benchmark::Counter(time_overlap);

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
      KOKKOS_ASSERT(time_fenced > 2.0 * time_not_fenced);
    }

    state.counters["Time fenced"]     = benchmark::Counter(time_fenced);
    state.counters["Time not fenced"] = benchmark::Counter(time_not_fenced);

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

    KOKKOS_ASSERT(h_result1() == h_result());
    KOKKOS_ASSERT(h_result2() == h_result());

    if (is_overlapping(space)) {
      KOKKOS_ASSERT(time_overlapped_reduce < 1.5 * time_no_overlapped_reduce);
    }

    state.counters["Time Reduce: NonOverlap"] =
        benchmark::Counter(time_no_overlapped_reduce);
    state.counters["Time Reduce: Time Overlap"] =
        benchmark::Counter(time_overlapped_reduce);
  }
}

BENCHMARK(OverlapTeamPolicy)
    ->ArgNames({"N", "M", "R"})
    ->Args({20, 1'000'000, 10});

}  // namespace Test
