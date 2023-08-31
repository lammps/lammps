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
#include "Benchmark_Context.hpp"
#include "PerfTest_Category.hpp"
#include <Kokkos_Random.hpp>
#include <utility>

#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
namespace Test {
template <class Scalar>
std::pair<double, Scalar> custom_reduction_test(int N, int R) {
  Kokkos::Random_XorShift64_Pool<> rand_pool(183291);
  Kokkos::View<Scalar*> a("A", N);
  Kokkos::fill_random(a, rand_pool, 1.0);

  Scalar max;

  int team_size = 32;
  if (team_size > Kokkos::DefaultExecutionSpace().concurrency())
    team_size = Kokkos::DefaultExecutionSpace().concurrency();
  // Warm up
  Kokkos::parallel_reduce(
      Kokkos::TeamPolicy<>(N / 1024, team_size),
      KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team,
                    Scalar& lmax) {
        Scalar team_max = Scalar(0);
        for (int rr = 0; rr < R; rr++) {
          int i = team.league_rank();
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(team, 32),
              [&](const int& j, Scalar& thread_max) {
                Scalar t_max = Scalar(0);
                Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange(team, 32),
                    [&](const int& k, Scalar& max_) {
                      const Scalar val = a((i * 32 + j) * 32 + k);
                      if (val > max_) max_ = val;
                      if ((k == 11) && (j == 17) && (i == 2)) max_ = 11.5;
                    },
                    Kokkos::Max<Scalar>(t_max));
                if (t_max > thread_max) thread_max = t_max;
              },
              Kokkos::Max<Scalar>(team_max));
        }
        if (team_max > lmax) lmax = team_max;
      },
      Kokkos::Max<Scalar>(max));

  // Timing
  Kokkos::Timer timer;
  Kokkos::parallel_reduce(
      Kokkos::TeamPolicy<>(N / 1024, team_size),
      KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team,
                    Scalar& lmax) {
        Scalar team_max = Scalar(0);
        for (int rr = 0; rr < R; rr++) {
          int i = team.league_rank();
          Kokkos::parallel_reduce(
              Kokkos::TeamThreadRange(team, 32),
              [&](const int& j, Scalar& thread_max) {
                Scalar t_max = Scalar(0);
                Kokkos::parallel_reduce(
                    Kokkos::ThreadVectorRange(team, 32),
                    [&](const int& k, Scalar& max_) {
                      const Scalar val = a((i * 32 + j) * 32 + k);
                      if (val > max_) max_ = val;
                      if ((k == 11) && (j == 17) && (i == 2)) max_ = 11.5;
                    },
                    Kokkos::Max<Scalar>(t_max));
                if (t_max > thread_max) thread_max = t_max;
              },
              Kokkos::Max<Scalar>(team_max));
        }
        if (team_max > lmax) lmax = team_max;
      },
      Kokkos::Max<Scalar>(max));

  return std::make_pair(timer.seconds(), max);
}

int get_N(benchmark::State& state) {
  return (Test::command_line_num_args() > 1)
             ? std::stoi(Test::command_line_arg(1))
             : state.range(0);
}

int get_R(benchmark::State& state) {
  return (Test::command_line_num_args() > 2)
             ? std::stoi(Test::command_line_arg(2))
             : state.range(1);
}

template <class Scalar>
static void CustomReduction(benchmark::State& state) {
  int N = get_N(state);
  int R = get_R(state);

  for (auto _ : state) {
    auto results = custom_reduction_test<double>(N, R);
    // data processed in gigabytes
    const double data_processed =
        N * R * sizeof(Scalar) / results.first / 1'000'000'000;

    state.SetIterationTime(results.first);
    state.counters[KokkosBenchmark::benchmark_fom("GB/s")] = benchmark::Counter(
        data_processed, benchmark::Counter::kIsIterationInvariantRate);
    state.counters["Max"] = benchmark::Counter(results.second);
  }
}

BENCHMARK(CustomReduction<double>)
    ->ArgNames({"N", "R"})
    ->Args({100'000, 1'000})
    ->UseManualTime();

}  // namespace Test
#endif
