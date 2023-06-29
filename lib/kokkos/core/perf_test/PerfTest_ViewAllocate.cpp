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

namespace Test {

static constexpr int N = 10;

template <class Layout>
static void ViewAllocate_Rank1(benchmark::State& state) {
  const int N8 = std::pow(state.range(0), 8);

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::View<double*, Layout> a("A1", N8);
    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

template <class Layout>
static void ViewAllocate_Rank2(benchmark::State& state) {
  const int N4 = std::pow(state.range(0), 4);

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::View<double**, Layout> a("A2", N4, N4);
    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

template <class Layout>
static void ViewAllocate_Rank3(benchmark::State& state) {
  const int N2 = std::pow(state.range(0), 2);
  const int N3 = std::pow(state.range(0), 3);

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::View<double***, Layout> a("A3", N3, N3, N2);
    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

template <class Layout>
static void ViewAllocate_Rank4(benchmark::State& state) {
  const int N2 = std::pow(state.range(0), 2);

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::View<double****, Layout> a("A4", N2, N2, N2, N2);
    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

template <class Layout>
static void ViewAllocate_Rank5(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::View<double*****, Layout> a("A5", N2, N2, N1, N1, N2);
    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

template <class Layout>
static void ViewAllocate_Rank6(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::View<double******, Layout> a("A6", N2, N1, N1, N1, N1, N2);
    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

template <class Layout>
static void ViewAllocate_Rank7(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::View<double*******, Layout> a("A7", N2, N1, N1, N1, N1, N1, N1);
    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

template <class Layout>
static void ViewAllocate_Rank8(benchmark::State& state) {
  const int N1 = state.range(0);

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::View<double********, Layout> a("A8", N1, N1, N1, N1, N1, N1, N1,
                                           N1);
    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

template <class Layout>
static void ViewAllocate_Raw(benchmark::State& state) {
  const int N8 = std::pow(state.range(0), 8);
  for (auto _ : state) {
    Kokkos::Timer timer;
    double* a_ptr =
        static_cast<double*>(Kokkos::kokkos_malloc("A", sizeof(double) * N8));
    Kokkos::parallel_for(
        N8, KOKKOS_LAMBDA(const int& i) { a_ptr[i] = 0.0; });
    Kokkos::fence();
    const auto time = timer.seconds();
    Kokkos::kokkos_free(a_ptr);

    state.SetIterationTime(time);
    // data processed in megabytes
    const double data_processed = 1 * N8 * sizeof(double) / 1'000'000;
    state.counters["MB"]        = benchmark::Counter(data_processed);
    state.counters[KokkosBenchmark::benchmark_fom("GB/s")] = benchmark::Counter(
        data_processed / 1'000, benchmark::Counter::kIsIterationInvariantRate);
  }
}

BENCHMARK(ViewAllocate_Rank1<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank1<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank2<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank2<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank3<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank3<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank4<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank4<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank5<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank5<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank6<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank6<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank7<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank7<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank8<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Rank8<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

#if defined(KOKKOS_ENABLE_CUDA_LAMBDA) || !defined(KOKKOS_ENABLE_CUDA)
BENCHMARK(ViewAllocate_Raw<Kokkos::LayoutLeft>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();

BENCHMARK(ViewAllocate_Raw<Kokkos::LayoutRight>)
    ->ArgName("N")
    ->Arg(N)
    ->UseManualTime();
#endif

}  // namespace Test
