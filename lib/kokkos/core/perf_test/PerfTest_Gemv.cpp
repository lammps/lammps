// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/**
 * @file PerfTest_Gemv.cpp
 *
 * This file implements a performance test of a naive GEMV implementation.
 * This was created to reproduce a performance problem observed on MI300A using
 * the HIP backend. With a relatively small number of rows (4000), the HIP
 * backend was choosing large block sizes, which led to very few active EUs on
 * the GPU.
 */

#include <Kokkos_Core.hpp>
#include <benchmark/benchmark.h>
#include "Benchmark_Context.hpp"

namespace Benchmark {

template <typename P, typename L>
void impl(benchmark::State& state) {
  const size_t M = state.range(0);
  const size_t N = state.range(1);

  Kokkos::View<int**, L> A("A", M, N);
  Kokkos::View<int*> y("y", M);
  Kokkos::View<int*> x("x", N);

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::parallel_for(
        P(0, y.extent(0)), KOKKOS_LAMBDA(const int i) {
          int sum = 0;
          for (int j = 0; j < x.extent_int(0); j++) {
            sum += A(i, j) * x(j);
          }
          y(i) = sum;
        });
    Kokkos::fence();
    state.SetIterationTime(timer.seconds());
  }
}

template <typename L>
static void GemvDefault(benchmark::State& state) {
  using P = Kokkos::RangePolicy<>;
  impl<P, L>(state);
}

template <int BS, typename L>
static void Gemv(benchmark::State& state) {
  using P = Kokkos::RangePolicy<Kokkos::LaunchBounds<BS, 1>>;
  impl<P, L>(state);
}

#define COMMON_ARGS()                 \
  ArgNames({"M", "N"})                \
      ->UseManualTime()               \
      ->Unit(benchmark::kMicrosecond) \
      ->Args({4'000, 10'000})         \
      ->Args({400'000, 100})

BENCHMARK(GemvDefault<Kokkos::LayoutLeft>)->COMMON_ARGS();
BENCHMARK(GemvDefault<Kokkos::LayoutRight>)->COMMON_ARGS();

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
BENCHMARK(Gemv<64, Kokkos::LayoutLeft>)->COMMON_ARGS();
BENCHMARK(Gemv<64, Kokkos::LayoutRight>)->COMMON_ARGS();
BENCHMARK(Gemv<256, Kokkos::LayoutLeft>)->COMMON_ARGS();
BENCHMARK(Gemv<256, Kokkos::LayoutRight>)->COMMON_ARGS();
BENCHMARK(Gemv<1024, Kokkos::LayoutLeft>)->COMMON_ARGS();
BENCHMARK(Gemv<1024, Kokkos::LayoutRight>)->COMMON_ARGS();
#endif

#undef COMMON_ARGS

}  // namespace Benchmark
