// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include <benchmark/benchmark.h>
#include "Benchmark_Context.hpp"

// Returns time for reduction in us
template <class F, class R>
void run_reduction(int N, F f, R& r) {
  Kokkos::parallel_reduce(N, f, r);
  // fence because r could be a View and thus
  // parallel_reduce can be asynchronous
  Kokkos::fence();
}

template <class T, size_t TestCase>
static void Test_Reduction(benchmark::State& state) {
  int N = state.range(0);

  // Prep some Views used for benchmark
  // avoids allocating/deallocating during benchmark phase
  Kokkos::View<T*> data("data", N);
  [[maybe_unused]] Kokkos::View<T> d_result("d_result");
  [[maybe_unused]] Kokkos::View<T, Kokkos::HostSpace> h_result("h_result");
  Kokkos::parallel_for(
      "FillData", N, KOKKOS_LAMBDA(int i) { data(i) = i; });

  // Simple kernel with little work
  auto f_simple = KOKKOS_LAMBDA(int i, T& val) { val += data[i]; };

  // Expensive kernel with more work per thread
  auto f_expensive = KOKKOS_LAMBDA(int i, T& val) {
    T result = data(i);
    for (int k = 0; k < 100; k++) {
      result = Kokkos::log(result);
      result = Kokkos::sqrt(result);
      result = Kokkos::exp(result);
      result = result * result;
    }
    val += result;
  };

  [[maybe_unused]] T s_result = 0;

  // Warmup to get scratch sizes set up properly
  if constexpr (TestCase == 0) run_reduction(N, f_simple, s_result);
  if constexpr (TestCase == 1) run_reduction(N, f_simple, h_result);
  if constexpr (TestCase == 2) run_reduction(N, f_simple, d_result);
  if constexpr (TestCase == 3) run_reduction(N, f_expensive, s_result);
  Kokkos::fence();

  for (auto _ : state) {
    // Simple kernel with scalar as result
    if constexpr (TestCase == 0) run_reduction(N, f_simple, s_result);
    // Simple kernel with host view as result
    if constexpr (TestCase == 1) run_reduction(N, f_simple, h_result);
    // Simple kernel with device view as result
    if constexpr (TestCase == 2) run_reduction(N, f_simple, d_result);

    // Expensive kernel with scalar result
    // Adding views as result doesn't add anything here since the runtime
    // outweighs the copy of the result
    if constexpr (TestCase == 3) run_reduction(N, f_expensive, s_result);
  }
}

template <class T>
static void ReductionCheapScalarResult(benchmark::State& state) {
  Test_Reduction<T, 0>(state);
}

template <class T>
static void ReductionCheapHostResult(benchmark::State& state) {
  Test_Reduction<T, 1>(state);
}

template <class T>
static void ReductionCheapDeviceResult(benchmark::State& state) {
  Test_Reduction<T, 2>(state);
}

template <class T>
static void ReductionExpensive(benchmark::State& state) {
  Test_Reduction<T, 3>(state);
}

BENCHMARK(ReductionCheapScalarResult<double>)
    ->Iterations(10)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(10)
    ->Range(10'000, 10'000'000);
BENCHMARK(ReductionCheapHostResult<double>)
    ->Iterations(10)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(10)
    ->Range(10'000, 10'000'000);
BENCHMARK(ReductionCheapDeviceResult<double>)
    ->Iterations(10)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(10)
    ->Range(10'000, 10'000'000);
BENCHMARK(ReductionExpensive<double>)
    ->Iterations(10)
    ->Unit(benchmark::kMicrosecond)
    ->RangeMultiplier(10)
    ->Range(10'000, 10'000'000);
