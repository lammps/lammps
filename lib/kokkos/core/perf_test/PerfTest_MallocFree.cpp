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

namespace Benchmark {

// when the time will be recorded
enum class When { after_malloc, after_touch, after_free };

static void Impl(benchmark::State& state, const bool touch, const When when) {
  const size_t N = state.range(0);
  for (auto _ : state) {
    Kokkos::Timer timer;
    char* a_ptr = static_cast<char*>(Kokkos::kokkos_malloc("A", N));
    if (When::after_malloc == when) {
      state.SetIterationTime(timer.seconds());
    }
    if (touch) {
      constexpr size_t STRIDE = 1024;  // stride for touching the allocation.
      // this is intended to be a safe value that would touch every "page", but
      // not saturate the memory bandwidth
      Kokkos::parallel_for(
          N / STRIDE,
          KOKKOS_LAMBDA(const size_t& i) { a_ptr[i * STRIDE] = i * STRIDE; });
      Kokkos::fence();
    }
    if (When::after_touch == when) {
      state.SetIterationTime(timer.seconds());
    }
    Kokkos::kokkos_free(a_ptr);
    if (When::after_free == when) {
      state.SetIterationTime(timer.seconds());
    }
  }

  state.counters[KokkosBenchmark::benchmark_fom("rate")] =
      benchmark::Counter(state.iterations(), benchmark::Counter::kIsRate);
}

static void Malloc(benchmark::State& state) {
  Impl(state, false, When::after_malloc);
}

static void MallocFree(benchmark::State& state) {
  Impl(state, false, When::after_free);
}

static void MallocTouch(benchmark::State& state) {
  Impl(state, true, When::after_touch);
}

static void MallocTouchFree(benchmark::State& state) {
  Impl(state, true, When::after_free);
}

#ifdef KOKKOS_IMPL_32BIT
constexpr int test_range = 30;
#else
#ifndef KOKKOS_ENABLE_LARGE_MEM_TESTS
constexpr int test_range = 31;
#else
constexpr int test_range = 32;
#endif
#endif

BENCHMARK(Malloc)
    ->ArgName("N")
    ->RangeMultiplier(16)
    ->Range(1, int64_t(1) << test_range)
    ->UseManualTime()
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(MallocFree)
    ->ArgName("N")
    ->RangeMultiplier(16)
    ->Range(1, int64_t(1) << test_range)
    ->UseManualTime()
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(MallocTouch)
    ->ArgName("N")
    ->RangeMultiplier(16)
    ->Range(1, int64_t(1) << test_range)
    ->UseManualTime()
    ->Unit(benchmark::kMicrosecond);

BENCHMARK(MallocTouchFree)
    ->ArgName("N")
    ->RangeMultiplier(16)
    ->Range(1, int64_t(1) << test_range)
    ->UseManualTime()
    ->Unit(benchmark::kMicrosecond);

}  // namespace Benchmark
