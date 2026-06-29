// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "PerfTestMDRange_Stream.hpp"

namespace Benchmark {

template <int Rank>
void MDRangePolicy_Copy(benchmark::State& state) {
  int N = static_cast<int>(state.range(0));

  MDStreamTest<TEST_EXECSPACE, Rank, double> stream_bench(N);
  stream_bench.test_copy(state);
}

template <int Rank>
void MDRangePolicy_Set(benchmark::State& state) {
  int N = static_cast<int>(state.range(0));

  MDStreamTest<TEST_EXECSPACE, Rank, double> stream_bench(N);
  stream_bench.test_set(state);
}

template <int Rank>
void MDRangePolicy_Scale(benchmark::State& state) {
  int N = static_cast<int>(state.range(0));

  MDStreamTest<TEST_EXECSPACE, Rank, double> stream_bench(N);
  stream_bench.test_scale(state);
}

template <int Rank>
void MDRangePolicy_Add(benchmark::State& state) {
  int N = static_cast<int>(state.range(0));

  MDStreamTest<TEST_EXECSPACE, Rank, double> stream_bench(N);
  stream_bench.test_add(state);
}

template <int Rank>
void MDRangePolicy_Triad(benchmark::State& state) {
  int N = static_cast<int>(state.range(0));

  MDStreamTest<TEST_EXECSPACE, Rank, double> stream_bench(N);
  stream_bench.test_triad(state);
}

// Small size for CPU backends, the problem size is computed as N^6
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL)
#define MDRANGE_BENCHMARK_ARG_SIZE 22
#else
#define MDRANGE_BENCHMARK_ARG_SIZE 16
#endif

// Macros to generate benchmarks
#define MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, RANKS) \
  BENCHMARK_TEMPLATE(BENCH_FUNCTION, RANKS)           \
      ->Arg(MDRANGE_BENCHMARK_ARG_SIZE)               \
      ->UseManualTime()                               \
      ->Unit(benchmark::kMillisecond);

#if defined(KOKKOS_ENABLE_BENCHMARKS_HEAVY)

// Generate benchmarks for ranks 1 to 6
#define MDRANGE_MAKE_BENCHMARK(BENCH_FUNCTION) \
  MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, 1)    \
  MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, 2)    \
  MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, 3)    \
  MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, 4)    \
  MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, 5)    \
  MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, 6)

MDRANGE_MAKE_BENCHMARK(MDRangePolicy_Set)
MDRANGE_MAKE_BENCHMARK(MDRangePolicy_Copy)
MDRANGE_MAKE_BENCHMARK(MDRangePolicy_Scale)
MDRANGE_MAKE_BENCHMARK(MDRangePolicy_Add)
MDRANGE_MAKE_BENCHMARK(MDRangePolicy_Triad)

#else

// Generate benchmarks for ranks 1, 3 and 6
#define MDRANGE_MAKE_BENCHMARK(BENCH_FUNCTION) \
  MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, 1)    \
  MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, 3)    \
  MDRANGE_BENCHMARK_ARGS(BENCH_FUNCTION, 6)

MDRANGE_MAKE_BENCHMARK(MDRangePolicy_Set)
MDRANGE_MAKE_BENCHMARK(MDRangePolicy_Triad)

#endif  // defined(KOKKOS_ENABLE_BENCHMARKS_HEAVY)

#undef MDRANGE_BENCHMARK_ARG_SIZE
#undef MDRANGE_BENCHMARK_ARGS
#undef MDRANGE_MAKE_BENCHMARK

}  // namespace Benchmark
