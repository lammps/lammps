// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <benchmark/benchmark.h>
// FIXME: Benchmark_Context.hpp should be moved to a common location
#include "../../core/perf_test/Benchmark_Context.hpp"

namespace Benchmark {

// Fills each entry of
// * a view of size N
// * with the sum of K random numbers
template <typename Pool>
static void Random(benchmark::State &state) {
  const size_t N     = state.range(0);
  const size_t K     = state.range(1);
  constexpr double I = 1024;

  Kokkos::View<double *> out("out", N);
  Pool random_pool(/*seed=*/12345);

  for ([[maybe_unused]] auto _ : state) {
    Kokkos::parallel_for(
        N, KOKKOS_LAMBDA(const int i) {
          auto generator = random_pool.get_state();
          double acc     = 0;

          for (size_t k = 0; k < K; ++k) {
            acc += generator.drand(I);
          }
          random_pool.free_state(generator);
          out(i) = acc;
        });
    Kokkos::fence();
  }

  state.counters[KokkosBenchmark::benchmark_fom("rate")] = benchmark::Counter(
      state.iterations() * N * K, benchmark::Counter::kIsRate);
}

static void Random64(benchmark::State &state) {
  return Random<Kokkos::Random_XorShift64_Pool<>>(state);
}

static void Random1024(benchmark::State &state) {
  return Random<Kokkos::Random_XorShift1024_Pool<>>(state);
}

#define RANDOM_ARGS()                      \
  ArgNames({"N", "K"})                     \
      ->ArgsProduct({{1 << 21}, {1, 256}}) \
      ->UseRealTime()                      \
      ->Unit(benchmark::kMicrosecond)

BENCHMARK(Random64)->RANDOM_ARGS();
BENCHMARK(Random1024)->RANDOM_ARGS();

#undef RANDOM_ARGS

}  // namespace Benchmark
