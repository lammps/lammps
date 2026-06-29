// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "Benchmark_Context.hpp"

namespace Benchmark {

template <typename DataType>
void ViewFirstTouch_ParallelFor(benchmark::State& state) {
  const int N    = state.range(0);
  using ViewType = Kokkos::View<DataType*>;

  for (auto _ : state) {
    ViewType v_a("A", N);
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::parallel_for(
        "ViewFirstTouch_ParallelFor", N, KOKKOS_LAMBDA(const int i) {
          v_a(i) = static_cast<DataType>(2) * v_a(i) + static_cast<DataType>(1);
        });
    Kokkos::fence();
    KokkosBenchmark::report_results(state, v_a, 2, timer.seconds());
  }
}

BENCHMARK_TEMPLATE(ViewFirstTouch_ParallelFor, double)
    ->ArgName("N")
    ->RangeMultiplier(8)
    ->Range(int64_t(1) << 6, int64_t(1) << 24)
    ->UseManualTime();

}  // namespace Benchmark
