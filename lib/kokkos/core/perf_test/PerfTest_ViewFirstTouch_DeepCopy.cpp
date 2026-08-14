// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "Benchmark_Context.hpp"

namespace Benchmark {

template <typename DataType>
void ViewFirstTouch_DeepCopy(benchmark::State& state) {
  const int N               = state.range(0);
  const DataType init_value = static_cast<DataType>(state.range(1));
  using ViewType            = Kokkos::View<DataType*>;
  ViewType v_a("A", N);

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::deep_copy(v_a, init_value);
    KokkosBenchmark::report_results(state, v_a, 2, timer.seconds());
  }
}

BENCHMARK_TEMPLATE(ViewFirstTouch_DeepCopy, double)
    ->ArgNames({"N", "init_value"})
    ->RangeMultiplier(8)
    ->Ranges({{int64_t(1) << 6, int64_t(1) << 24}, {0, 1}})
    ->UseManualTime();

}  // namespace Benchmark
