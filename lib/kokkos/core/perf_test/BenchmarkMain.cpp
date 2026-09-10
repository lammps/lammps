// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <benchmark/benchmark.h>

#include "Benchmark_Context.hpp"
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include "PerfTest_Category.hpp"

int main(int argc, char** argv) {
  Kokkos::initialize(argc, argv);
  benchmark::Initialize(&argc, argv);
  // FIXME: seconds as default time unit leads to precision loss
  benchmark::SetDefaultTimeUnit(benchmark::kSecond);
  KokkosBenchmark::add_benchmark_context(true);

  (void)Test::command_line_num_args(argc);
  (void)Test::command_line_arg(0, argv);

  benchmark::RunSpecifiedBenchmarks();

  benchmark::Shutdown();
  Kokkos::finalize();
  return 0;
}
