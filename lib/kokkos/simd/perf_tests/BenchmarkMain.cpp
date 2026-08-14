// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <benchmark/benchmark.h>
#include <Kokkos_Core.hpp>

#include <Benchmark_Context.hpp>
#if defined(KOKKOS_IMPL_SIMD_HOST_PERFTEST)
#include <PerfTest_Host.hpp>
#else
#include <PerfTest_Device.hpp>
#endif

int main(int argc, char** argv) {
  Kokkos::ScopeGuard guard(argc, argv);

  {
#if defined(KOKKOS_IMPL_SIMD_HOST_PERFTEST)
    register_host_benchmarks();
#else
    register_device_benchmarks();
#endif

    benchmark::Initialize(&argc, argv);
    benchmark::SetDefaultTimeUnit(benchmark::kMillisecond);
    KokkosBenchmark::add_benchmark_context(true);

    benchmark::RunSpecifiedBenchmarks();

    benchmark::Shutdown();
  }

  return 0;
}
