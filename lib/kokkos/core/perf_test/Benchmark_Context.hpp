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

#ifndef KOKKOS_CORE_PERFTEST_BENCHMARK_CONTEXT_HPP
#define KOKKOS_CORE_PERFTEST_BENCHMARK_CONTEXT_HPP

#include <string>

#include <benchmark/benchmark.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Version_Info.hpp>

namespace KokkosBenchmark {

/**
 * \brief Gather all context information and add it to benchmark context data
 */
void add_benchmark_context(bool verbose = false);

/**
 * \brief Mark the label as a figure of merit.
 */
inline std::string benchmark_fom(const std::string& label) {
  return "FOM: " + label;
}

/**
 * \brief Report throughput and amount of data processed for simple View
 * operations
 */
template <class ViewType>
void report_results(benchmark::State& state, ViewType view, int data_ratio,
                    double time) {
  // data processed in megabytes
  const double data_processed = data_ratio * view.size() *
                                sizeof(typename ViewType::value_type) /
                                1'000'000;

  state.SetIterationTime(time);
  state.counters["MB"] = benchmark::Counter(data_processed);
  state.counters[KokkosBenchmark::benchmark_fom("GB/s")] = benchmark::Counter(
      data_processed / 1'000, benchmark::Counter::kIsIterationInvariantRate);
}

}  // namespace KokkosBenchmark

#endif
