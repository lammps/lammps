// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <type_traits>

#include <benchmark/benchmark.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.std_algorithms;
#else
#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#endif
#include <Kokkos_Timer.hpp>
// FIXME: Benchmark_Context.hpp should be moved to a common location
#include "../../core/perf_test/Benchmark_Context.hpp"

namespace {

namespace KE = Kokkos::Experimental;

using ExecSpace     = Kokkos::DefaultExecutionSpace;
using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

// A tag struct to identify when inclusive scan with the implicit sum
// based binary operation needs to be called.
template <class ValueType>
struct ImpSumBinOp;

template <class ValueType>
struct SumFunctor {
  KOKKOS_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return (a + b);
  }
};

template <class ValueType>
struct MaxFunctor {
  KOKKOS_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    if (a > b)
      return a;
    else
      return b;
  }
};

// Helper to obtain last element of a view
template <class T>
T obtain_last_elem(const Kokkos::View<T*, ExecSpace>& v) {
  T last_element;
  Kokkos::deep_copy(last_element, Kokkos::subview(v, v.extent(0) - 1));
  return last_element;
}

// Helper to allocate input and output views
template <class T>
auto prepare_views(const std::size_t kProbSize) {
  Kokkos::View<T*, ExecSpace> in{"input", kProbSize};
  Kokkos::View<T*, ExecSpace> out{"output", kProbSize};

  auto h_in = Kokkos::create_mirror_view(in);

  for (std::size_t i = 0; i < kProbSize; ++i) {
    h_in(i) = i;
  }

  Kokkos::deep_copy(in, h_in);

  return std::make_tuple(in, out, h_in);
}

// Perform scan with a reference implementation
template <class T, class ViewType, class ScanFunctor = SumFunctor<T>>
T ref_scan(const ViewType& h_in, ScanFunctor scan_functor = ScanFunctor()) {
  std::size_t view_size = h_in.extent(0);

  Kokkos::View<T*, HostExecSpace> h_out("output", view_size);

  // FIXME: We have GCC 8.4.0 based check in our ORNL Jenkins CI.
  // std::inclusive_scan is available only from GCC 9.3. Since, GCC 9.1
  // std::inclusive_scan that takes execution policy is available. However,
  // there is error with <execution> header before GCC 10.1.
  h_out(0) = h_in(0);

  for (std::size_t i = 1; i < view_size; ++i) {
    h_out(i) = scan_functor(h_in(i), h_out(i - 1));
  }

  return h_out(view_size - 1);
}

// Inclusive Scan with default binary operation (sum) or user provided functor
// Note: The nature of the functor must be compatible with the
// elements in the input and output views
template <class T, template <class> class ScanFunctor = ImpSumBinOp>
auto inclusive_scan(const Kokkos::View<T*, ExecSpace>& in,
                    const Kokkos::View<T*, ExecSpace>& out, T res_check) {
  ExecSpace().fence();
  Kokkos::Timer timer;

  if constexpr (std::is_same_v<ScanFunctor<T>, ImpSumBinOp<T>>) {
    KE::inclusive_scan("Default scan", ExecSpace(), KE::cbegin(in),
                       KE::cend(in), KE::begin(out));
  } else {
    KE::inclusive_scan("Scan using a functor", ExecSpace(), KE::cbegin(in),
                       KE::cend(in), KE::begin(out), ScanFunctor<T>());
  }

  ExecSpace().fence();
  double time_scan = timer.seconds();

  T res_scan  = obtain_last_elem(out);
  bool passed = (res_check == res_scan);

  return std::make_tuple(time_scan, passed);
}

// Benchmark: Inclusive Scan with default binary operation (sum)
// or user provided functor
template <class T, template <class> class ScanFunctor = ImpSumBinOp>
void BM_inclusive_scan(benchmark::State& state) {
  const std::size_t kProbSize = state.range(0);

  auto [in, out, h_in] = prepare_views<T>(kProbSize);

  T res_check;

  if constexpr (std::is_same_v<ScanFunctor<T>, ImpSumBinOp<T>>) {
    res_check = ref_scan<T>(h_in);
  } else {
    res_check = ref_scan<T>(h_in, ScanFunctor<T>());
  }

  double time_scan = 0.;
  bool passed      = false;

  for (auto _ : state) {
    if constexpr (std::is_same_v<ScanFunctor<T>, ImpSumBinOp<T>>) {
      std::tie(time_scan, passed) = inclusive_scan<T>(in, out, res_check);
    } else {
      std::tie(time_scan, passed) =
          inclusive_scan<T, ScanFunctor>(in, out, res_check);
    }

    KokkosBenchmark::report_results(state, in, 2, time_scan);
    state.counters["Passed"] = passed;
  }
}

constexpr std::size_t PROB_SIZE = 100'000'000;

}  // anonymous namespace

// FIXME: Add logic to pass min. warm-up time. Also, the value should be set
// by the user. Say, via the environment variable BENCHMARK_MIN_WARMUP_TIME.

BENCHMARK(BM_inclusive_scan<std::uint64_t>)->Arg(PROB_SIZE)->UseManualTime();
BENCHMARK(BM_inclusive_scan<std::int64_t>)->Arg(PROB_SIZE)->UseManualTime();
BENCHMARK(BM_inclusive_scan<double>)->Arg(PROB_SIZE)->UseManualTime();
BENCHMARK(BM_inclusive_scan<std::uint64_t, SumFunctor>)
    ->Arg(PROB_SIZE)
    ->UseManualTime();
BENCHMARK(BM_inclusive_scan<std::int64_t, SumFunctor>)
    ->Arg(PROB_SIZE)
    ->UseManualTime();
BENCHMARK(BM_inclusive_scan<double, SumFunctor>)
    ->Arg(PROB_SIZE)
    ->UseManualTime();
BENCHMARK(BM_inclusive_scan<std::uint64_t, MaxFunctor>)
    ->Arg(PROB_SIZE)
    ->UseManualTime();
BENCHMARK(BM_inclusive_scan<std::int64_t, MaxFunctor>)
    ->Arg(PROB_SIZE)
    ->UseManualTime();
BENCHMARK(BM_inclusive_scan<double, MaxFunctor>)
    ->Arg(PROB_SIZE)
    ->UseManualTime();
