// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/**
 * @file PerfTest_Stream.cpp
 * @brief Implementation of STREAM benchmark operations for Kokkos.
 *
 * @details This file provides a set of memory bandwidth benchmarks based on the
 * STREAM benchmark suite. It implements the five core STREAM operations (Set,
 * Copy, Scale, Add, and Triad) using Kokkos parallel primitives. It includes
 * validation.
 *
 * The implementation strives to use as few Kokkos features as possible, thus
 * validation is performed on the host rather than via parallel_reduce.
 */

#include <Kokkos_Core.hpp>
#include <benchmark/benchmark.h>
#include "Benchmark_Context.hpp"

namespace {

using StreamType                   = double;
constexpr static StreamType A_INIT = 1.0;
constexpr static StreamType B_INIT = 2.0;
constexpr static StreamType C_INIT = 3.0;
constexpr static StreamType SCALAR = 4.0;

template <unsigned MemTraits>
using StreamView = Kokkos::View<StreamType*, Kokkos::MemoryTraits<MemTraits>>;

// different than benchmarks/stream, which uses int
// wide index types are common as GPU memory grows
using StreamIndex = int64_t;
using Policy      = Kokkos::RangePolicy<Kokkos::IndexType<StreamIndex>>;

template <typename V>
void perform_set(const V& a, typename V::const_value_type scalar) {
  Kokkos::parallel_for(
      "set", Policy(0, a.extent(0)),
      KOKKOS_LAMBDA(const StreamIndex i) { a[i] = scalar; });

  Kokkos::fence();
}

template <typename V>
void perform_copy(const V& a, const V& b) {
  Kokkos::parallel_for(
      "copy", Policy(0, a.extent(0)),
      KOKKOS_LAMBDA(const StreamIndex i) { b[i] = a[i]; });

  Kokkos::fence();
}

template <typename V>
void perform_scale(const V& b, const V& c,
                   typename V::const_value_type scalar) {
  Kokkos::parallel_for(
      "scale", Policy(0, b.extent(0)),
      KOKKOS_LAMBDA(const StreamIndex i) { b[i] = scalar * c[i]; });

  Kokkos::fence();
}

template <typename V>
void perform_add(const V& a, const V& b, const V& c) {
  Kokkos::parallel_for(
      "add", Policy(0, a.extent(0)),
      KOKKOS_LAMBDA(const StreamIndex i) { c[i] = a[i] + b[i]; });

  Kokkos::fence();
}

template <typename V>
void perform_triad(const V& a, const V& b, const V& c,
                   typename V::const_value_type scalar) {
  Kokkos::parallel_for(
      "triad", Policy(0, a.extent(0)),
      KOKKOS_LAMBDA(const StreamIndex i) { a[i] = b[i] + scalar * c[i]; });

  Kokkos::fence();
}

template <typename V>
int validate_array(V& a_dev, typename V::const_value_type expected) {
  using scalar_type = typename V::non_const_value_type;

  const auto a =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, a_dev);

  scalar_type error = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    error += std::abs(a[i] - expected);
  }
  const scalar_type avgError = error / (scalar_type)a.size();

  // all values here are pretty easy for float types to represent
  // so let's make the tolerances very tight.
  return std::abs(avgError / expected) >
         Kokkos::Experimental::epsilon_v<scalar_type>;
}

template <unsigned MemTraits>
static void StreamSet(benchmark::State& state) {
  const size_t N8                 = std::pow(state.range(0), 8);
  static constexpr int DATA_RATIO = 1;

  StreamView<MemTraits> a(Kokkos::view_alloc(Kokkos::WithoutInitializing, "a"),
                          N8);

  for (auto _ : state) {
    Kokkos::Timer timer;
    perform_set(a, SCALAR);
    KokkosBenchmark::report_results(state, a, DATA_RATIO, timer.seconds());
  }

  if (validate_array(a, SCALAR)) {
    state.SkipWithError("validation failure");
  }
}

template <unsigned MemTraits>
static void StreamCopy(benchmark::State& state) {
  const size_t N8                 = std::pow(state.range(0), 8);
  static constexpr int DATA_RATIO = 2;

  StreamView<MemTraits> a("a", N8), b("b", N8);

  perform_set(a, A_INIT);

  for (auto _ : state) {
    Kokkos::Timer timer;
    perform_copy(a, b);
    KokkosBenchmark::report_results(state, a, DATA_RATIO, timer.seconds());
  }

  if (validate_array(b, A_INIT)) {
    state.SkipWithError("validation failure");
  }
}

template <unsigned MemTraits>
static void StreamScale(benchmark::State& state) {
  const size_t N8                 = std::pow(state.range(0), 8);
  static constexpr int DATA_RATIO = 2;

  StreamView<MemTraits> a("a", N8), b("b", N8);

  perform_set(b, B_INIT);

  for (auto _ : state) {
    Kokkos::Timer timer;
    perform_scale(a, b, SCALAR);
    KokkosBenchmark::report_results(state, b, DATA_RATIO, timer.seconds());
  }

  if (validate_array(a, B_INIT * SCALAR)) {
    state.SkipWithError("validation failure");
  }
}

template <unsigned MemTraits>
static void StreamAdd(benchmark::State& state) {
  const size_t N8                 = std::pow(state.range(0), 8);
  static constexpr int DATA_RATIO = 3;

  StreamView<MemTraits> a("a", N8), b("b", N8), c("c", N8);

  perform_set(a, A_INIT);
  perform_set(b, B_INIT);
  perform_set(c, C_INIT);

  for (auto _ : state) {
    Kokkos::Timer timer;
    perform_add(a, b, c);
    KokkosBenchmark::report_results(state, c, DATA_RATIO, timer.seconds());
  }

  if (validate_array(c, A_INIT + B_INIT)) {
    state.SkipWithError("validation failure");
  }
}

template <unsigned MemTraits>
static void StreamTriad(benchmark::State& state) {
  const size_t N8                 = std::pow(state.range(0), 8);
  static constexpr int DATA_RATIO = 3;

  StreamView<MemTraits> a("a", N8), b("b", N8), c("c", N8);

  perform_set(a, A_INIT);
  perform_set(b, B_INIT);
  perform_set(c, C_INIT);

  for (auto _ : state) {
    Kokkos::Timer timer;
    perform_triad(a, b, c, SCALAR);
    KokkosBenchmark::report_results(state, a, DATA_RATIO, timer.seconds());
  }

  if (validate_array(a, B_INIT + SCALAR * C_INIT)) {
    state.SkipWithError("validation failure");
  }
}

// skips a benchmark with an error from thrown exceptions
template <void (*bm)(benchmark::State&)>
static void or_skip(benchmark::State& state) {
  try {
    bm(state);
  } catch (const std::runtime_error& e) {
    state.SkipWithError(e.what());
  }
}

// As of May 2025, 10^8 doubles is larger than caches, but not so large as
// to be inconvenient. Also run 11^8 for a quick check of convergence.
#define STREAM_ARGS(label)            \
  Name(label)                         \
      ->ArgName("N")                  \
      ->Arg(10)                       \
      ->Arg(11)                       \
      ->Unit(benchmark::kMillisecond) \
      ->UseManualTime()

// clang-format off
// clang-format formatted these lines inconsistently, making it hard to
// see the common pattern
BENCHMARK(or_skip<StreamSet<0>>)
    ->STREAM_ARGS("StreamSet");

BENCHMARK(or_skip<StreamSet<Kokkos::Restrict>>)
    ->STREAM_ARGS("StreamSet<Restrict>");

BENCHMARK(or_skip<StreamCopy<0>>)
    ->STREAM_ARGS("StreamCopy");

BENCHMARK(or_skip<StreamCopy<Kokkos::Restrict>>)
    ->STREAM_ARGS("StreamCopy<Restrict>");

BENCHMARK(or_skip<StreamScale<0>>)
    ->STREAM_ARGS("StreamScale");

BENCHMARK(or_skip<StreamScale<Kokkos::Restrict>>)
    ->STREAM_ARGS("StreamScale<Restrict>");

BENCHMARK(or_skip<StreamAdd<0>>)
    ->STREAM_ARGS("StreamAdd");

BENCHMARK(or_skip<StreamAdd<Kokkos::Restrict>>)
    ->STREAM_ARGS("StreamAdd<Restrict>");

BENCHMARK(or_skip<StreamTriad<0>>)
    ->STREAM_ARGS("StreamTriad");

BENCHMARK(or_skip<StreamTriad<Kokkos::Restrict>>)
    ->STREAM_ARGS("StreamTriad<Restrict>");
// clang-format on

#undef STREAM_ARGS

}  // namespace
