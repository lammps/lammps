// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CORE_PERFTEST_BENCHMARK_VIEW_COPY_HPP
#define KOKKOS_CORE_PERFTEST_BENCHMARK_VIEW_COPY_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include <benchmark/benchmark.h>

#include "Benchmark_Context.hpp"
#include <cmath>
#include <stdexcept>

namespace Test {

static constexpr int DATA_RATIO = 2;

template <class ViewTypeA, class ViewTypeB>
void deepcopy_view(ViewTypeA& a, ViewTypeB& b, benchmark::State& state) {
  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::deep_copy(a, b);
    KokkosBenchmark::report_results(state, a, DATA_RATIO, timer.seconds());
  }
}

template <
    class LayoutA, class LayoutB,
    class MemorySpaceA = typename Kokkos::DefaultExecutionSpace::memory_space,
    class MemorySpaceB = typename Kokkos::DefaultExecutionSpace::memory_space>
static void ViewDeepCopy_Rank1(benchmark::State& state) {
  const int N8 = std::pow(state.range(0), 8);

  Kokkos::View<double*, LayoutA, MemorySpaceA> a("A1", N8);
  Kokkos::View<double*, LayoutB, MemorySpaceB> b("B1", N8);

  deepcopy_view(a, b, state);
}

template <class LayoutA, class LayoutB>
static void ViewDeepCopy_Rank2(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;
  const int N4 = N2 * N2;

  Kokkos::View<double**, LayoutA> a("A2", N4, N4);
  Kokkos::View<double**, LayoutB> b("B2", N4, N4);

  deepcopy_view(a, b, state);
}

template <class LayoutA, class LayoutB>
static void ViewDeepCopy_Rank3(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;
  const int N3 = N2 * N1;

  Kokkos::View<double***, LayoutA> a("A3", N3, N3, N2);
  Kokkos::View<double***, LayoutB> b("B3", N3, N3, N2);

  deepcopy_view(a, b, state);
}

template <class LayoutA, class LayoutB>
static void ViewDeepCopy_Rank4(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double****, LayoutA> a("A4", N2, N2, N2, N2);
  Kokkos::View<double****, LayoutB> b("B4", N2, N2, N2, N2);

  deepcopy_view(a, b, state);
}

template <class LayoutA, class LayoutB>
static void ViewDeepCopy_Rank5(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double*****, LayoutA> a("A5", N2, N2, N1, N1, N2);
  Kokkos::View<double*****, LayoutB> b("B5", N2, N2, N1, N1, N2);

  deepcopy_view(a, b, state);
}

template <class LayoutA, class LayoutB>
static void ViewDeepCopy_Rank6(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double******, LayoutA> a("A6", N2, N1, N1, N1, N1, N2);
  Kokkos::View<double******, LayoutB> b("B6", N2, N1, N1, N1, N1, N2);

  deepcopy_view(a, b, state);
}

template <class LayoutA, class LayoutB>
static void ViewDeepCopy_Rank7(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double*******, LayoutA> a("A7", N2, N1, N1, N1, N1, N1, N1);
  Kokkos::View<double*******, LayoutB> b("B7", N2, N1, N1, N1, N1, N1, N1);

  deepcopy_view(a, b, state);
}

template <class LayoutA, class LayoutB>
static void ViewDeepCopy_Rank8(benchmark::State& state) {
  const int N1 = state.range(0);

  Kokkos::View<double********, LayoutA> a("A8", N1, N1, N1, N1, N1, N1, N1, N1);
  Kokkos::View<double********, LayoutB> b("B8", N1, N1, N1, N1, N1, N1, N1, N1);

  deepcopy_view(a, b, state);
}

template <class LayoutA, class LayoutB>
static void ViewDeepCopy_Raw(benchmark::State& state) {
  const int N8 = std::pow(state.range(0), 8);

  Kokkos::View<double*, LayoutA> a("A1", N8);
  Kokkos::View<double*, LayoutB> b("B1", N8);
  double* const a_ptr       = a.data();
  const double* const b_ptr = b.data();

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::parallel_for(
        N8, KOKKOS_LAMBDA(const int& i) { a_ptr[i] = b_ptr[i]; });
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, DATA_RATIO, timer.seconds());
  }
}

template <typename DstMemorySpace, typename SrcMemorySpace>
static void ViewDeepCopy_Rank1Strided(benchmark::State& state) {
  const size_t N8 = std::pow(state.range(0), 8);

  // This benchmark allocates more data in order to measure a deep_copy
  // of the same size as the contiguous benchmarks, so in cases where they
  // can be run, this one may fail to allocate data (e.g., on a small CI runner)
  try {
    // allocate 2x the size since layout only has 1/2 the elements
    Kokkos::View<double*, DstMemorySpace> a("A1", N8 * 2);
    Kokkos::View<double*, SrcMemorySpace> b("B1", N8 * 2);

    Kokkos::LayoutStride layout(N8 / 2, 2);
    Kokkos::View<double*, Kokkos::LayoutStride, DstMemorySpace> a_stride(
        a.data(), layout);
    Kokkos::View<double*, Kokkos::LayoutStride, SrcMemorySpace> b_stride(
        b.data(), layout);
    deepcopy_view(a_stride, b_stride, state);
  } catch (const std::runtime_error& e) {
    state.SkipWithError(e.what());
  }
}

}  // namespace Test

#endif
