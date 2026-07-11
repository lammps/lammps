// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "Benchmark_Context.hpp"

#include <cmath>
#include <stdexcept>

namespace Test {

static constexpr int N = 10;

template <class ViewType>
void fill_view(ViewType& a, typename ViewType::const_value_type& val,
               benchmark::State& state) {
  for (auto _ : state) {
    Kokkos::fence();

    Kokkos::Timer timer;
    Kokkos::deep_copy(a, val);
    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

template <
    class Layout,
    class MemorySpace = typename Kokkos::DefaultExecutionSpace::memory_space>
static void ViewFill_Rank1(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;
  const int N4 = N2 * N2;
  const int N8 = N4 * N4;

  Kokkos::View<double*, Layout, MemorySpace> a("A1", N8);
  fill_view(a, 1.1, state);
}

template <class Layout>
static void ViewFill_Rank2(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;
  const int N4 = N2 * N2;

  Kokkos::View<double**, Layout> a("A2", N4, N4);
  fill_view(a, 1.1, state);
}

template <class Layout>
static void ViewFill_Rank3(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;
  const int N3 = N2 * N1;

  Kokkos::View<double***, Layout> a("A3", N3, N3, N2);
  fill_view(a, 1.1, state);
}

template <class Layout>
static void ViewFill_Rank4(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double****, Layout> a("A4", N2, N2, N2, N2);
  fill_view(a, 1.1, state);
}

template <class Layout>
static void ViewFill_Rank5(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double*****, Layout> a("A5", N2, N2, N1, N1, N2);
  fill_view(a, 1.1, state);
}

template <class Layout>
static void ViewFill_Rank6(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double******, Layout> a("A6", N2, N1, N1, N1, N1, N2);
  fill_view(a, 1.1, state);
}

template <class Layout>
static void ViewFill_Rank7(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double*******, Layout> a("A7", N2, N1, N1, N1, N1, N1, N1);
  fill_view(a, 1.1, state);
}

template <class Layout>
static void ViewFill_Rank8(benchmark::State& state) {
  const int N1 = state.range(0);

  Kokkos::View<double********, Layout> a("A8", N1, N1, N1, N1, N1, N1, N1, N1);
  fill_view(a, 1.1, state);
}

template <class Layout>
static void ViewFill_Raw(benchmark::State& state) {
  const int N8 = std::pow(state.range(0), 8);

  Kokkos::View<double*, Layout> a("A1", N8);
  double* a_ptr = a.data();

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::parallel_for(
        N8, KOKKOS_LAMBDA(const int& i) { a_ptr[i] = 1.1; });
    Kokkos::fence();

    KokkosBenchmark::report_results(state, a, 1, timer.seconds());
  }
}

[[maybe_unused]] static void ViewFill_Rank1Strided(benchmark::State& state) {
  const size_t N8 = std::pow(state.range(0), 8);

  // This benchmark allocates more data in order to measure a view fill
  // of the same size as the contiguous benchmarks, so in cases where they
  // can be run, this one may fail to allocate data (e.g., on a small CI runner)
  try {
    // allocate 2x the size since layout only has 1/2 the elements
    Kokkos::View<double*> a("A1", N8 * 2);

    Kokkos::LayoutStride layout(N8 / 2, 2);
    Kokkos::View<double*, Kokkos::LayoutStride> a_stride(a.data(), layout);

    fill_view(a_stride, 1.1, state);

  } catch (const std::runtime_error& e) {
    state.SkipWithError(e.what());
  }
}

}  // namespace Test
