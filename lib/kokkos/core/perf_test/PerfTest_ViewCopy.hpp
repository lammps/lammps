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

#ifndef KOKKOS_CORE_PERFTEST_BENCHMARK_VIEW_COPY_HPP
#define KOKKOS_CORE_PERFTEST_BENCHMARK_VIEW_COPY_HPP

#include <Kokkos_Core.hpp>

#include <benchmark/benchmark.h>

#include "Benchmark_Context.hpp"
#include <cmath>

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

template <class LayoutA, class LayoutB>
static void ViewDeepCopy_Rank1(benchmark::State& state) {
  const int N8 = std::pow(state.range(0), 8);

  Kokkos::View<double*, LayoutA> a("A1", N8);
  Kokkos::View<double*, LayoutB> b("B1", N8);

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

}  // namespace Test

#endif
