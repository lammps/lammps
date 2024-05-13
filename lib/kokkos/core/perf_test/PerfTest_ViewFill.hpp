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

#include "Benchmark_Context.hpp"

#include <cmath>

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

template <class Layout>
static void ViewFill_Rank1(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;
  const int N4 = N2 * N2;
  const int N8 = N4 * N4;

  Kokkos::View<double*, Layout> a("A1", N8);
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

}  // namespace Test
