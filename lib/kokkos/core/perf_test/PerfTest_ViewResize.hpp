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

#include <Kokkos_Core.hpp>
#include <benchmark/benchmark.h>
#include <cmath>
#include "Benchmark_Context.hpp"

namespace Test {

static constexpr int R = 10;
static constexpr int N = 10;

template <class Layout>
static void ViewResize_Rank1(benchmark::State& state) {
  const int N8 = std::pow(state.range(0), 8);
  Kokkos::View<double*, Layout> a("A1", N8);
  Kokkos::View<double*, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(a_, int(N8 * 1.1));
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_Rank2(benchmark::State& state) {
  const int N4 = std::pow(state.range(0), 4);
  Kokkos::View<double**, Layout> a("A2", N4, N4);
  Kokkos::View<double**, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(a_, int(N4 * 1.1), N4);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_Rank3(benchmark::State& state) {
  const int N2 = std::pow(state.range(0), 2);
  const int N3 = std::pow(state.range(0), 3);
  Kokkos::View<double***, Layout> a("A3", N3, N3, N2);
  Kokkos::View<double***, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(a_, int(N3 * 1.1), N3, N2);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_Rank4(benchmark::State& state) {
  const int N2 = std::pow(state.range(0), 2);
  Kokkos::View<double****, Layout> a("A4", N2, N2, N2, N2);
  Kokkos::View<double****, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(a_, int(N2 * 1.1), N2, N2, N2);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_Rank5(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double*****, Layout> a("A5", N2, N2, N1, N1, N2);
  Kokkos::View<double*****, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(a_, int(N2 * 1.1), N2, N1, N1, N2);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_Rank6(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double******, Layout> a("A6", N2, N1, N1, N1, N1, N2);
  Kokkos::View<double******, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(a_, int(N2 * 1.1), N1, N1, N1, N1, N2);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_Rank7(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double*******, Layout> a("A7", N2, N1, N1, N1, N1, N1, N1);
  Kokkos::View<double*******, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(a_, int(N2 * 1.1), N1, N1, N1, N1, N1, N1);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_Rank8(benchmark::State& state) {
  const int N1 = state.range(0);

  Kokkos::View<double********, Layout> a("A8", N1, N1, N1, N1, N1, N1, N1, N1);
  Kokkos::View<double********, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(a_, int(N1 * 1.1), N1, N1, N1, N1, N1, N1, N1);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_NoInit_Rank1(benchmark::State& state) {
  const int N8 = std::pow(state.range(0), 8);
  Kokkos::View<double*, Layout> a("A1", N8);
  Kokkos::View<double*, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(Kokkos::WithoutInitializing, a_, int(N8 * 1.1));
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_NoInit_Rank2(benchmark::State& state) {
  const int N4 = std::pow(state.range(0), 4);
  Kokkos::View<double**, Layout> a("A2", N4, N4);
  Kokkos::View<double**, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(Kokkos::WithoutInitializing, a_, int(N4 * 1.1), N4);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_NoInit_Rank3(benchmark::State& state) {
  const int N2 = std::pow(state.range(0), 2);
  const int N3 = std::pow(state.range(0), 3);
  Kokkos::View<double***, Layout> a("A3", N3, N3, N2);
  Kokkos::View<double***, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(Kokkos::WithoutInitializing, a_, int(N3 * 1.1), N3, N2);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_NoInit_Rank4(benchmark::State& state) {
  const int N2 = std::pow(state.range(0), 2);
  Kokkos::View<double****, Layout> a("A4", N2, N2, N2, N2);
  Kokkos::View<double****, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(Kokkos::WithoutInitializing, a_, int(N2 * 1.1), N2, N2, N2);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_NoInit_Rank5(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double*****, Layout> a("A5", N2, N2, N1, N1, N2);
  Kokkos::View<double*****, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(Kokkos::WithoutInitializing, a_, int(N2 * 1.1), N2, N1, N1,
                   N2);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_NoInit_Rank6(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double******, Layout> a("A6", N2, N1, N1, N1, N1, N2);
  Kokkos::View<double******, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(Kokkos::WithoutInitializing, a_, int(N2 * 1.1), N1, N1, N1,
                   N1, N2);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_NoInit_Rank7(benchmark::State& state) {
  const int N1 = state.range(0);
  const int N2 = N1 * N1;

  Kokkos::View<double*******, Layout> a("A7", N2, N1, N1, N1, N1, N1, N1);
  Kokkos::View<double*******, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(Kokkos::WithoutInitializing, a_, int(N2 * 1.1), N1, N1, N1,
                   N1, N1, N1);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_NoInit_Rank8(benchmark::State& state) {
  const int N1 = state.range(0);

  Kokkos::View<double********, Layout> a("A8", N1, N1, N1, N1, N1, N1, N1, N1);
  Kokkos::View<double********, Layout> a_(a);

  for (auto _ : state) {
    Kokkos::fence();
    Kokkos::Timer timer;
    Kokkos::resize(Kokkos::WithoutInitializing, a_, int(N1 * 1.1), N1, N1, N1,
                   N1, N1, N1, N1);
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

template <class Layout>
static void ViewResize_NoInit_Raw(benchmark::State& state) {
  const int N8 = std::pow(state.range(0), 8);
  Kokkos::View<double*, Layout> a("A1", N8);
  double* a_ptr = a.data();

  for (auto _ : state) {
    Kokkos::Timer timer;
    Kokkos::View<double*, Layout> a1(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "A1"), int(N8 * 1.1));
    double* a1_ptr = a1.data();
    Kokkos::parallel_for(
        N8, KOKKOS_LAMBDA(const int& i) { a1_ptr[i] = a_ptr[i]; });
    Kokkos::fence();
    KokkosBenchmark::report_results(state, a, 2, timer.seconds());
  }
}

}  // namespace Test
