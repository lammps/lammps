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

// export OMP_PROC_BIND=spread ; export OMP_PLACES=threads
// c++  -O2 -g -DNDEBUG  -fopenmp
// ../core/perf_test/test_atomic_minmax_simple.cpp -I../core/src/ -I. -o
// test_atomic_minmax_simple.x  containers/src/libkokkoscontainers.a
// core/src/libkokkoscore.a -ldl && OMP_NUM_THREADS=1
// ./test_atomic_minmax_simple.x 10000000

#include <benchmark/benchmark.h>

#include "Benchmark_Context.hpp"
#include "PerfTest_Category.hpp"
#include <Kokkos_Core.hpp>

using exec_space = Kokkos::DefaultExecutionSpace;

constexpr int LENGTH = 1'000'000;

template <typename T>
Kokkos::View<T*, exec_space> prepare_input(const int length, const T value) {
  Kokkos::View<T*, exec_space> input("input", length);
  Kokkos::parallel_for(
      length, KOKKOS_LAMBDA(const int i) { input(i) = value; });
  Kokkos::fence();
  return input;
}

int get_length(benchmark::State& state) {
  return (Test::command_line_num_args() == 2)
             ? std::stoi(Test::command_line_arg(1))
             : state.range(0);
}

template <typename T>
int check_errors_replacement(Kokkos::View<T*, exec_space> view) {
  int errors = 0;
  Kokkos::parallel_reduce(
      view.size(),
      KOKKOS_LAMBDA(const int i, int& inner) { inner += (view(i) != (T)i); },
      errors);
  Kokkos::fence();
  return errors;
}

template <typename T>
double atomic_min_replacement(Kokkos::View<T*, exec_space> input) {
  const int length = input.size();
  Kokkos::Timer timer;
  Kokkos::parallel_for(
      length, KOKKOS_LAMBDA(const int i) {
        (void)Kokkos::atomic_fetch_min(&(input(i)), (T)i);
      });
  Kokkos::fence();
  return timer.seconds();
}

template <typename T>
static void Atomic_MinReplacements(benchmark::State& state) {
  const int length = get_length(state);
  auto inp         = prepare_input(length, std::numeric_limits<T>::max());

  for (auto _ : state) {
    const auto time   = atomic_min_replacement(inp);
    const auto errors = check_errors_replacement(inp);

    // report results
    state.SetIterationTime(time);
    if (errors > 0) {
      state.counters["Errors"] = benchmark::Counter(errors);
    }
  }
}

template <typename T>
double atomic_max_replacement(Kokkos::View<T*, exec_space> input) {
  const int length = input.size();
  Kokkos::Timer timer;
  Kokkos::parallel_for(
      length, KOKKOS_LAMBDA(const int i) {
        (void)Kokkos::atomic_max_fetch(&(input(i)), (T)i);
      });
  Kokkos::fence();
  return timer.seconds();
}

template <typename T>
static void Atomic_MaxReplacements(benchmark::State& state) {
  const auto length = get_length(state);
  auto inp          = prepare_input(length, std::numeric_limits<T>::lowest());

  for (auto _ : state) {
    const auto time   = atomic_max_replacement(inp);
    const auto errors = check_errors_replacement(inp);

    // report results
    state.SetIterationTime(time);
    if (errors > 0) {
      state.counters["Errors"] = benchmark::Counter(errors);
    }
  }
}

template <typename T>
int check_errors_early_exit(Kokkos::View<T*, exec_space> view, const T ref) {
  int errors = 0;
  Kokkos::parallel_reduce(
      view.size(),
      KOKKOS_LAMBDA(const int i, int& inner) { inner += (view(i) != ref); },
      errors);
  Kokkos::fence();
  return errors;
}

template <typename T>
static void Atomic_MaxEarlyExits(benchmark::State& state) {
  const auto length = get_length(state);
  auto inp          = prepare_input(length, std::numeric_limits<T>::max());

  for (auto _ : state) {
    const auto time = atomic_max_replacement(inp);
    const auto errors =
        check_errors_early_exit(inp, std::numeric_limits<T>::max());

    // report results
    state.SetIterationTime(time);
    if (errors > 0) {
      state.counters["Errors"] = benchmark::Counter(errors);
    }
  }
}

template <typename T>
static void Atomic_MinEarlyExits(benchmark::State& state) {
  const auto length = get_length(state);
  auto inp          = prepare_input(length, std::numeric_limits<T>::lowest());

  for (auto _ : state) {
    const auto time = atomic_min_replacement(inp);
    const auto errors =
        check_errors_early_exit(inp, std::numeric_limits<T>::lowest());

    // report results
    state.SetIterationTime(time);
    if (errors > 0) {
      state.counters["Errors"] = benchmark::Counter(errors);
    }
  }
}

template <typename T>
void report_errors_contentious_replacement(benchmark::State& state,
                                           const T final, const T first,
                                           const T expected) {
  state.counters["Errors"]   = benchmark::Counter(1);
  state.counters["Final"]    = benchmark::Counter(final);
  state.counters["First"]    = benchmark::Counter(first);
  state.counters["Expected"] = benchmark::Counter(expected);
}

template <typename T>
double atomic_contentious_max_replacement(benchmark::State& state,
                                          Kokkos::View<T*, exec_space> input,
                                          const int con_length) {
  const auto max = std::numeric_limits<T>::max();
  T current      = 0;

  Kokkos::Timer timer;
  Kokkos::parallel_reduce(
      con_length,
      KOKKOS_LAMBDA(const int i, T& inner) {
        inner = Kokkos::atomic_max_fetch(&(input(0)), inner + 1);
        if (i == con_length - 1) {
          Kokkos::atomic_max_fetch(&(input(0)), max);
          inner = max;
        }
      },
      Kokkos::Max<T>(current));
  Kokkos::fence();
  const auto time = timer.seconds();

  if (current < max) {
    report_errors_contentious_replacement(state, current, input(0), max);
  }

  return time;
}

template <typename T>
static void Atomic_ContentiousMaxReplacements(benchmark::State& state) {
  const auto length = get_length(state);
  auto inp          = prepare_input(1, std::numeric_limits<T>::lowest());

  for (auto _ : state) {
    const auto time = atomic_contentious_max_replacement(state, inp, length);

    state.SetIterationTime(time);
  }
}

template <typename T>
double atomic_contentious_min_replacement(benchmark::State& state,
                                          Kokkos::View<T*, exec_space> input,
                                          const int con_length) {
  const auto min = std::numeric_limits<T>::lowest();
  T current      = 0;

  Kokkos::Timer timer;
  Kokkos::parallel_reduce(
      con_length,
      KOKKOS_LAMBDA(const int i, T& inner) {
        inner = Kokkos::atomic_min_fetch(&(input(0)), inner - 1);
        if (i == con_length - 1) {
          Kokkos::atomic_min_fetch(&(input(0)), min);
          inner = min;
        }
      },
      Kokkos::Min<T>(current));
  Kokkos::fence();
  const auto time = timer.seconds();

  if (current > min) {
    report_errors_contentious_replacement(state, current, input(0), min);
  }

  return time;
}

template <typename T>
static void Atomic_ContentiousMinReplacements(benchmark::State& state) {
  const auto length = get_length(state);
  auto inp          = prepare_input(1, std::numeric_limits<T>::max());

  for (auto _ : state) {
    const auto time = atomic_contentious_max_replacement(state, inp, length);

    state.SetIterationTime(time);
  }
}

// int
BENCHMARK(Atomic_MinReplacements<int>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxReplacements<int>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxEarlyExits<int>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MinEarlyExits<int>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMaxReplacements<int>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMinReplacements<int>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);
///////////////////////////////////////////////////////////////////////

// long
BENCHMARK(Atomic_MinReplacements<long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxReplacements<long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxEarlyExits<long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MinEarlyExits<long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMaxReplacements<long>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMinReplacements<long>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);
///////////////////////////////////////////////////////////////////////

// long long
BENCHMARK(Atomic_MinReplacements<long long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxReplacements<long long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxEarlyExits<long long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MinEarlyExits<long long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMaxReplacements<long long>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMinReplacements<long long>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);
///////////////////////////////////////////////////////////////////////

// unsigned int
BENCHMARK(Atomic_MinReplacements<unsigned int>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxReplacements<unsigned int>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxEarlyExits<unsigned int>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MinEarlyExits<unsigned int>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMaxReplacements<unsigned int>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMinReplacements<unsigned int>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);
///////////////////////////////////////////////////////////////////////

// unsigned long
BENCHMARK(Atomic_MinReplacements<unsigned long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxReplacements<unsigned long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxEarlyExits<unsigned long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MinEarlyExits<unsigned long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMaxReplacements<unsigned long>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMinReplacements<unsigned long>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);
///////////////////////////////////////////////////////////////////////

// unsigned long long
BENCHMARK(Atomic_MinReplacements<unsigned long long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxReplacements<unsigned long long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxEarlyExits<unsigned long long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MinEarlyExits<unsigned long long>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMaxReplacements<unsigned long long>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMinReplacements<unsigned long long>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);
///////////////////////////////////////////////////////////////////////

// float
BENCHMARK(Atomic_MinReplacements<float>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxReplacements<float>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxEarlyExits<float>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MinEarlyExits<float>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMaxReplacements<float>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMinReplacements<float>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);
///////////////////////////////////////////////////////////////////////

// double
BENCHMARK(Atomic_MinReplacements<double>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxReplacements<double>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MaxEarlyExits<double>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_MinEarlyExits<double>)
    ->ArgName("Length")
    ->Arg(LENGTH)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMaxReplacements<double>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);

BENCHMARK(Atomic_ContentiousMinReplacements<double>)
    ->ArgName("Length")
    ->Arg(LENGTH / 5)
    ->UseManualTime()
    ->Iterations(10);
