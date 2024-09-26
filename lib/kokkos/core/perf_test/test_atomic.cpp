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

#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <benchmark/benchmark.h>
#include "Benchmark_Context.hpp"

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>

using exec_space = Kokkos::DefaultExecutionSpace;

template <class T, class DEVICE_TYPE>
struct ZeroFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = typename Kokkos::View<T, execution_space>;
  using h_type          = typename Kokkos::View<T, execution_space>::HostMirror;
  type data;
  KOKKOS_INLINE_FUNCTION
  void operator()(int) const { data() = 0; }
};

//---------------------------------------------------
//--------------atomic_fetch_add---------------------
//---------------------------------------------------

template <class T, class DEVICE_TYPE>
struct AddFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;
  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()(int) const { Kokkos::atomic_fetch_add(&data(), (T)1); }
};

template <class T>
T AddLoop(int loop) {
  struct ZeroFunctor<T, exec_space> f_zero;
  typename ZeroFunctor<T, exec_space>::type data("Data");
  typename ZeroFunctor<T, exec_space>::h_type h_data("HData");
  f_zero.data = data;
  Kokkos::parallel_for(1, f_zero);
  exec_space().fence();

  struct AddFunctor<T, exec_space> f_add;
  f_add.data = data;
  Kokkos::parallel_for(loop, f_add);
  exec_space().fence();

  Kokkos::deep_copy(h_data, data);
  T val = h_data();
  return val;
}

template <class T, class DEVICE_TYPE>
struct AddNonAtomicFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;
  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()(int) const { data() += (T)1; }
};

template <class T>
T AddLoopNonAtomic(int loop) {
  struct ZeroFunctor<T, exec_space> f_zero;
  typename ZeroFunctor<T, exec_space>::type data("Data");
  typename ZeroFunctor<T, exec_space>::h_type h_data("HData");

  f_zero.data = data;
  Kokkos::parallel_for(1, f_zero);
  exec_space().fence();

  struct AddNonAtomicFunctor<T, exec_space> f_add;
  f_add.data = data;
  Kokkos::parallel_for(loop, f_add);
  exec_space().fence();

  Kokkos::deep_copy(h_data, data);
  T val = h_data();

  return val;
}

template <class T>
T AddLoopSerial(int loop) {
  T* data = new T[1];
  data[0] = 0;

  for (int i = 0; i < loop; i++) *data += (T)1;

  T val = *data;
  delete[] data;
  return val;
}

template <class T, class DEVICE_TYPE>
struct CASFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;
  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()(int) const {
    T old = data();
    T newval, assumed;
    do {
      assumed = old;
      newval  = assumed + (T)1;
      old     = Kokkos::atomic_compare_exchange(&data(), assumed, newval);
    } while (old != assumed);
  }
};

template <class T>
T CASLoop(int loop) {
  struct ZeroFunctor<T, exec_space> f_zero;
  typename ZeroFunctor<T, exec_space>::type data("Data");
  typename ZeroFunctor<T, exec_space>::h_type h_data("HData");
  f_zero.data = data;
  Kokkos::parallel_for(1, f_zero);
  exec_space().fence();

  struct CASFunctor<T, exec_space> f_cas;
  f_cas.data = data;
  Kokkos::parallel_for(loop, f_cas);
  exec_space().fence();

  Kokkos::deep_copy(h_data, data);
  T val = h_data();

  return val;
}

template <class T, class DEVICE_TYPE>
struct CASNonAtomicFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;
  type data;

  KOKKOS_INLINE_FUNCTION
  void operator()(int) const {
    volatile T assumed;
    volatile T newval;
    bool fail = 1;
    do {
      assumed = data();
      newval  = assumed + (T)1;
      if (data() == assumed) {
        data() = newval;
        fail   = 0;
      }
    } while (fail);
  }
};

template <class T>
T CASLoopNonAtomic(int loop) {
  struct ZeroFunctor<T, exec_space> f_zero;
  typename ZeroFunctor<T, exec_space>::type data("Data");
  typename ZeroFunctor<T, exec_space>::h_type h_data("HData");
  f_zero.data = data;
  Kokkos::parallel_for(1, f_zero);
  exec_space().fence();

  struct CASNonAtomicFunctor<T, exec_space> f_cas;
  f_cas.data = data;
  Kokkos::parallel_for(loop, f_cas);
  exec_space().fence();

  Kokkos::deep_copy(h_data, data);
  T val = h_data();

  return val;
}

template <class T>
T CASLoopSerial(int loop) {
  T* data = new T[1];
  data[0] = 0;

  for (int i = 0; i < loop; i++) {
    T assumed;
    T newval;
    T old;
    do {
      assumed = *data;
      newval  = assumed + (T)1;
      old     = *data;
      *data   = newval;
    } while (!(assumed == old));
  }

  T val = *data;
  delete[] data;
  return val;
}

template <class T, class DEVICE_TYPE>
struct ExchFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;
  type data, data2;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    T old = Kokkos::atomic_exchange(&data(), (T)i);
    Kokkos::atomic_fetch_add(&data2(), old);
  }
};

template <class T>
T ExchLoop(int loop) {
  struct ZeroFunctor<T, exec_space> f_zero;
  typename ZeroFunctor<T, exec_space>::type data("Data");
  typename ZeroFunctor<T, exec_space>::h_type h_data("HData");
  f_zero.data = data;
  Kokkos::parallel_for(1, f_zero);
  exec_space().fence();

  typename ZeroFunctor<T, exec_space>::type data2("Data");
  typename ZeroFunctor<T, exec_space>::h_type h_data2("HData");
  f_zero.data = data2;
  Kokkos::parallel_for(1, f_zero);
  exec_space().fence();

  struct ExchFunctor<T, exec_space> f_exch;
  f_exch.data  = data;
  f_exch.data2 = data2;
  Kokkos::parallel_for(loop, f_exch);
  exec_space().fence();

  Kokkos::deep_copy(h_data, data);
  Kokkos::deep_copy(h_data2, data2);
  T val = h_data() + h_data2();

  return val;
}

template <class T, class DEVICE_TYPE>
struct ExchNonAtomicFunctor {
  using execution_space = DEVICE_TYPE;
  using type            = Kokkos::View<T, execution_space>;
  type data, data2;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    T old  = data();
    data() = (T)i;
    data2() += old;
  }
};

template <class T>
T ExchLoopNonAtomic(int loop) {
  struct ZeroFunctor<T, exec_space> f_zero;
  typename ZeroFunctor<T, exec_space>::type data("Data");
  typename ZeroFunctor<T, exec_space>::h_type h_data("HData");
  f_zero.data = data;
  Kokkos::parallel_for(1, f_zero);
  exec_space().fence();

  typename ZeroFunctor<T, exec_space>::type data2("Data");
  typename ZeroFunctor<T, exec_space>::h_type h_data2("HData");
  f_zero.data = data2;
  Kokkos::parallel_for(1, f_zero);
  exec_space().fence();

  struct ExchNonAtomicFunctor<T, exec_space> f_exch;
  f_exch.data  = data;
  f_exch.data2 = data2;
  Kokkos::parallel_for(loop, f_exch);
  exec_space().fence();

  Kokkos::deep_copy(h_data, data);
  Kokkos::deep_copy(h_data2, data2);
  T val = h_data() + h_data2();

  return val;
}

template <class T>
T ExchLoopSerial(int loop) {
  T* data  = new T[1];
  T* data2 = new T[1];
  data[0]  = 0;
  data2[0] = 0;
  for (int i = 0; i < loop; i++) {
    T old = *data;
    *data = (T)i;
    *data2 += old;
  }

  T val = *data2 + *data;
  delete[] data;
  delete[] data2;
  return val;
}

template <class T>
T LoopVariant(int loop, int test) {
  switch (test) {
    case 1: return AddLoop<T>(loop);
    case 2: return CASLoop<T>(loop);
    case 3: return ExchLoop<T>(loop);
  }
  return 0;
}

template <class T>
T LoopVariantSerial(int loop, int test) {
  switch (test) {
    case 1: return AddLoopSerial<T>(loop);
    case 2: return CASLoopSerial<T>(loop);
    case 3: return ExchLoopSerial<T>(loop);
  }
  return 0;
}

template <class T>
T LoopVariantNonAtomic(int loop, int test) {
  switch (test) {
    case 1: return AddLoopNonAtomic<T>(loop);
    case 2: return CASLoopNonAtomic<T>(loop);
    case 3: return ExchLoopNonAtomic<T>(loop);
  }
  return 0;
}

template <class T>
void Loop(benchmark::State& state, int test) {
  int loop = state.range(0);

  LoopVariant<T>(loop, test);

  Kokkos::Timer timer;
  T res       = LoopVariant<T>(loop, test);
  double time = timer.seconds();

  timer.reset();
  T resNonAtomic       = LoopVariantNonAtomic<T>(loop, test);
  double timeNonAtomic = timer.seconds();

  timer.reset();
  T resSerial       = LoopVariantSerial<T>(loop, test);
  double timeSerial = timer.seconds();

  time *= 1e6 / loop;
  timeNonAtomic *= 1e6 / loop;
  timeSerial *= 1e6 / loop;

  bool passed = (resSerial == res);

  state.counters["Passed"]           = benchmark::Counter(passed);
  state.counters["Value serial"]     = benchmark::Counter(resSerial);
  state.counters["Value atomic"]     = benchmark::Counter(res);
  state.counters["Value non-atomic"] = benchmark::Counter(resNonAtomic);
  state.counters["Time serial"]      = benchmark::Counter(timeSerial);
  state.counters["Time atomic"]      = benchmark::Counter(time);
  state.counters["Time non-atomic"]  = benchmark::Counter(timeNonAtomic);
  state.counters["Size of type"]     = benchmark::Counter(sizeof(T));
}

template <class T>
static void Test_Atomic(benchmark::State& state) {
  for (auto _ : state) {
    Loop<T>(state, 1);
    Loop<T>(state, 2);
    Loop<T>(state, 3);
  }
}

static constexpr int LOOP = 100'000;

BENCHMARK(Test_Atomic<int>)->Arg(30'000)->Iterations(10);
BENCHMARK(Test_Atomic<long int>)->Arg(LOOP)->Iterations(10);
BENCHMARK(Test_Atomic<long long int>)->Arg(LOOP)->Iterations(10);
BENCHMARK(Test_Atomic<unsigned int>)->Arg(LOOP)->Iterations(10);
BENCHMARK(Test_Atomic<unsigned long int>)->Arg(LOOP)->Iterations(10);
BENCHMARK(Test_Atomic<unsigned long long int>)->Arg(LOOP)->Iterations(10);
BENCHMARK(Test_Atomic<float>)->Arg(LOOP)->Iterations(10);
BENCHMARK(Test_Atomic<double>)->Arg(LOOP)->Iterations(10);
