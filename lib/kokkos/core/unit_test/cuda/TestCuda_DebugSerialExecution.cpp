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
#include <TestCuda_Category.hpp>

namespace Test {

using ViewType = Kokkos::View<double*>;

struct TestForFunctor {
  ViewType a;
  ViewType b;

  TestForFunctor(int N) : a(ViewType("A", N)), b(ViewType("B", N)) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { a(i) = b(i); }

  double time_par_for() {
    Kokkos::Timer timer;
    Kokkos::parallel_for("CudaDebugSerialExecution::par_for", a.extent(0),
                         *this);
    Kokkos::fence();
    return timer.seconds();
  }
};

struct TestRedFunctor {
  ViewType a;
  ViewType b;

  TestRedFunctor(int N) : a(ViewType("A", N)), b(ViewType("B", N)) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, double& val) const { val += a(i) * b(i); }

  double time_par_red() {
    Kokkos::Timer timer;
    double dot;
    Kokkos::parallel_reduce("CudaDebugSerialExecution::par_red", a.extent(0),
                            *this, dot);
    Kokkos::fence();
    return timer.seconds();
  }
};

struct TestScanFunctor {
  ViewType a;
  ViewType b;

  TestScanFunctor(int N) : a(ViewType("A", N)), b(ViewType("B", N)) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, double& val, bool final) const {
    val += b(i);
    if (final) a(i) = val;
  }

  double time_par_scan() {
    Kokkos::Timer timer;
    double dot;
    Kokkos::parallel_scan("CudaDebugSerialExecution::par_scan", a.extent(0),
                          *this, dot);
    Kokkos::fence();
    return timer.seconds();
  }
};

TEST(cuda, debug_serial_execution) {
  double time_par_for_1, time_par_for_2, time_par_for_serial;
  double time_par_red_1, time_par_red_2, time_par_red_serial;
  double time_par_scan_1, time_par_scan_2, time_par_scan_serial;

  int N = 10000000;
  {
    TestForFunctor f(N);
    f.time_par_for();
    time_par_for_1 = f.time_par_for();
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    kokkos_impl_cuda_set_serial_execution(true);
#endif
    time_par_for_serial = f.time_par_for();
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    kokkos_impl_cuda_set_serial_execution(false);
#endif
    time_par_for_2 = f.time_par_for();

    bool passed_par_for =
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
        (time_par_for_serial > time_par_for_1 * 20.0) &&
        (time_par_for_serial > time_par_for_2 * 20.0);
#else
        (time_par_for_serial < time_par_for_1 * 2.0) &&
        (time_par_for_serial < time_par_for_2 * 2.0);
#endif
    if (!passed_par_for)
      printf("Time For1: %lf For2: %lf ForSerial: %lf\n", time_par_for_1,
             time_par_for_2, time_par_for_serial);
    ASSERT_TRUE(passed_par_for);
  }
  {
    TestRedFunctor f(N);
    f.time_par_red();
    time_par_red_1 = f.time_par_red();
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    kokkos_impl_cuda_set_serial_execution(true);
#endif
    time_par_red_serial = f.time_par_red();
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    kokkos_impl_cuda_set_serial_execution(false);
#endif
    time_par_red_2 = f.time_par_red();

    bool passed_par_red =
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
        (time_par_red_serial > time_par_red_1 * 2.0) &&
        (time_par_red_serial > time_par_red_2 * 2.0);
#else
        (time_par_red_serial < time_par_red_1 * 2.0) &&
        (time_par_red_serial < time_par_red_2 * 2.0);
#endif
    if (!passed_par_red)
      printf("Time Red1: %lf Red2: %lf RedSerial: %lf\n", time_par_red_1,
             time_par_red_2, time_par_red_serial);
    ASSERT_TRUE(passed_par_red);
  }
  {
    TestScanFunctor f(N);
    f.time_par_scan();
    time_par_scan_1 = f.time_par_scan();
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    kokkos_impl_cuda_set_serial_execution(true);
#endif
    time_par_scan_serial = f.time_par_scan();
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
    kokkos_impl_cuda_set_serial_execution(false);
#endif
    time_par_scan_2 = f.time_par_scan();

    bool passed_par_scan =
#ifdef KOKKOS_IMPL_DEBUG_CUDA_SERIAL_EXECUTION
        (time_par_scan_serial > time_par_scan_1 * 2.0) &&
        (time_par_scan_serial > time_par_scan_2 * 2.0);
#else
        (time_par_scan_serial < time_par_scan_1 * 2.0) &&
        (time_par_scan_serial < time_par_scan_2 * 2.0);
#endif
    if (!passed_par_scan)
      printf("Time Scan1: %lf Scan2: %lf ScanSerial: %lf\n", time_par_scan_1,
             time_par_scan_2, time_par_scan_serial);
    ASSERT_TRUE(passed_par_scan);
  }
}

}  // namespace Test
