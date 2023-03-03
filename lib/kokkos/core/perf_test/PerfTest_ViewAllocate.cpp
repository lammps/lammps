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
#include <gtest/gtest.h>
#include <cstdio>
#include <PerfTest_Category.hpp>

namespace Test {

template <class Layout>
void run_allocateview_tests(int N, int R) {
  const int N1 = N;
  const int N2 = N * N;
  const int N3 = N2 * N;
  const int N4 = N2 * N2;
  const int N8 = N4 * N4;

  double time1, time2, time3, time4, time5, time6, time7, time8,
      time_raw = 100000.0;
  {
    Kokkos::Timer timer;
    for (int r = 0; r < R; r++) {
      Kokkos::View<double*, Layout> a("A1", N8);
    }
    time1 = timer.seconds() / R;
  }
  {
    Kokkos::Timer timer;
    for (int r = 0; r < R; r++) {
      Kokkos::View<double**, Layout> a("A2", N4, N4);
    }
    time2 = timer.seconds() / R;
  }
  {
    Kokkos::Timer timer;
    for (int r = 0; r < R; r++) {
      Kokkos::View<double***, Layout> a("A3", N3, N3, N2);
    }
    time3 = timer.seconds() / R;
  }
  {
    Kokkos::Timer timer;
    for (int r = 0; r < R; r++) {
      Kokkos::View<double****, Layout> a("A4", N2, N2, N2, N2);
    }
    time4 = timer.seconds() / R;
  }
  {
    Kokkos::Timer timer;
    for (int r = 0; r < R; r++) {
      Kokkos::View<double*****, Layout> a("A5", N2, N2, N1, N1, N2);
    }
    time5 = timer.seconds() / R;
  }
  {
    Kokkos::Timer timer;
    for (int r = 0; r < R; r++) {
      Kokkos::View<double******, Layout> a("A6", N2, N1, N1, N1, N1, N2);
    }
    time6 = timer.seconds() / R;
  }
  {
    Kokkos::Timer timer;
    for (int r = 0; r < R; r++) {
      Kokkos::View<double*******, Layout> a("A7", N2, N1, N1, N1, N1, N1, N1);
    }
    time7 = timer.seconds() / R;
  }
  {
    Kokkos::Timer timer;
    for (int r = 0; r < R; r++) {
      Kokkos::View<double********, Layout> a("A8", N1, N1, N1, N1, N1, N1, N1,
                                             N1);
    }
    time8 = timer.seconds() / R;
  }
#if defined(KOKKOS_ENABLE_CUDA_LAMBDA) || !defined(KOKKOS_ENABLE_CUDA)
  {
    Kokkos::Timer timer;
    for (int r = 0; r < R; r++) {
      double* a_ptr =
          static_cast<double*>(Kokkos::kokkos_malloc("A", sizeof(double) * N8));
      Kokkos::parallel_for(
          N8, KOKKOS_LAMBDA(const int& i) { a_ptr[i] = 0.0; });
      Kokkos::fence();
      Kokkos::kokkos_free(a_ptr);
    }
    time_raw = timer.seconds() / R;
  }
#endif
  double size = 1.0 * N8 * 8 / 1024 / 1024;
  printf("   Raw:   %lf s   %lf MB   %lf GB/s\n", time_raw, size,
         size / 1024 / time_raw);
  printf("   Rank1: %lf s   %lf MB   %lf GB/s\n", time1, size,
         size / 1024 / time1);
  printf("   Rank2: %lf s   %lf MB   %lf GB/s\n", time2, size,
         size / 1024 / time2);
  printf("   Rank3: %lf s   %lf MB   %lf GB/s\n", time3, size,
         size / 1024 / time3);
  printf("   Rank4: %lf s   %lf MB   %lf GB/s\n", time4, size,
         size / 1024 / time4);
  printf("   Rank5: %lf s   %lf MB   %lf GB/s\n", time5, size,
         size / 1024 / time5);
  printf("   Rank6: %lf s   %lf MB   %lf GB/s\n", time6, size,
         size / 1024 / time6);
  printf("   Rank7: %lf s   %lf MB   %lf GB/s\n", time7, size,
         size / 1024 / time7);
  printf("   Rank8: %lf s   %lf MB   %lf GB/s\n", time8, size,
         size / 1024 / time8);
}

TEST(default_exec, ViewCreate) {
  printf("Create View Performance for LayoutLeft:\n");
  run_allocateview_tests<Kokkos::LayoutLeft>(10, 1);
  printf("Create View Performance for LayoutRight:\n");
  run_allocateview_tests<Kokkos::LayoutRight>(10, 1);
}

}  // namespace Test
