/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
      double* a_ptr = (double*)Kokkos::kokkos_malloc("A", sizeof(double) * N8);
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
