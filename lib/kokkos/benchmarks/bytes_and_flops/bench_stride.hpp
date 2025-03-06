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

#define UNROLL 1
#include "bench_unroll_stride.hpp"
#undef UNROLL
#define UNROLL 2
#include "bench_unroll_stride.hpp"
#undef UNROLL
#define UNROLL 3
#include "bench_unroll_stride.hpp"
#undef UNROLL
#define UNROLL 4
#include "bench_unroll_stride.hpp"
#undef UNROLL
#define UNROLL 5
#include "bench_unroll_stride.hpp"
#undef UNROLL
#define UNROLL 6
#include "bench_unroll_stride.hpp"
#undef UNROLL
#define UNROLL 7
#include "bench_unroll_stride.hpp"
#undef UNROLL
#define UNROLL 8
#include "bench_unroll_stride.hpp"
#undef UNROLL

template <class Scalar>
struct RunStride<Scalar, STRIDE> {
  static void run_1(int N, int K, int R, int F, int T, int S, int B, int I) {
    Run<Scalar, 1, STRIDE>::run(N, K, R, F, T, S, B, I);
  }
  static void run_2(int N, int K, int R, int F, int T, int S, int B, int I) {
    Run<Scalar, 2, STRIDE>::run(N, K, R, F, T, S, B, I);
  }
  static void run_3(int N, int K, int R, int F, int T, int S, int B, int I) {
    Run<Scalar, 3, STRIDE>::run(N, K, R, F, T, S, B, I);
  }
  static void run_4(int N, int K, int R, int F, int T, int S, int B, int I) {
    Run<Scalar, 4, STRIDE>::run(N, K, R, F, T, S, B, I);
  }
  static void run_5(int N, int K, int R, int F, int T, int S, int B, int I) {
    Run<Scalar, 5, STRIDE>::run(N, K, R, F, T, S, B, I);
  }
  static void run_6(int N, int K, int R, int F, int T, int S, int B, int I) {
    Run<Scalar, 6, STRIDE>::run(N, K, R, F, T, S, B, I);
  }
  static void run_7(int N, int K, int R, int F, int T, int S, int B, int I) {
    Run<Scalar, 7, STRIDE>::run(N, K, R, F, T, S, B, I);
  }
  static void run_8(int N, int K, int R, int F, int T, int S, int B, int I) {
    Run<Scalar, 8, STRIDE>::run(N, K, R, F, T, S, B, I);
  }

  static void run(int N, int K, int R, int U, int F, int T, int S, int B,
                  int I) {
    if (U == 1) {
      run_1(N, K, R, F, T, S, B, I);
    }
    if (U == 2) {
      run_2(N, K, R, F, T, S, B, I);
    }
    if (U == 3) {
      run_3(N, K, R, F, T, S, B, I);
    }
    if (U == 4) {
      run_4(N, K, R, F, T, S, B, I);
    }
    if (U == 5) {
      run_5(N, K, R, F, T, S, B, I);
    }
    if (U == 6) {
      run_6(N, K, R, F, T, S, B, I);
    }
    if (U == 7) {
      run_7(N, K, R, F, T, S, B, I);
    }
    if (U == 8) {
      run_8(N, K, R, F, T, S, B, I);
    }
  }
};
