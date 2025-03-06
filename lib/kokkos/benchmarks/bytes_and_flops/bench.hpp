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
#include <Kokkos_Timer.hpp>

template <class Scalar, int Unroll, int Stride>
struct Run {
  static void run(int N, int K, int R, int F, int T, int S, int B, int I);
};

template <class Scalar, int Stride>
struct RunStride {
  static void run_1(int N, int K, int R, int F, int T, int S, int B, int I);
  static void run_2(int N, int K, int R, int F, int T, int S, int B, int I);
  static void run_3(int N, int K, int R, int F, int T, int S, int B, int I);
  static void run_4(int N, int K, int R, int F, int T, int S, int B, int I);
  static void run_5(int N, int K, int R, int F, int T, int S, int B, int I);
  static void run_6(int N, int K, int R, int F, int T, int S, int B, int I);
  static void run_7(int N, int K, int R, int F, int T, int S, int B, int I);
  static void run_8(int N, int K, int R, int F, int T, int S, int B, int I);
  static void run(int N, int K, int R, int U, int F, int T, int S, int B,
                  int I);
};

#define STRIDE 1
#include "bench_stride.hpp"
#undef STRIDE
#define STRIDE 2
#include "bench_stride.hpp"
#undef STRIDE
#define STRIDE 4
#include "bench_stride.hpp"
#undef STRIDE
#define STRIDE 8
#include "bench_stride.hpp"
#undef STRIDE
#define STRIDE 16
#include "bench_stride.hpp"
#undef STRIDE
#define STRIDE 32
#include "bench_stride.hpp"
#undef STRIDE

template <class Scalar>
void run_stride_unroll(int N, int K, int R, int D, int U, int F, int T, int S,
                       int B, int I) {
  if (D == 1) RunStride<Scalar, 1>::run(N, K, R, U, F, T, S, B, I);
  if (D == 2) RunStride<Scalar, 2>::run(N, K, R, U, F, T, S, B, I);
  if (D == 4) RunStride<Scalar, 4>::run(N, K, R, U, F, T, S, B, I);
  if (D == 8) RunStride<Scalar, 8>::run(N, K, R, U, F, T, S, B, I);
  if (D == 16) RunStride<Scalar, 16>::run(N, K, R, U, F, T, S, B, I);
  if (D == 32) RunStride<Scalar, 32>::run(N, K, R, U, F, T, S, B, I);
}
