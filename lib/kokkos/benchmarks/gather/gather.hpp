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

template <class Scalar, int UNROLL>
struct RunGather {
  static void run(int N, int K, int D, int R, int F);
};

#define UNROLL 1
#include "gather_unroll.hpp"
#undef UNROLL
#define UNROLL 2
#include "gather_unroll.hpp"
#undef UNROLL
#define UNROLL 3
#include "gather_unroll.hpp"
#undef UNROLL
#define UNROLL 4
#include "gather_unroll.hpp"
#undef UNROLL
#define UNROLL 5
#include "gather_unroll.hpp"
#undef UNROLL
#define UNROLL 6
#include "gather_unroll.hpp"
#undef UNROLL
#define UNROLL 7
#include "gather_unroll.hpp"
#undef UNROLL
#define UNROLL 8
#include "gather_unroll.hpp"
#undef UNROLL

template <class Scalar>
void run_gather_test(int N, int K, int D, int R, int U, int F) {
  if (U == 1) RunGather<Scalar, 1>::run(N, K, D, R, F);
  if (U == 2) RunGather<Scalar, 2>::run(N, K, D, R, F);
  if (U == 3) RunGather<Scalar, 3>::run(N, K, D, R, F);
  if (U == 4) RunGather<Scalar, 4>::run(N, K, D, R, F);
  if (U == 5) RunGather<Scalar, 5>::run(N, K, D, R, F);
  if (U == 6) RunGather<Scalar, 6>::run(N, K, D, R, F);
  if (U == 7) RunGather<Scalar, 7>::run(N, K, D, R, F);
  if (U == 8) RunGather<Scalar, 8>::run(N, K, D, R, F);
}
