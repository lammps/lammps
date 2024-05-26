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
#include "bench.hpp"
#include <cstdlib>

extern template void run_stride_unroll<float>(int, int, int, int, int, int, int,
                                              int, int, int);
extern template void run_stride_unroll<double>(int, int, int, int, int, int,
                                               int, int, int, int);
extern template void run_stride_unroll<int32_t>(int, int, int, int, int, int,
                                                int, int, int, int);
extern template void run_stride_unroll<int64_t>(int, int, int, int, int, int,
                                                int, int, int, int);

int main(int argc, char* argv[]) {
  Kokkos::initialize();

  if (argc < 10) {
    printf("Arguments: N K R D U F T S B I\n");
    printf("  P:   Precision (1==float, 2==double, 3==int32_t, 4==int64_t)\n");
    printf("  N,K: dimensions of the 2D array to allocate\n");
    printf("  R:   how often to loop through the K dimension with each team\n");
    printf("  D:   distance between loaded elements (stride)\n");
    printf("  U:   how many independent flops to do per load\n");
    printf(
        "  F:   how many times to repeat the U unrolled operations before "
        "reading next element\n");
    printf("  T:   team size\n");
    printf(
        "  S:   shared memory per team (used to control occupancy on GPUs)\n");
    printf(
        "  B:   units for reported memory bandwidths (2=GiB, 10=GB, "
        "default=2)\n");
    printf("  I:   iterations of the kernel to time over (default=10)\n");
    printf("Example Input GPU:\n");
    printf("  Bandwidth Bound : 2 100000 1024 1 1 1 1 256 6000\n");
    printf("  Cache Bound     : 2 100000 1024 64 1 1 1 512 20000\n");
    printf("  Compute Bound   : 2 100000 1024 1 1 8 64 256 6000\n");
    printf("  Load Slots Used : 2 20000 256 32 16 1 1 256 6000\n");
    printf("  Inefficient Load: 2 20000 256 32 2 1 1 256 20000\n");
    Kokkos::finalize();
    return 0;
  }

  int P = std::stoi(argv[1]);
  int N = std::stoi(argv[2]);
  int K = std::stoi(argv[3]);
  int R = std::stoi(argv[4]);
  int D = std::stoi(argv[5]);
  int U = std::stoi(argv[6]);
  int F = std::stoi(argv[7]);
  int T = std::stoi(argv[8]);
  int S = std::stoi(argv[9]);

  int B = 2;
  if (argc >= 11) {
    B = std::atoi(argv[10]);
  }

  int I = 10;
  if (argc >= 12) {
    I = std::atoi(argv[11]);
  }

  if (U > 8) {
    printf("U must be 1-8\n");
    return 0;
  }
  if ((D != 1) && (D != 2) && (D != 4) && (D != 8) && (D != 16) && (D != 32)) {
    printf("D must be one of 1,2,4,8,16,32\n");
    return 0;
  }
  if ((P < 1) || (P > 4)) {
    printf("P must be one of 1,2,3,4\n");
    return 0;
  }

  if ((B != 2) && (B != 10)) {
    printf("B must be one of 2,10\n");
    return 0;
  }

  if (I < 1) {
    printf("I must be >= 1\n");
    return 0;
  }

  if (P == 1) {
    run_stride_unroll<float>(N, K, R, D, U, F, T, S, B, I);
  }
  if (P == 2) {
    run_stride_unroll<double>(N, K, R, D, U, F, T, S, B, I);
  }
  if (P == 3) {
    run_stride_unroll<int32_t>(N, K, R, D, U, F, T, S, B, I);
  }
  if (P == 4) {
    run_stride_unroll<int64_t>(N, K, R, D, U, F, T, S, B, I);
  }

  Kokkos::finalize();
}
