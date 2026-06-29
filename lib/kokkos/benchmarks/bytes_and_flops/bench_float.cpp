// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include "bench.hpp"

template void run_stride_unroll<float>(int N, int K, int R, int D, int U, int F,
                                       int T, int S, int B, int I);
