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
#include <cstdio>

//
// First reduction (parallel_reduce) example:
//   1. Start up Kokkos
//   2. Execute a parallel_reduce loop in the default execution space,
//      using a C++11 lambda to define the loop body
//   3. Shut down Kokkos
//
// This example only builds if C++11 is enabled.  Compare this example
// to 02_simple_reduce, which uses a functor to define the loop body
// of the parallel_reduce.
//

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  const int n = 10;

  // Compute the sum of squares of integers from 0 to n-1, in
  // parallel, using Kokkos.  This time, use a lambda instead of a
  // functor.  The lambda takes the same arguments as the functor's
  // operator().
  int sum = 0;
// The KOKKOS_LAMBDA macro replaces the capture-by-value clause [=].
// It also handles any other syntax needed for CUDA.
// We also need to protect the usage of a lambda against compiling
// with a backend which doesn't support it (i.e. Cuda 6.5/7.0).
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  Kokkos::parallel_reduce(
      n, KOKKOS_LAMBDA(const int i, int& lsum) { lsum += i * i; }, sum);
#endif
  printf(
      "Sum of squares of integers from 0 to %i, "
      "computed in parallel, is %i\n",
      n - 1, sum);

  // Compare to a sequential loop.
  int seqSum = 0;
  for (int i = 0; i < n; ++i) {
    seqSum += i * i;
  }
  printf(
      "Sum of squares of integers from 0 to %i, "
      "computed sequentially, is %i\n",
      n - 1, seqSum);
  Kokkos::finalize();
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  return (sum == seqSum) ? 0 : -1;
#else
  return 0;
#endif
}
