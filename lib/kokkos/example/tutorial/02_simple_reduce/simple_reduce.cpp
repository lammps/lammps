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
//      using a functor to define the loop body
//   3. Shut down Kokkos
//
// Compare this example to 02_simple_reduce_lambda, which uses a C++11
// lambda to define the loop body of the parallel_reduce.
//

// Reduction functor for computing the sum of squares.
//
// More advanced reduction examples will show how to control the
// reduction's "join" operator.  If the join operator is not provided,
// it defaults to binary operator+ (adding numbers together).
struct squaresum {
  // Specify the type of the reduction value with a "value_type"
  // alias.  In this case, the reduction value has type int.
  using value_type = int;

  // The reduction functor's operator() looks a little different than
  // the parallel_for functor's operator().  For the reduction, we
  // pass in both the loop index i, and the intermediate reduction
  // value lsum.  The latter MUST be passed in by nonconst reference.
  // (If the reduction type is an array like int[], indicating an
  // array reduction result, then the second argument is just int[].)
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, int& lsum) const {
    lsum += i * i;  // compute the sum of squares
  }
};

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  const int n = 10;

  // Compute the sum of squares of integers from 0 to n-1, in
  // parallel, using Kokkos.
  int sum = 0;
  Kokkos::parallel_reduce(n, squaresum(), sum);
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
  return (sum == seqSum) ? 0 : -1;
}
