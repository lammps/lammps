/*
//@HEADER
// ************************************************************************
// 
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

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
  // typedef.  In this case, the reduction value has type int.
  typedef int value_type;

  // The reduction functor's operator() looks a little different than
  // the parallel_for functor's operator().  For the reduction, we
  // pass in both the loop index i, and the intermediate reduction
  // value lsum.  The latter MUST be passed in by nonconst reference.
  // (If the reduction type is an array like int[], indicating an
  // array reduction result, then the second argument is just int[].)
  KOKKOS_INLINE_FUNCTION
  void operator () (const int i, int& lsum) const {
    lsum += i*i; // compute the sum of squares
  }
};

int main (int argc, char* argv[]) {
  Kokkos::initialize (argc, argv);
  const int n = 10;

  // Compute the sum of squares of integers from 0 to n-1, in
  // parallel, using Kokkos.
  int sum = 0;
  Kokkos::parallel_reduce (n, squaresum (), sum);
  printf ("Sum of squares of integers from 0 to %i, "
          "computed in parallel, is %i\n", n - 1, sum);

  // Compare to a sequential loop.
  int seqSum = 0;
  for (int i = 0; i < n; ++i) {
    seqSum += i*i;
  }
  printf ("Sum of squares of integers from 0 to %i, "
          "computed sequentially, is %i\n", n - 1, seqSum);
  Kokkos::finalize ();
  return (sum == seqSum) ? 0 : -1;
}

