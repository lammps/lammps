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
//      using a C++11 lambda to define the loop body
//   3. Shut down Kokkos
//
// This example only builds if C++11 is enabled.  Compare this example
// to 02_simple_reduce, which uses a functor to define the loop body
// of the parallel_reduce.
//

int main (int argc, char* argv[]) {
  Kokkos::initialize (argc, argv);
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
  Kokkos::parallel_reduce (n, KOKKOS_LAMBDA (const int i, int& lsum) {
      lsum += i*i;
    }, sum);
  #endif
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
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  return (sum == seqSum) ? 0 : -1;
#else
  return 0;
#endif
}

