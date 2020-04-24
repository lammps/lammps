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
#include <cstdio>
#include <typeinfo>

//
// "Hello world" parallel_for example:
//   1. Start up Kokkos
//   2. Execute a parallel for loop in the default execution space,
//      using a C++11 lambda to define the loop body
//   3. Shut down Kokkos
//
// This example only builds if C++11 is enabled.  Compare this example
// to 01_hello_world, which uses functors (explicitly defined classes)
// to define the loop body of the parallel_for.  Both functors and
// lambdas have their places.
//

int main(int argc, char* argv[]) {
  // You must call initialize() before you may call Kokkos.
  //
  // With no arguments, this initializes the default execution space
  // (and potentially its host execution space) with default
  // parameters.  You may also pass in argc and argv, analogously to
  // MPI_Init().  It reads and removes command-line arguments that
  // start with "--kokkos-".
  Kokkos::initialize(argc, argv);

  // Print the name of Kokkos' default execution space.  We're using
  // typeid here, so the name might get a bit mangled by the linker,
  // but you should still be able to figure out what it is.
  printf("Hello World on Kokkos execution space %s\n",
         typeid(Kokkos::DefaultExecutionSpace).name());

  // Run lambda on the default Kokkos execution space in parallel,
  // with a parallel for loop count of 15.  The lambda's argument is
  // an integer which is the parallel for's loop index.  As you learn
  // about different kinds of parallelism, you will find out that
  // there are other valid argument types as well.
  //
  // For a single level of parallelism, we prefer that you use the
  // KOKKOS_LAMBDA macro.  If CUDA is disabled, this just turns into
  // [=].  That captures variables from the surrounding scope by
  // value.  Do NOT capture them by reference!  If CUDA is enabled,
  // this macro may have a special definition that makes the lambda
  // work correctly with CUDA.  Compare to the KOKKOS_INLINE_FUNCTION
  // macro, which has a special meaning if CUDA is enabled.
  //
  // The following parallel_for would look like this if we were using
  // OpenMP by itself, instead of Kokkos:
  //
  // #pragma omp parallel for
  // for (int i = 0; i < 15; ++i) {
  //   printf ("Hello from i = %i\n", i);
  // }
  //
  // You may notice that the printed numbers do not print out in
  // order.  Parallel for loops may execute in any order.
  // We also need to protect the usage of a lambda against compiling
  // with a backend which doesn't support it (i.e. Cuda 6.5/7.0).
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  Kokkos::parallel_for(
      15, KOKKOS_LAMBDA(const int i) {
        // printf works in a CUDA parallel kernel; std::ostream does not.
        printf("Hello from i = %i\n", i);
      });
#endif
  // You must call finalize() after you are done using Kokkos.
  Kokkos::finalize();
}
