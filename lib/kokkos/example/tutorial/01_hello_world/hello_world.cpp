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
#include <typeinfo>

//
// "Hello world" parallel_for example:
//   1. Start up Kokkos
//   2. Execute a parallel for loop in the default execution space,
//      using a functor to define the loop body
//   3. Shut down Kokkos
//
// If Kokkos was built with C++11 enabled, try comparing this example
// to 01_hello_world_lambda.  The latter uses C++11 lambdas (anonymous
// functions) to define the loop body of the parallel_for.  That makes
// the code much more concise and readable.  On the other hand,
// breaking out the loop body into an explicit functor makes it easier
// to test the loop independently of the parallel pattern.
//

// Functor that defines the parallel_for's loop body.
//
// A "functor" is just a class or struct with a public operator()
// instance method.
struct hello_world {
  // If a functor has an "execution_space" (or "execution_space", for
  // backwards compatibility) public typedef, parallel_* will only run
  // the functor in that execution space.  That's a good way to mark a
  // functor as specific to an execution space.  If the functor lacks
  // this typedef, parallel_for will run it in the default execution
  // space, unless you tell it otherwise (that's an advanced topic;
  // see "execution policies").

  // The functor's operator() defines the loop body.  It takes an
  // integer argument which is the parallel for loop index.  Other
  // arguments are possible; see the "hierarchical parallelism" part
  // of the tutorial.
  //
  // The operator() method must be const, and must be marked with the
  // KOKKOS_INLINE_FUNCTION macro.  If building with CUDA, this macro
  // will mark your method as suitable for running on the CUDA device
  // (as well as on the host).  If not building with CUDA, the macro
  // is unnecessary but harmless.
  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {
    printf ("Hello from i = %i\n", i);
  }
};

int main (int argc, char* argv[]) {
  // You must call initialize() before you may call Kokkos.
  //
  // With no arguments, this initializes the default execution space
  // (and potentially its host execution space) with default
  // parameters.  You may also pass in argc and argv, analogously to
  // MPI_Init().  It reads and removes command-line arguments that
  // start with "--kokkos-".
  Kokkos::initialize (argc, argv);

  // Print the name of Kokkos' default execution space.  We're using
  // typeid here, so the name might get a bit mangled by the linker,
  // but you should still be able to figure out what it is.
  printf ("Hello World on Kokkos execution space %s\n",
          typeid (Kokkos::DefaultExecutionSpace).name ());

  // Run the above functor on the default Kokkos execution space in
  // parallel, with a parallel for loop count of 15.
  //
  // The Kokkos::DefaultExecutionSpace typedef gives the default
  // execution space.  Depending on how Kokkos was configured, this
  // could be OpenMP, Threads, Cuda, Serial, or even some other
  // execution space.
  //
  // The following line of code would look like this in OpenMP:
  //
  // #pragma omp parallel for
  // for (int i = 0; i < 15; ++i) {
  //   printf ("Hello from i = %i\n", i);
  // }
  //
  // You may notice that the printed numbers do not print out in
  // order.  Parallel for loops may execute in any order.
  Kokkos::parallel_for ("HelloWorld",15, hello_world ());

  // You must call finalize() after you are done using Kokkos.
  Kokkos::finalize ();
}

