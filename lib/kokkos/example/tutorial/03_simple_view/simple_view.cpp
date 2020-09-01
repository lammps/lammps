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

//
// First Kokkos::View (multidimensional array) example:
//   1. Start up Kokkos
//   2. Allocate a Kokkos::View
//   3. Execute a parallel_for and a parallel_reduce over that View's data
//   4. Shut down Kokkos
//
// Compare this example to 03_simple_view_lambda, which uses C++11
// lambdas to define the loop bodies of the parallel_for and
// parallel_reduce.
//

#include <Kokkos_Core.hpp>
#include <cstdio>

// A Kokkos::View is an array of zero or more dimensions.  The number
// of dimensions is specified at compile time, as part of the type of
// the View.  This array has two dimensions.  The first one
// (represented by the asterisk) is a run-time dimension, and the
// second (represented by [3]) is a compile-time dimension.  Thus,
// this View type is an N x 3 array of type double, where N is
// specified at run time in the View's constructor.
//
// The first dimension of the View is the dimension over which it is
// efficient for Kokkos to parallelize.
using view_type = Kokkos::View<double * [3]>;

// parallel_for functor that fills the View given to its constructor.
// The View must already have been allocated.
struct InitView {
  view_type a;

  // Views have "view semantics."  This means that they behave like
  // pointers, not like std::vector.  Their copy constructor and
  // operator= only do shallow copies.  Thus, you can pass View
  // objects around by "value"; they won't do a deep copy unless you
  // explicitly ask for a deep copy.
  InitView(view_type a_) : a(a_) {}

  // Fill the View with some data.  The parallel_for loop will iterate
  // over the View's first dimension N.
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    // Acesss the View just like a Fortran array.  The layout depends
    // on the View's memory space, so don't rely on the View's
    // physical memory layout unless you know what you're doing.
    a(i, 0) = 1.0 * i;
    a(i, 1) = 1.0 * i * i;
    a(i, 2) = 1.0 * i * i * i;
  }
};

// Reduction functor that reads the View given to its constructor.
struct ReduceFunctor {
  view_type a;

  // Constructor takes View by "value"; this does a shallow copy.
  ReduceFunctor(view_type a_) : a(a_) {}

  // If you write a functor to do a reduction, you must specify the
  // type of the reduction result via a public 'value_type' alias.
  using value_type = double;

  KOKKOS_INLINE_FUNCTION
  void operator()(int i, double& lsum) const {
    lsum += a(i, 0) * a(i, 1) / (a(i, 2) + 0.1);
  }
};

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
    const int N = 10;

    // Allocate the View.  The first dimension is a run-time parameter
    // N.  We set N = 10 here.  The second dimension is a compile-time
    // parameter, 3.  We don't specify it here because we already set it
    // by declaring the type of the View.
    //
    // Views get initialized to zero by default.  This happens in
    // parallel, using the View's memory space's default execution
    // space.  Parallel initialization ensures first-touch allocation.
    // There is a way to shut off default initialization.
    //
    // You may NOT allocate a View inside of a parallel_{for, reduce,
    // scan}.  Treat View allocation as a "thread collective."
    //
    // The string "A" is just the label; it only matters for debugging.
    // Different Views may have the same label.
    view_type a("A", N);

    Kokkos::parallel_for(N, InitView(a));
    double sum = 0;
    Kokkos::parallel_reduce(N, ReduceFunctor(a), sum);
    printf("Result: %f\n", sum);
  }  // use this scope to ensure the lifetime of "A" ends before finalize
  Kokkos::finalize();
}
