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
#include <impl/Kokkos_Timer.hpp>
#include <cstdio>

// These two View types are both 2-D arrays of double.  However, they
// have different layouts in memory.  left_type has "layout left,"
// which means "column major," the same as in Fortran, the BLAS, or
// LAPACK.  right_type has "layout right," which means "row major,"
// the same as in C, C++, or Java.
typedef Kokkos::View<double**, Kokkos::LayoutLeft> left_type;
typedef Kokkos::View<double**, Kokkos::LayoutRight> right_type;
// This is a one-dimensional View, so the layout matters less.
// However, it still has a layout!  Since its layout is not specified
// explicitly in the type, its layout is a function of the memory
// space.  For example, the default Cuda layout is LayoutLeft, and the
// default Host layout is LayoutRight.
typedef Kokkos::View<double*> view_type;

// parallel_for functor that fills the given View with some data.  It
// expects to access the View by rows in parallel: each call i of
// operator() accesses a row.
template <class ViewType>
struct init_view {
  ViewType a;
  init_view(ViewType a_) : a(a_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename ViewType::size_type i) const {
    // On CPUs this loop could be vectorized so j should do stride 1
    // access on a for optimal performance. I.e. a should be LayoutRight.
    // On GPUs threads should do coalesced loads and stores. That means
    // that i should be the stride one access for optimal performance.
    for (typename ViewType::size_type j = 0; j < a.extent(1); ++j) {
      a(i, j) = 1.0 * a.extent(0) * i + 1.0 * j;
    }
  }
};

// Compute a contraction of v1 and v2 into a:
//
//   a(i) := sum_j (v1(i,j) * v2(j,i))
//
// Since the functor is templated on the ViewTypes itself it doesn't matter what
// there layouts are. That means you can use different layouts on different
// architectures.
template <class ViewType1, class ViewType2>
struct contraction {
  view_type a;
  typename ViewType1::const_type v1;
  typename ViewType2::const_type v2;
  contraction(view_type a_, ViewType1 v1_, ViewType2 v2_)
      : a(a_), v1(v1_), v2(v2_) {}

  // As with the initialization functor the performance of this operator
  // depends on the architecture and the chosen data layouts.
  // On CPUs optimal would be to vectorize the inner loop, so j should be the
  // stride 1 access. That means v1 should be LayoutRight and v2 LayoutLeft.
  // In order to get coalesced access on GPUs where i corresponds closely to
  // the thread Index, i must be the stride 1 dimension. That means v1 should be
  // LayoutLeft and v2 LayoutRight.
  KOKKOS_INLINE_FUNCTION
  void operator()(const view_type::size_type i) const {
    for (view_type::size_type j = 0; j < v1.extent(1); ++j) {
      a(i) = v1(i, j) * v2(j, i);
    }
  }
};

// Compute a dot product. This is used for result verification.
struct dot {
  view_type a;
  dot(view_type a_) : a(a_) {}
  typedef double value_type;  // Specify type for reduction target, lsum
  KOKKOS_INLINE_FUNCTION
  void operator()(const view_type::size_type i, double& lsum) const {
    lsum += a(i) * a(i);
  }
};

int main(int narg, char* arg[]) {
  // When initializing Kokkos, you may pass in command-line arguments,
  // just like with MPI_Init().  Kokkos reserves the right to remove
  // arguments from the list that start with '--kokkos-'.
  Kokkos::initialize(narg, arg);

  {
    int size = 10000;
    view_type a("A", size);

    // Define two views with LayoutLeft and LayoutRight.
    left_type l("L", size, 10000);
    right_type r("R", size, 10000);

    // Initialize the data in the views.
    Kokkos::parallel_for(size, init_view<left_type>(l));
    Kokkos::parallel_for(size, init_view<right_type>(r));
    Kokkos::fence();

    // Measure time to execute the contraction kernel when giving it a
    // LayoutLeft view for v1 and a LayoutRight view for v2. This should be
    // fast on GPUs and slow on CPUs
    Kokkos::Timer time1;
    Kokkos::parallel_for(size, contraction<left_type, right_type>(a, l, r));
    Kokkos::fence();
    double sec1 = time1.seconds();

    double sum1 = 0;
    Kokkos::parallel_reduce(size, dot(a), sum1);
    Kokkos::fence();

    // Measure time to execute the contraction kernel when giving it a
    // LayoutRight view for v1 and a LayoutLeft view for v2. This should be
    // fast on CPUs and slow on GPUs
    Kokkos::Timer time2;
    Kokkos::parallel_for(size, contraction<right_type, left_type>(a, r, l));
    Kokkos::fence();
    double sec2 = time2.seconds();

    double sum2 = 0;
    Kokkos::parallel_reduce(size, dot(a), sum2);

    // Kokkos' reductions are deterministic.
    // The results should always be equal.
    printf("Result Left/Right %f Right/Left %f (equal result: %i)\n", sec1,
           sec2, sum2 == sum1);
  }

  Kokkos::finalize();
}
