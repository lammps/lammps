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

/// @Kokkos_Feature_Level_Required:14
// Incremental test for MDRange reduction .
// Reduction is tested with scalar, view and a customized reduction.

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

namespace Test {
using value_type = double;
const int N      = 10;
const int M      = 10;

// A structure for complex number.
struct MyComplex {
  value_type _re, _im;

  MyComplex() = default;

  KOKKOS_INLINE_FUNCTION
  MyComplex(value_type re, value_type im) : _re(re), _im(im) {}

  KOKKOS_INLINE_FUNCTION
  MyComplex(const MyComplex& src) : _re(src._re), _im(src._im) {}

  KOKKOS_INLINE_FUNCTION
  void operator+=(const MyComplex& src) {
    _re += src._re;
    _im += src._im;
  }
};

template <class ExecSpace>
struct TestMDRangeReduce {
  // 1D  View of double
  using View_1D = Kokkos::View<value_type*, ExecSpace>;

  // 2D  View of double
  using View_2D = Kokkos::View<value_type**, ExecSpace>;

  // Index Type for the iterator
  using int_index = Kokkos::IndexType<int>;

  // An MDRangePolicy for 2 nested loops
  using MDPolicyType_2D =
      Kokkos::MDRangePolicy<ExecSpace, Kokkos::Rank<2>, int_index>;

  //  1D - complex View
  using Complex_View_1D = Kokkos::View<MyComplex*, ExecSpace>;

  // Reduction when ExecPolicy = MDRangePolicy and ReducerArgument =
  // scalar/1-element view
  void reduce_MDRange() {
    View_2D d_data("d_data", N, M);

    MDPolicyType_2D mdPolicy_2D({0, 0}, {N, M});

    // Store the reduced value.
    value_type d_result = 0.0, h_result = 0.0;
    Kokkos::View<value_type, ExecSpace> d_resultView("result View");

    // Compute reference solution on the host.
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j) h_result += i * j;
    h_result *= 0.5;

    // Fill data.
    Kokkos::parallel_for(
        mdPolicy_2D, KOKKOS_LAMBDA(const int i, const int j) {
          d_data(i, j) = i * j * 0.5;
        });

    // Parallel reduce on a scalar.
    Kokkos::parallel_reduce(
        mdPolicy_2D,
        KOKKOS_LAMBDA(const int i, const int j, value_type& update_value) {
          update_value += d_data(i, j);
        },
        d_result);

    // Parallel reduce on a view.
    Kokkos::parallel_reduce(
        mdPolicy_2D,
        KOKKOS_LAMBDA(const int i, const int j, value_type& update_value) {
          update_value += d_data(i, j);
        },
        d_resultView);

    // Check correctness.
    ASSERT_EQ(h_result, d_result);

    // Copy view back to host.
    value_type view_result = 0.0;
    Kokkos::deep_copy(view_result, d_resultView);
    ASSERT_EQ(h_result, view_result);
  }

  // Custom Reduction
  void reduce_custom() {
    Complex_View_1D d_data("complex array", N);
    MyComplex result(0.0, 0.0);
    int sum = 0;

    // Fill data
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, N), KOKKOS_LAMBDA(const int i) {
          d_data(i) = MyComplex(i * 0.5, -i * 0.5);
        });

    // Reduction for complex number.
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<ExecSpace>(0, N),
        KOKKOS_LAMBDA(const int i, MyComplex& update_value) {
          update_value += d_data(i);
        },
        result);

    // Correctness Check
    for (int i = 0; i < N; ++i) sum += i;

    ASSERT_EQ(result._re, sum * 0.5);
    ASSERT_EQ(result._im, -sum * 0.5);
  }
};

// Reductions tests for MDRange policy and customized reduction.
TEST(TEST_CATEGORY, incr_14_MDrangeReduce) {
  TestMDRangeReduce<TEST_EXECSPACE> test;
  test.reduce_MDRange();
// FIXME_OPENMPTARGET: custom reductions are not yet supported in the
// OpenMPTarget backend.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
  test.reduce_custom();
#endif
}

}  // namespace Test
