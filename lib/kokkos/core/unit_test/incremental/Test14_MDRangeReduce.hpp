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
  MyComplex& operator=(const MyComplex& src) {
    _re = src._re;
    _im = src._im;
    return *this;
  }

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

// FIXME_OPENACC: scalar reduction variable on the device is not yet supported.
#if !defined(KOKKOS_ENABLE_OPENACC)
    // Parallel reduce on a view.
    Kokkos::parallel_reduce(
        mdPolicy_2D,
        KOKKOS_LAMBDA(const int i, const int j, value_type& update_value) {
          update_value += d_data(i, j);
        },
        d_resultView);
#endif

    // Check correctness.
    ASSERT_EQ(h_result, d_result);

// FIXME_OPENACC: scalar reduction variable on the device is not yet supported.
#if !defined(KOKKOS_ENABLE_OPENACC)
    // Copy view back to host.
    value_type view_result = 0.0;
    Kokkos::deep_copy(view_result, d_resultView);
    ASSERT_EQ(h_result, view_result);
#endif
  }

// FIXME_OPENACC: custom reductions are not yet supported in the
// OpenACC backend.
#if !defined(KOKKOS_ENABLE_OPENACC)
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
#endif
};

// Reductions tests for MDRange policy and customized reduction.
TEST(TEST_CATEGORY, incr_14_MDrangeReduce) {
  TestMDRangeReduce<TEST_EXECSPACE> test;
  test.reduce_MDRange();
// FIXME_OPENMPTARGET: custom reductions are not yet supported in the
// OpenMPTarget backend.
// FIXME_OPENACC: custom reductions are not yet supported in the
// OpenACC backend.
#if !defined(KOKKOS_ENABLE_OPENMPTARGET)
#if !defined(KOKKOS_ENABLE_OPENACC)
  test.reduce_custom();
#endif
#endif
}

}  // namespace Test
