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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace Test {

TEST(TEST_CATEGORY, view_remap) {
  enum { N0 = 3, N1 = 2, N2 = 8, N3 = 9 };

#ifdef KOKKOS_ENABLE_CUDA
#define EXECSPACE                                                     \
  std::conditional<std::is_same<TEST_EXECSPACE, Kokkos::Cuda>::value, \
                   Kokkos::CudaHostPinnedSpace, TEST_EXECSPACE>::type
#else
#ifdef KOKKOS_ENABLE_ROCM
#define EXECSPACE                                                      \
  std::conditional<                                                    \
      std::is_same<TEST_EXECSPACE, Kokkos::Experimental::ROCm>::value, \
      Kokkos::Experimental::ROCmHostPinnedSpace, TEST_EXECSPACE>::type
#else
#if defined(KOKKOS_ENABLE_OPENMPTARGET)
#define EXECSPACE Kokkos::HostSpace
#else
#define EXECSPACE TEST_EXECSPACE
#endif
#endif
#endif

  typedef Kokkos::View<double * [N1][N2][N3], Kokkos::LayoutRight, EXECSPACE>
      output_type;

  typedef Kokkos::View<int* * [N2][N3], Kokkos::LayoutLeft, EXECSPACE>
      input_type;

  typedef Kokkos::View<int * [N0][N2][N3], Kokkos::LayoutLeft, EXECSPACE>
      diff_type;

  output_type output("output", N0);
  input_type input("input", N0, N1);
  diff_type diff("diff", N0);

  Kokkos::fence();
  int value = 0;

  for (size_t i3 = 0; i3 < N3; ++i3)
    for (size_t i2 = 0; i2 < N2; ++i2)
      for (size_t i1 = 0; i1 < N1; ++i1)
        for (size_t i0 = 0; i0 < N0; ++i0) {
          input(i0, i1, i2, i3) = ++value;
        }

  Kokkos::fence();
  // Kokkos::deep_copy( diff, input ); // Throw with incompatible shape.
  Kokkos::deep_copy(output, input);
  Kokkos::fence();

  value = 0;

  for (size_t i3 = 0; i3 < N3; ++i3)
    for (size_t i2 = 0; i2 < N2; ++i2)
      for (size_t i1 = 0; i1 < N1; ++i1)
        for (size_t i0 = 0; i0 < N0; ++i0) {
          ++value;
          ASSERT_EQ(value, ((int)output(i0, i1, i2, i3)));
        }
}

TEST(TEST_CATEGORY, view_mirror_nonconst) {
  Kokkos::View<int*, TEST_EXECSPACE> d_view("d_view", 10);
  Kokkos::View<const int*, TEST_EXECSPACE> d_view_const = d_view;
  auto h_view = Kokkos::create_mirror(d_view_const);
  Kokkos::deep_copy(h_view, d_view_const);
  auto h_view2 = Kokkos::create_mirror(Kokkos::HostSpace(), d_view_const);
  Kokkos::deep_copy(h_view2, d_view_const);
}

template <typename DataType, typename... Extents>
void test_left_stride(Extents... extents) {
  using view_type =
      Kokkos::View<DataType, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  view_type view("view", extents...);
  size_t expected_stride = 1;
  size_t all_strides[view_type::rank + 1];
  view.stride(all_strides);
  for (int i = 0; i < view_type::rank; ++i) {
    ASSERT_EQ(view.stride(i), expected_stride);
    ASSERT_EQ(all_strides[i], expected_stride);
    expected_stride *= view.extent(i);
  }
}

template <typename DataType, typename... Extents>
void test_right_stride(Extents... extents) {
  using view_type =
      Kokkos::View<DataType, Kokkos::LayoutRight, Kokkos::HostSpace>;
  view_type view("view", extents...);
  size_t expected_stride = 1;
  size_t all_strides[view_type::rank + 1];
  view.stride(all_strides);
  for (int ri = 0; ri < view_type::rank; ++ri) {
    auto i = view_type::rank - 1 - ri;
    ASSERT_EQ(view.stride(i), expected_stride);
    ASSERT_EQ(all_strides[i], expected_stride);
    expected_stride *= view.extent(i);
  }
}

template <typename DataType, typename... Extents>
void test_stride(Extents... extents) {
  test_right_stride<DataType>(extents...);
  test_left_stride<DataType>(extents...);
}

TEST(TEST_CATEGORY, view_stride_method) {
  test_stride<double[3]>();
  test_stride<double*>(3);
  test_stride<double[3][7][13]>();
  test_stride<double***>(3, 7, 13);
  // factorial(8) = 40320
  test_stride<double[1][2][3][4][5][6][7][8]>();
  test_stride<double********>(1, 2, 3, 4, 5, 6, 7, 8);
}

inline void test_anonymous_space() {
  /* apparently TEST_EXECSPACE is sometimes a memory space. */
  using ExecSpace = TEST_EXECSPACE::execution_space;
  int host_array[10];
  Kokkos::View<int[10], Kokkos::AnonymousSpace> host_anon_stat_view(host_array);
  Kokkos::View<int*, Kokkos::AnonymousSpace> host_anon_dyn_view(host_array, 10);
  Kokkos::View<int*, Kokkos::HostSpace> host_view("host_view", 10);
  Kokkos::View<int*, Kokkos::AnonymousSpace> host_anon_assign_view = host_view;
  for (int i = 0; i < 10; ++i) {
    host_anon_stat_view(i) = host_anon_dyn_view(i) = 142;
    host_anon_assign_view(i)                       = 142;
  }
  Kokkos::View<int**, Kokkos::LayoutRight, ExecSpace> d_view("d_view", 100, 10);
#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecSpace, int>(0, 100), KOKKOS_LAMBDA(int i) {
        int* ptr = &(d_view(i, 0));
        Kokkos::View<int[10], Kokkos::AnonymousSpace> d_anon_stat_view(ptr);
        Kokkos::View<int*, Kokkos::AnonymousSpace> d_anon_dyn_view(ptr, 10);
        auto sub = Kokkos::subview(d_view, i, Kokkos::ALL());
        Kokkos::View<int*, Kokkos::AnonymousSpace> d_anon_assign_view = sub;
        for (int j = 0; j < 10; ++j) {
          d_anon_stat_view(j) = 50;
          d_anon_assign_view(j) += 50;
          d_anon_dyn_view(j) += 42;
        }
      });
  Kokkos::fence();
#endif
}

TEST(TEST_CATEGORY, anonymous_space) { test_anonymous_space(); }

template <class ExecSpace>
struct TestViewOverloadResolution {
  // Overload based on value_type and rank
  static int foo(Kokkos::View<const double**, ExecSpace> /*a*/) { return 1; }
  static int foo(Kokkos::View<const int**, ExecSpace> /*a*/) { return 2; }
  static int foo(Kokkos::View<const double***, ExecSpace> /*a*/) { return 3; }

  // Overload based on compile time dimensions
  static int bar(Kokkos::View<double * [3], ExecSpace> /*a*/) { return 4; }
  static int bar(Kokkos::View<double * [4], ExecSpace> /*a*/) { return 5; }

  static void test_function_overload() {
    Kokkos::View<double**, typename ExecSpace::execution_space::array_layout,
                 ExecSpace>
        a("A", 10, 3);
    int data_type_1 = foo(a);
    int data_type_3 =
        foo(Kokkos::View<const double**,
                         typename ExecSpace::execution_space::array_layout,
                         ExecSpace>(a));
    Kokkos::View<double***, typename ExecSpace::execution_space::array_layout,
                 ExecSpace>
        b("B", 10, 3, 4);
    int data_type_2 = foo(b);
    Kokkos::View<double * [3],
                 typename ExecSpace::execution_space::array_layout, ExecSpace>
        c(a);
    int static_extent = bar(c);
    ASSERT_EQ(1, data_type_1);
    ASSERT_EQ(3, data_type_2);
    ASSERT_EQ(1, data_type_3);
    ASSERT_EQ(4, static_extent);
  }
};

TEST(TEST_CATEGORY, view_overload_resolution) {
  TestViewOverloadResolution<TEST_EXECSPACE>::test_function_overload();
}
}  // namespace Test

#include <TestViewIsAssignable.hpp>
