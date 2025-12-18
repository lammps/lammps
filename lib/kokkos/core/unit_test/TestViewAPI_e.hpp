// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <sstream>
#include <iostream>

namespace Test {

TEST(TEST_CATEGORY, view_remap) {
  enum { N0 = 3, N1 = 2, N2 = 8, N3 = 9 };

#if defined(KOKKOS_ENABLE_CUDA)
#define EXECSPACE                                                  \
  std::conditional_t<std::is_same_v<TEST_EXECSPACE, Kokkos::Cuda>, \
                     Kokkos::CudaHostPinnedSpace, TEST_EXECSPACE>
#elif defined(KOKKOS_ENABLE_HIP)
#define EXECSPACE                                                 \
  std::conditional_t<std::is_same_v<TEST_EXECSPACE, Kokkos::HIP>, \
                     Kokkos::HIPHostPinnedSpace, TEST_EXECSPACE>
#elif defined(KOKKOS_ENABLE_SYCL)
#define EXECSPACE                                                  \
  std::conditional_t<std::is_same_v<TEST_EXECSPACE, Kokkos::SYCL>, \
                     Kokkos::SYCLHostUSMSpace, TEST_EXECSPACE>
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
#define EXECSPACE Kokkos::HostSpace
#else
#define EXECSPACE TEST_EXECSPACE
#endif

  using output_type =
      Kokkos::View<double* [N1][N2][N3], Kokkos::LayoutRight, EXECSPACE>;

  using input_type =
      Kokkos::View<int** [N2][N3], Kokkos::LayoutLeft, EXECSPACE>;

  using diff_type =
      Kokkos::View<int* [N0][N2][N3], Kokkos::LayoutLeft, EXECSPACE>;

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
  for (size_t i = 0; i < view_type::rank; ++i) {
    ASSERT_EQ(view.stride(i), expected_stride);
    ASSERT_EQ(all_strides[i], expected_stride);
    expected_stride *= view.extent(i);
  }
  ASSERT_EQ(all_strides[view_type::rank()], expected_stride);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  for (size_t i = view_type::rank(); i < size_t(8); ++i) {
    ASSERT_EQ(view.stride(i), view.stride(view_type::rank() - 1) *
                                  view.extent(view_type::rank() - 1));
  }
#endif
}

template <typename DataType, typename... Extents>
void test_right_stride(Extents... extents) {
  using view_type =
      Kokkos::View<DataType, Kokkos::LayoutRight, Kokkos::HostSpace>;
  view_type view("view", extents...);
  size_t expected_stride = 1;
  size_t all_strides[view_type::rank + 1];
  view.stride(all_strides);
  for (size_t ri = 0; ri < view_type::rank; ++ri) {
    auto i = view_type::rank - 1 - ri;
    ASSERT_EQ(view.stride(i), expected_stride);
    ASSERT_EQ(all_strides[i], expected_stride);
    expected_stride *= view.extent(i);
  }
  ASSERT_EQ(all_strides[view_type::rank()], expected_stride);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  for (size_t i = view_type::rank(); i < size_t(8); ++i) {
    ASSERT_EQ(view.stride(i), size_t(1));
  }
#endif
}

template <typename DataType, typename... Extents>
void test_stride_stride(Extents... extents) {
  using view_type =
      Kokkos::View<DataType, Kokkos::LayoutStride, Kokkos::HostSpace>;

  auto test = [](auto view_org) {
    view_type view(view_org);

    size_t all_strides[view_type::rank() + 1];
    view.stride(all_strides);

    size_t max_stride     = 0;
    size_t max_stride_idx = 0;
    for (size_t i = 0; i < view_type::rank(); ++i) {
      ASSERT_EQ(view.stride(i), view_org.stride(i));
      ASSERT_EQ(all_strides[i], view_org.stride(i));
      if (view.stride(i) > max_stride) {
        max_stride     = view.stride(i);
        max_stride_idx = i;
      }
    }
    ASSERT_EQ(all_strides[view_type::rank()],
              max_stride * view.extent(max_stride_idx));
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
    for (size_t i = view_type::rank(); i < size_t(8); ++i) {
      ASSERT_EQ(view.stride(i), size_t(0));
    }
#endif
  };

  Kokkos::View<DataType, Kokkos::LayoutRight, Kokkos::HostSpace> view_right(
      "view", extents...);
  test(view_right);
  Kokkos::View<DataType, Kokkos::LayoutRight, Kokkos::HostSpace> view_left(
      "view", extents...);
  test(view_left);
}

template <typename DataType, typename... Extents>
void test_stride(Extents... extents) {
  test_right_stride<DataType>(extents...);
  test_left_stride<DataType>(extents...);
  test_stride_stride<DataType>(extents...);
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
}

TEST(TEST_CATEGORY, anonymous_space) { test_anonymous_space(); }

template <class ExecSpace>
struct TestViewOverloadResolution {
  // Overload based on value_type and rank
  static int foo(Kokkos::View<const double**, ExecSpace> /*a*/) { return 1; }
  static int foo(Kokkos::View<const int**, ExecSpace> /*a*/) { return 2; }
  static int foo(Kokkos::View<const double***, ExecSpace> /*a*/) { return 3; }

  // Overload based on compile time dimensions
  static int bar(Kokkos::View<double* [3], ExecSpace> /*a*/) { return 4; }
  static int bar(Kokkos::View<double* [4], ExecSpace> /*a*/) { return 5; }

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
    Kokkos::View<double* [3], typename ExecSpace::execution_space::array_layout,
                 ExecSpace>
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

template <typename MemorySpace>
struct TestViewAllocationLargeRank {
  using ViewType = Kokkos::View<char********, MemorySpace>;

  KOKKOS_FUNCTION void operator()(int) const {
    size_t idx = v.extent(0) - 1;
    auto& lhs  = v(idx, idx, idx, idx, idx, idx, idx, idx);
    lhs        = 42;  // This is where it segfaulted
  }

  ViewType v;
};

// Orin and other Jetson devices come with smaller memory capacity, 8GB total
TEST(TEST_CATEGORY, view_allocation_large_rank) {
#ifdef KOKKOS_ARCH_AMPERE87
  GTEST_SKIP() << "skipping for Jetson devices that have only 8GB memory";
#endif
#ifdef KOKKOS_IMPL_32BIT
  GTEST_SKIP() << "skipping for 32-bit builds";
#endif
  using ExecutionSpace = typename TEST_EXECSPACE::execution_space;
  using MemorySpace    = typename TEST_EXECSPACE::memory_space;
  constexpr int dim    = 15;
  using FunctorType    = TestViewAllocationLargeRank<MemorySpace>;
  typename FunctorType::ViewType v("v", dim, dim, dim, dim, dim, dim, dim, dim);

  Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0, 1),
                       FunctorType{v});
  typename FunctorType::ViewType v_single(v.data() + v.size() - 1, 1, 1, 1, 1,
                                          1, 1, 1, 1);
  auto result =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, v_single);
  ASSERT_EQ(result(0, 0, 0, 0, 0, 0, 0, 0), 42);
}

template <typename ExecSpace, typename ViewType>
struct TestViewShmemSizeOnDevice {
  using ViewTestType = Kokkos::View<size_t, ExecSpace>;

  TestViewShmemSizeOnDevice(size_t d1_, size_t d2_, size_t d3_)
      : d1(d1_), d2(d2_), d3(d3_), shmemSize("shmemSize") {}

  KOKKOS_FUNCTION void operator()(const int&) const {
    auto shmem  = ViewType::shmem_size(d1, d2, d3);
    shmemSize() = shmem;
  }

  size_t d1, d2, d3;
  ViewTestType shmemSize;
};

TEST(TEST_CATEGORY, view_shmem_size_on_device) {
  using ExecSpace = typename TEST_EXECSPACE::execution_space;
  using ViewType  = Kokkos::View<int64_t***, ExecSpace>;

  constexpr size_t d1 = 5;
  constexpr size_t d2 = 7;
  constexpr size_t d3 = 11;

  TestViewShmemSizeOnDevice<ExecSpace, ViewType> testShmemSize(d1, d2, d3);

  Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace>(0, 1), testShmemSize);

  auto size = ViewType::shmem_size(d1, d2, d3);

  auto shmemSizeHost = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace(), testShmemSize.shmemSize);

  ASSERT_EQ(size, shmemSizeHost());
}

}  // namespace Test

#include <TestViewIsAssignable.hpp>
