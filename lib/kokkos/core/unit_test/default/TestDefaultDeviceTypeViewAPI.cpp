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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <TestDefaultDeviceType_Category.hpp>

template <size_t... Ds>
using _sizes = std::integer_sequence<size_t, Ds...>;

template <class>
struct TestViewAPI;
template <class DataType, class Layout, size_t... DynamicSizes,
          size_t... AllSizes>
struct TestViewAPI<
    std::tuple<DataType, Layout, std::integer_sequence<size_t, DynamicSizes...>,
               std::integer_sequence<size_t, AllSizes...>>>
    : public ::testing::Test {
  using data_type   = DataType;
  using layout_type = Layout;
  using space_type  = Kokkos::DefaultExecutionSpace;
  using traits_type =
      Kokkos::MemoryTraits<0>;  // maybe we want to add that later to the matrix
  using view_type =
      Kokkos::View<data_type, layout_type, space_type, traits_type>;
  using alloc_layout_type =
      std::conditional_t<std::is_same<layout_type, Kokkos::LayoutStride>::value,
                         Kokkos::LayoutLeft, layout_type>;
  using d_alloc_type = Kokkos::View<data_type, alloc_layout_type, space_type>;
  using h_alloc_type = typename Kokkos::View<data_type, alloc_layout_type,
                                             space_type>::HostMirror;

  // add a +1 to avoid zero length static array
  size_t dyn_sizes[sizeof...(DynamicSizes) + 1] = {DynamicSizes..., 1};
  size_t all_sizes[sizeof...(AllSizes) + 1]     = {AllSizes..., 1};

  constexpr static size_t expected_rank = sizeof...(AllSizes);

  inline view_type create_view() const {
    return d_alloc_type("TestViewAPI", DynamicSizes...);
  }
};

using Kokkos::LayoutLeft;
using Kokkos::LayoutRight;
using Kokkos::LayoutStride;

using compatible_extents_test_types = ::testing::Types<
    // LayoutLeft
    std::tuple<int, LayoutLeft, _sizes<>, _sizes<>>,
    std::tuple<int[5], LayoutLeft, _sizes<>, _sizes<5>>,
    std::tuple<int*, LayoutLeft, _sizes<5>, _sizes<5>>,
    std::tuple<int[5][10], LayoutLeft, _sizes<>, _sizes<5, 10>>,
    std::tuple<int * [10], LayoutLeft, _sizes<5>, _sizes<5, 10>>,
    std::tuple<int**, LayoutLeft, _sizes<5, 10>, _sizes<5, 10>>,
    std::tuple<int[5][10][15], LayoutLeft, _sizes<>, _sizes<5, 10, 15>>,
    std::tuple<int * [10][15], LayoutLeft, _sizes<5>, _sizes<5, 10, 15>>,
    std::tuple<int* * [15], LayoutLeft, _sizes<5, 10>, _sizes<5, 10, 15>>,
    std::tuple<int***, LayoutLeft, _sizes<5, 10, 15>, _sizes<5, 10, 15>>,
    // LayoutRight
    std::tuple<int, LayoutRight, _sizes<>, _sizes<>>,
    std::tuple<int[5], LayoutRight, _sizes<>, _sizes<5>>,
    std::tuple<int*, LayoutRight, _sizes<5>, _sizes<5>>,
    std::tuple<int[5][10], LayoutRight, _sizes<>, _sizes<5, 10>>,
    std::tuple<int * [10], LayoutRight, _sizes<5>, _sizes<5, 10>>,
    std::tuple<int**, LayoutRight, _sizes<5, 10>, _sizes<5, 10>>,
    std::tuple<int[5][10][15], LayoutRight, _sizes<>, _sizes<5, 10, 15>>,
    std::tuple<int * [10][15], LayoutRight, _sizes<5>, _sizes<5, 10, 15>>,
    std::tuple<int* * [15], LayoutRight, _sizes<5, 10>, _sizes<5, 10, 15>>,
    std::tuple<int***, LayoutRight, _sizes<5, 10, 15>, _sizes<5, 10, 15>>,
    // LayoutStride
    std::tuple<int, LayoutStride, _sizes<>, _sizes<>>,
    std::tuple<int[5], LayoutStride, _sizes<>, _sizes<5>>,
    std::tuple<int*, LayoutStride, _sizes<5>, _sizes<5>>,
    std::tuple<int[5][10], LayoutStride, _sizes<>, _sizes<5, 10>>,
    std::tuple<int * [10], LayoutStride, _sizes<5>, _sizes<5, 10>>,
    std::tuple<int**, LayoutStride, _sizes<5, 10>, _sizes<5, 10>>,
    std::tuple<int[5][10][15], LayoutStride, _sizes<>, _sizes<5, 10, 15>>,
    std::tuple<int * [10][15], LayoutStride, _sizes<5>, _sizes<5, 10, 15>>,
    std::tuple<int* * [15], LayoutStride, _sizes<5, 10>, _sizes<5, 10, 15>>,
    std::tuple<int***, LayoutStride, _sizes<5, 10, 15>, _sizes<5, 10, 15>>,
    // Degenerated Sizes
    std::tuple<int*, LayoutLeft, _sizes<0>, _sizes<0>>,
    std::tuple<int * [10], LayoutLeft, _sizes<0>, _sizes<0, 10>>,
    std::tuple<int* * [15], LayoutLeft, _sizes<0, 0>, _sizes<0, 0, 15>>,
    std::tuple<int*, LayoutRight, _sizes<0>, _sizes<0>>,
    std::tuple<int * [10], LayoutRight, _sizes<0>, _sizes<0, 10>>,
    std::tuple<int* * [15], LayoutRight, _sizes<0, 0>, _sizes<0, 0, 15>>,
    std::tuple<int*, LayoutStride, _sizes<0>, _sizes<0>>,
    std::tuple<int * [10], LayoutStride, _sizes<0>, _sizes<0, 10>>,
    std::tuple<int* * [15], LayoutStride, _sizes<0, 0>, _sizes<0, 0, 15>>>;

TYPED_TEST_SUITE(TestViewAPI, compatible_extents_test_types, );

TYPED_TEST(TestViewAPI, sizes) {
  using view_t = typename TestFixture::view_type;
  auto a       = this->create_view();
  static_assert(view_t::rank == TestFixture::expected_rank,
                "TestViewAPI: Error: rank mismatch");
  size_t expected_span = 1;
  // avoid pointless comparison of unsigned integer with zero warning
  if constexpr (view_t::rank > 0) {
    for (size_t r = 0; r < view_t::rank; r++)
      expected_span *= this->all_sizes[r];
  }

  EXPECT_EQ(expected_span, a.span());
  if constexpr (view_t::rank > 0) {
    for (size_t r = 0; r < view_t::rank; r++) {
      EXPECT_EQ(this->all_sizes[r], a.extent(r));
      EXPECT_EQ(this->all_sizes[r], size_t(a.extent_int(r)));
    }
  }
}
