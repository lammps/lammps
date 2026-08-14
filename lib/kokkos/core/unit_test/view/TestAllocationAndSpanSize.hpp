// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif

#include <gtest/gtest.h>

template <class ViewT, std::integral... Sizes>
void test_required_span_size_single_rank(size_t expected_size,
                                         std::string label, Sizes... sizes) {
  ViewT view(label, sizes...);
#ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  // Lets get the required size two ways: based on mapping + accessor and using
  // the required_allocation_size
  size_t span_size = view.mapping().required_span_size();
  size_t span_size_based_bytes =
      span_size * sizeof(typename ViewT::element_type);
  ASSERT_EQ(span_size_based_bytes, expected_size);
#endif
  size_t req_allocation_size = ViewT::required_allocation_size(sizes...);
  ASSERT_EQ(req_allocation_size, expected_size);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_5
  if constexpr (sizeof...(Sizes) == ViewT::rank()) {
    size_t req_allocation_size_extra =
        ViewT::required_allocation_size(sizes..., 3);
    ASSERT_EQ(req_allocation_size_extra, expected_size);
  }
#endif
}

template <class Layout>
void test_required_span_size_layout() {
  // 8 is the size of the element types we store.
  // While Foo::Bar uses float, the accessor and reference type associated
  // with Foo::Bar actually make the allocation store double!
  test_required_span_size_single_rank<
      Kokkos::View<double, Layout, TEST_EXECSPACE>>(8, "A");
  test_required_span_size_single_rank<
      Kokkos::View<double*, Layout, TEST_EXECSPACE>>(7 * 8, "A", 7);
  test_required_span_size_single_rank<
      Kokkos::View<double[7], Layout, TEST_EXECSPACE>>(7 * 8, "A", 7);
  test_required_span_size_single_rank<
      Kokkos::View<double[7], Layout, TEST_EXECSPACE>>(7 * 8, "A");
  test_required_span_size_single_rank<
      Kokkos::View<double**, Layout, TEST_EXECSPACE>>(7 * 11 * 8, "A", 7, 11);
  test_required_span_size_single_rank<
      Kokkos::View<double* [11], Layout, TEST_EXECSPACE>>(7 * 11 * 8, "A", 7,
                                                          11);
  test_required_span_size_single_rank<
      Kokkos::View<double* [11], Layout, TEST_EXECSPACE>>(7 * 11 * 8, "A", 7);
  test_required_span_size_single_rank<
      Kokkos::View<double*******, Layout, TEST_EXECSPACE>>(
      7 * 11 * 13 * 17 * 19 * 2 * 3 * 8, "A", 7, 11, 13, 17, 19, 2, 3);
  test_required_span_size_single_rank<
      Kokkos::View<double***** [2][3], Layout, TEST_EXECSPACE>>(
      7 * 11 * 13 * 17 * 19 * 2 * 3 * 8, "A", 7, 11, 13, 17, 19, 2, 3);
  test_required_span_size_single_rank<
      Kokkos::View<double***** [2][3], Layout, TEST_EXECSPACE>>(
      7 * 11 * 13 * 17 * 19 * 2 * 3 * 8, "A", 7, 11, 13, 17, 19);
}

TEST(TEST_CATEGORY, view_required_span_size) {
  test_required_span_size_layout<Kokkos::LayoutLeft>();
  test_required_span_size_layout<Kokkos::LayoutRight>();
}
