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

#include <TestViewAPI.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_layout_left_with_stride) {
  Kokkos::LayoutLeft ll(10, 20);
  ll.stride = 15;
  Kokkos::View<int**, Kokkos::LayoutLeft> a("A", ll);
  ASSERT_EQ(static_cast<int>(a.extent(0)), 10);
  ASSERT_EQ(static_cast<int>(a.extent(1)), 20);
  ASSERT_EQ(static_cast<int>(a.stride(0)), 1);
  ASSERT_EQ(static_cast<int>(a.stride(1)), 15);

  auto ll2 = a.layout();
  ASSERT_EQ(static_cast<int>(ll2.dimension[0]), 10);
  ASSERT_EQ(static_cast<int>(ll2.dimension[1]), 20);
  ASSERT_EQ(static_cast<int>(ll2.stride), 15);
}

TEST(TEST_CATEGORY, view_layout_right_with_stride) {
  Kokkos::LayoutRight lr(10, 20);
  lr.stride = 25;
  Kokkos::View<int**, Kokkos::LayoutRight> a("A", lr);
  ASSERT_EQ(static_cast<int>(a.extent(0)), 10);
  ASSERT_EQ(static_cast<int>(a.extent(1)), 20);
  ASSERT_EQ(static_cast<int>(a.stride(0)), 25);
  ASSERT_EQ(static_cast<int>(a.stride(1)), 1);

  auto lr2 = a.layout();
  ASSERT_EQ(static_cast<int>(lr2.dimension[0]), 10);
  ASSERT_EQ(static_cast<int>(lr2.dimension[1]), 20);
  ASSERT_EQ(static_cast<int>(lr2.stride), 25);
}

TEST(TEST_CATEGORY, view_api_b) {
  TestViewAPI<double, TEST_EXECSPACE>::run_test_view_operator_a();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_mirror();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_scalar();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_contruction_from_layout();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_contruction_from_layout_2();
}

}  // namespace Test
