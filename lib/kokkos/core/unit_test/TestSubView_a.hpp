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

#ifndef KOKKOS_TEST_SUBVIEW_A_HPP
#define KOKKOS_TEST_SUBVIEW_A_HPP
#include <TestViewSubview.hpp>

namespace Test {

#ifndef KOKKOS_ENABLE_SYCL
TEST(TEST_CATEGORY, view_subview_auto_1d_left) {
  TestViewSubview::test_auto_1d<Kokkos::LayoutLeft, TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_subview_auto_1d_right) {
  TestViewSubview::test_auto_1d<Kokkos::LayoutRight, TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_subview_auto_1d_stride) {
  TestViewSubview::test_auto_1d<Kokkos::LayoutStride, TEST_EXECSPACE>();
}
#endif

TEST(TEST_CATEGORY, view_subview_assign_strided) {
  TestViewSubview::test_1d_strided_assignment<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_subview_left_0) {
  TestViewSubview::test_left_0<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_subview_left_1) {
  TestViewSubview::test_left_1<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_subview_left_2) {
  TestViewSubview::test_left_2<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_subview_left_3) {
  TestViewSubview::test_left_3<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_subview_right_0) {
  TestViewSubview::test_right_0<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_subview_right_1) {
  TestViewSubview::test_right_1<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_subview_right_3) {
  TestViewSubview::test_right_3<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, view_static_tests) {
  TestViewSubview::TestSubviewStaticSizes<TEST_EXECSPACE,
                                          Kokkos::LayoutLeft>()();
  TestViewSubview::TestSubviewStaticSizes<TEST_EXECSPACE,
                                          Kokkos::LayoutRight>()();
  TestViewSubview::TestExtentsStaticTests<TEST_EXECSPACE>();
}

}  // namespace Test
#endif
