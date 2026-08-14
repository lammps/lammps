// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

TEST(TEST_CATEGORY_DEATH, view_subview_wrong_extents) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

#ifndef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
  GTEST_SKIP() << "only enforced when debug bound checks is enabled";
#ifdef KOKKOS_COMPILER_NVHPC
  __builtin_unreachable();
#endif
#endif

  TestViewSubview::test_subview_extents<1, TEST_EXECSPACE>();
  TestViewSubview::test_subview_extents<2, TEST_EXECSPACE>();
  TestViewSubview::test_subview_extents<3, TEST_EXECSPACE>();
  TestViewSubview::test_subview_extents<4, TEST_EXECSPACE>();
  TestViewSubview::test_subview_extents<5, TEST_EXECSPACE>();
  TestViewSubview::test_subview_extents<6, TEST_EXECSPACE>();
  TestViewSubview::test_subview_extents<7, TEST_EXECSPACE>();
  TestViewSubview::test_subview_extents<8, TEST_EXECSPACE>();
}

}  // namespace Test
#endif
