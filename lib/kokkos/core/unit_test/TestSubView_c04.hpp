// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SUBVIEW_C04_HPP
#define KOKKOS_TEST_SUBVIEW_C04_HPP
#include <TestViewSubview.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_subview_2d_from_3d) {
  TestViewSubview::test_2d_subview_3d<TEST_EXECSPACE>();
}

}  // namespace Test
#endif
