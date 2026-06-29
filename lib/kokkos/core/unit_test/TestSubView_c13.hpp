// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SUBVIEW_C13_HPP
#define KOKKOS_TEST_SUBVIEW_C13_HPP
#include <TestViewSubview.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_test_unmanaged_subview_reset) {
  TestViewSubview::test_unmanaged_subview_reset<TEST_EXECSPACE>();
}

}  // namespace Test
#endif
