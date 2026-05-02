// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SUBVIEW_C14_HPP
#define KOKKOS_TEST_SUBVIEW_C14_HPP
#include <TestViewSubview.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_subview_memory_traits_construction) {
  TestViewSubview::test_subview_memory_traits_construction();
}

}  // namespace Test
#endif
