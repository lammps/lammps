// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SUBVIEW_C09_HPP
#define KOKKOS_TEST_SUBVIEW_C09_HPP
#include <TestViewSubview.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_subview_3d_from_5d_left_randomaccess) {
  TestViewSubview::test_3d_subview_5d_left<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::RandomAccess> >();
}

}  // namespace Test
#endif
