// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SUBVIEW_C06_HPP
#define KOKKOS_TEST_SUBVIEW_C06_HPP
#include <TestViewSubview.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_subview_2d_from_3d_randomaccess) {
  TestViewSubview::test_2d_subview_3d<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::RandomAccess> >();
}

}  // namespace Test
#endif
