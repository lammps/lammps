// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SUBVIEW_B_HPP
#define KOKKOS_TEST_SUBVIEW_B_HPP
#include <TestViewSubview.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_subview_layoutleft_to_layoutleft) {
  TestViewSubview::test_layoutleft_to_layoutleft<TEST_EXECSPACE>();
  TestViewSubview::test_layoutleft_to_layoutleft<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::Atomic> >();
  TestViewSubview::test_layoutleft_to_layoutleft<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::RandomAccess> >();
}

TEST(TEST_CATEGORY, view_subview_layoutright_to_layoutright) {
  TestViewSubview::test_layoutright_to_layoutright<TEST_EXECSPACE>();
  TestViewSubview::test_layoutright_to_layoutright<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::Atomic> >();
  TestViewSubview::test_layoutright_to_layoutright<
      TEST_EXECSPACE, Kokkos::MemoryTraits<Kokkos::RandomAccess> >();
}

}  // namespace Test
#endif
