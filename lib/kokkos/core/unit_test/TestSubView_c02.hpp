// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_SUBVIEW_C02_HPP
#define KOKKOS_TEST_SUBVIEW_C02_HPP
#include <TestViewSubview.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_subview_1d_assign_atomic) {
  TestViewSubview::test_1d_assign<TEST_EXECSPACE,
                                  Kokkos::MemoryTraits<Kokkos::Atomic> >();
}

}  // namespace Test
#endif
