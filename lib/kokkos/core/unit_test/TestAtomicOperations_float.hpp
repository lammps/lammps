// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestAtomicOperations.hpp>

namespace Test {
TEST(TEST_CATEGORY, atomic_operations_float) {
  const int start = -5;
  const int end   = 11;
  for (int i = start; i < end; ++i) {
    for (int t = 0; t < 7; t++)
      // FIXME_32BIT disable division test for 32bit where we have accuracy
      // issues with division atomics still compile it though
      if (t != 5 || sizeof(void*) == 8) {
        ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestNonIntegralType<
                     double, TEST_EXECSPACE>(i, end - i + start, t)));
      }
  }
}
}  // namespace Test
