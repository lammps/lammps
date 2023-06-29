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

#include <TestAtomicOperations.hpp>

namespace Test {
TEST(TEST_CATEGORY, atomic_operations_unsigned) {
  const int start = 1;  // Avoid zero for division.
  const int end   = 11;
  for (int i = start; i < end; ++i) {
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 1)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 2)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 3)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 4)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 5)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 6)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 7)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 8)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 9)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 11)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 12)));
    ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                 unsigned int, TEST_EXECSPACE>(start, end - i, 13)));
    ASSERT_TRUE(
        (TestAtomicOperations::AtomicOperationsTestUnsignedIntegralType<
            unsigned int, TEST_EXECSPACE>(start, end - i, 1)));  // Wrapping Inc
    ASSERT_TRUE(
        (TestAtomicOperations::AtomicOperationsTestUnsignedIntegralType<
            unsigned int, TEST_EXECSPACE>(start, end - i, 2)));  // Wrapping Dec
  }
}
}  // namespace Test
