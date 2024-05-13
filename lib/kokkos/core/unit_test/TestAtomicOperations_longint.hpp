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
TEST(TEST_CATEGORY, atomic_operations_long) {
  const int start = -5;
  const int end   = 11;
  for (int i = start; i < end; ++i) {
    for (int t = 0; t < 16; t++)
      ASSERT_TRUE((TestAtomicOperations::AtomicOperationsTestIntegralType<
                   long int, TEST_EXECSPACE>(i, end - i + start, t)));
  }
}
}  // namespace Test
