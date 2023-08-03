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

/// @Kokkos_Feature_Level_Required:2
// Unit test for atomic exchange, atomic add and atomic sub.
// Atomic exchange test : we interchange value1 with value2 and check for
// correctness. Atomic add test : we add value2 to value1 and check for
// correctness. Atomic sub test : we subtract value2 from value1 and check for
// correctmess.

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

namespace Test {

struct TestIncrAtomic {
  using value_type  = double;
  value_type value1 = 1.5, value2 = 0.5;

  void testExchange() {
    value_type ret_value = Kokkos::atomic_exchange(&value1, value2);

    ASSERT_EQ(value1, 0.5);
    ASSERT_EQ(ret_value, 1.5);
  }

  void testAdd() {
    Kokkos::atomic_add(&value1, value2);

    ASSERT_EQ(value1, 2.0);
  }

  void testSub() {
    Kokkos::atomic_sub(&value1, value2);

    ASSERT_EQ(value1, 1.0);
  }
};

TEST(TEST_CATEGORY, IncrTest_02_AtomicExchange) {
  TestIncrAtomic test;
  test.testExchange();
}

TEST(TEST_CATEGORY, IncrTest_02_AtomicAdd) {
  TestIncrAtomic test;
  test.testAdd();
}

TEST(TEST_CATEGORY, IncrTest_02_AtomicSub) {
  TestIncrAtomic test;
  test.testSub();
}

}  // namespace Test
