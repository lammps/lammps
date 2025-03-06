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

#include <gtest/gtest.h>
#include <impl/Kokkos_CheckedIntegerOps.hpp>
#include <limits>

namespace {

TEST(TEST_CATEGORY, checked_integer_operations_multiply_overflow) {
  {
    auto result      = 1u;
    auto is_overflow = Kokkos::Impl::multiply_overflow(1u, 2u, result);
    EXPECT_EQ(result, 2u);
    EXPECT_FALSE(is_overflow);
  }
  {
    auto result      = 1u;
    auto is_overflow = Kokkos::Impl::multiply_overflow(
        std::numeric_limits<unsigned>::max(), 2u, result);
    EXPECT_TRUE(is_overflow);
  }
}

TEST(TEST_CATEGORY_DEATH, checked_integer_operations_multiply_overflow_abort) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  {
    auto result = Kokkos::Impl::multiply_overflow_abort(1u, 2u);
    EXPECT_EQ(result, 2u);
  }
  {
    ASSERT_DEATH(Kokkos::Impl::multiply_overflow_abort(
                     std::numeric_limits<unsigned>::max(), 2u),
                 "Arithmetic overflow detected.");
  }
}

}  // namespace
