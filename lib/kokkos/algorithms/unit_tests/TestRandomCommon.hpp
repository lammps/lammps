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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TESTRANDOM_COMMON_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TESTRANDOM_COMMON_HPP

#include <TestRandom.hpp>

namespace Test {

TEST(TEST_CATEGORY, Random_XorShift64) {
  test_random_xorshift64<TEST_EXECSPACE>();
}
TEST(TEST_CATEGORY, Random_XorShift1024_0) {
  test_random_xorshift1024<TEST_EXECSPACE>();
}
}  // namespace Test

#endif
