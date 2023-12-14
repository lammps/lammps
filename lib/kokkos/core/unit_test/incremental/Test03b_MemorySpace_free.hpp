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

/// @Kokkos_Feature_Level_Required:3
// Unit test for Kokkos free.
// We constantly allocate and free the memory.
// If the kokkos_free does not free the allocated memory,
// we will exceed the available space.

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

namespace Test {

using value_type = double;

// Allocate M number of value_type elements N number of times.
const int N = 100000;
const int M = 100000;

template <class ExecSpace>
struct TestIncrMemorySpace_free {
  using memory_space = typename ExecSpace::memory_space;

  void test_free() {
    for (int i = 0; i < N; ++i) {
      auto *data = static_cast<value_type *>(
          Kokkos::kokkos_malloc<memory_space>("data", M * sizeof(value_type)));

      ASSERT_NE(data, nullptr);

      Kokkos::kokkos_free<memory_space>(data);
    }
  }
};

TEST(TEST_CATEGORY, IncrTest_03b_memspace_free) {
  TestIncrMemorySpace_free<TEST_EXECSPACE> test;
  test.test_free();
}

}  // namespace Test
