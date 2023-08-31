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

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

/// @Kokkos_Feature_Level_Required:3
// Unit Test for Kokkos malloc.
// Allocate memory to a pointer and check if the allocation has not returned a
// null pointer.

namespace Test {

using value_type       = double;
const int num_elements = 10;

template <class ExecSpace>
struct TestIncrMemorySpace_malloc {
  using memory_space = typename ExecSpace::memory_space;

  void test_malloc() {
    // Allocate memory
    auto *data = static_cast<value_type *>(Kokkos::kokkos_malloc<memory_space>(
        "data", num_elements * sizeof(value_type)));

    // Check if the allocated memory has not returned a NULL
    ASSERT_NE(data, nullptr);

    // Free the allocated memory
    Kokkos::kokkos_free<memory_space>(data);
  }
};

TEST(TEST_CATEGORY, IncrTest_03a_memspace_malloc) {
  TestIncrMemorySpace_malloc<TEST_EXECSPACE> test;
  test.test_malloc();
}

}  // namespace Test
