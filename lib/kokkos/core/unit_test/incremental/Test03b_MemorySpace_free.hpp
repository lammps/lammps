// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

/// @Kokkos_Feature_Level_Required:3
// Unit test for Kokkos free.
// We constantly allocate and free the memory.
// If the kokkos_free does not free the allocated memory,
// we will exceed the available space.

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
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
