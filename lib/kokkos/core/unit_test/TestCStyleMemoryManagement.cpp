// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include <TestDefaultDeviceType_Category.hpp>

#include <gtest/gtest.h>

namespace {

TEST(defaultdevicetype, c_style_memory_management_malloc_realloc_and_free) {
  int* data = static_cast<int*>(Kokkos::kokkos_malloc(100 * sizeof(int)));
  ASSERT_NE(data, nullptr);

  data = static_cast<int*>(Kokkos::kokkos_realloc(data, 120 * sizeof(int)));
  ASSERT_NE(data, nullptr);

  Kokkos::kokkos_free(data);
}

TEST(defaultdevicetype, c_style_memory_management_malloc_zero_byte_and_free) {
  int* data2 = static_cast<int*>(Kokkos::kokkos_malloc(0));
  ASSERT_EQ(data2, nullptr);

  Kokkos::kokkos_free(data2);
}

}  // namespace
