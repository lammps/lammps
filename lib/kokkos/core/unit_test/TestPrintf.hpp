// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

template <class ExecutionSpace>
void test_kokkos_printf() {
  ::testing::internal::CaptureStdout();
  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int) { Kokkos::printf("Print an integer: %d", 2); });
  Kokkos::fence();
  auto const captured = ::testing::internal::GetCapturedStdout();
  std::string expected_string("Print an integer: 2");
  ASSERT_EQ(captured, expected_string);
}

TEST(TEST_CATEGORY, kokkos_printf) { test_kokkos_printf<TEST_EXECSPACE>(); }
