// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>
#include <cstdlib>
#include <Kokkos_Core.hpp>

TEST(ProfilingDeathTest, EmptyProfilingLibraryErrors) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  ASSERT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
      },
      "Error: No profiling interface symbols were found in profiling library");
}
