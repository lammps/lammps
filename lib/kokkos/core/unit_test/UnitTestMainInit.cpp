// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

int main(int argc, char *argv[]) {
  Kokkos::initialize(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);

  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
