// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestReducers.hpp>

namespace Test {
TEST(TEST_CATEGORY, reducers_double) {
#ifdef KOKKOS_ENABLE_OPENMPTARGET  // FIXME_OPENMPTARGET
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>) {
    GTEST_SKIP()
        << "skipping since this leads to illegal memory access on device";
  }
#endif
  TestReducers<double, TEST_EXECSPACE>::execute_float();
}
}  // namespace Test
