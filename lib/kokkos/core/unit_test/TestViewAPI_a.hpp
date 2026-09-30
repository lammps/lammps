// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestViewAPI.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_api_a) {
  TestViewAPI<double, TEST_EXECSPACE>::run_test();
}

}  // namespace Test
