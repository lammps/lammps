// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestReducers.hpp>

namespace Test {
TEST(TEST_CATEGORY, reducers_int) {
  TestReducers<int, TEST_EXECSPACE>::execute_integer();
}

TEST(TEST_CATEGORY, reducers_unsigned_int) {
  TestReducers<unsigned int, TEST_EXECSPACE>::execute_integer();
}

}  // namespace Test
