// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestViewAPI.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_api_c) {
  TestViewAPI<double, TEST_EXECSPACE>::run_test_refcount_exception();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_deep_copy_empty();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_view_operator_b();
}

}  // namespace Test
