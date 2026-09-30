// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestViewAPI.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_api_d) {
  TestViewAPI<double, TEST_EXECSPACE>::run_test_const();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_subview();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_subview_strided();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_vector();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_view_operator_c();
}

}  // namespace Test
