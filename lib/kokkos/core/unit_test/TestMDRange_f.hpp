// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestMDRange.hpp>

namespace Test {

TEST(TEST_CATEGORY, mdrange_scalar) {
  TestMDRange_ReduceScalar<TEST_EXECSPACE>::test_scalar_reduce(12, 11);
}

}  // namespace Test
