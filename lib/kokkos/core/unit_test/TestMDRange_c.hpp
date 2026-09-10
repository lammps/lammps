// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestMDRange.hpp>

namespace Test {

TEST(TEST_CATEGORY, mdrange_1d) {
  TestMDRange_1D<TEST_EXECSPACE>::test_construct_policies(127);
  TestMDRange_1D<TEST_EXECSPACE>::test_for1(127);
  TestMDRange_1D<TEST_EXECSPACE>::test_reduce1(127);
}

TEST(TEST_CATEGORY, mdrange_2d) {
  TestMDRange_2D<TEST_EXECSPACE>::test_reduce2(100, 100);
  TestMDRange_2D<TEST_EXECSPACE>::test_for2(100, 100);
}

TEST(TEST_CATEGORY, mdrange_array_reduce) {
  TestMDRange_ReduceArray_1D<TEST_EXECSPACE>::test_arrayreduce1(7);
  TestMDRange_ReduceArray_2D<TEST_EXECSPACE>::test_arrayreduce2(4, 5);
  TestMDRange_ReduceArray_3D<TEST_EXECSPACE>::test_arrayreduce3(4, 5, 10);
}

}  // namespace Test
