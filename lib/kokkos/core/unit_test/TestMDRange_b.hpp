// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestMDRange.hpp>

namespace Test {

TEST(TEST_CATEGORY, mdrange_6d) {
#ifndef KOKKOS_ENABLE_OPENMPTARGET
  // FIXME_OPENMPTARGET requires MDRange parallel_reduce
  TestMDRange_6D<TEST_EXECSPACE>::test_reduce6(100, 10, 10, 10, 5, 5);
#endif
  TestMDRange_6D<TEST_EXECSPACE>::test_for6(10, 10, 10, 10, 5, 5);
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL)
  const int size_x = 2 << 19;  // 2^20
  TestMDRange_6D<TEST_EXECSPACE>::test_for6_eval_once(size_x, 1, 1, 1, 1, 1);
  TestMDRange_6D<TEST_EXECSPACE>::test_for6_eval_once(1, size_x, 1, 1, 1, 1);
  TestMDRange_6D<TEST_EXECSPACE>::test_for6_eval_once(1, 1, size_x, 1, 1, 1);
  TestMDRange_6D<TEST_EXECSPACE>::test_for6_eval_once(1, 1, 1, size_x, 1, 1);
  TestMDRange_6D<TEST_EXECSPACE>::test_for6_eval_once(1, 1, 1, 1, size_x, 1);
  TestMDRange_6D<TEST_EXECSPACE>::test_for6_eval_once(1, 1, 1, 1, 1, size_x);
#endif
}

}  // namespace Test
