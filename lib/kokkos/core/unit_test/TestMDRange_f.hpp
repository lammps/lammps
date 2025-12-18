// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestMDRange.hpp>

namespace Test {

// FIXME_OPENMPTARGET requires MDRange parallel_reduce
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, mdrange_scalar) {
  TestMDRange_ReduceScalar<TEST_EXECSPACE>::test_scalar_reduce(12, 11);
}
#endif

}  // namespace Test
