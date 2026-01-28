// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestAtomicOperations.hpp>

using namespace TestAtomicOperations;

namespace Test {
TEST(TEST_CATEGORY, atomic_operations_complexfloat) {
  // FIXME_OPENMPTARGET - causes runtime failure with CrayClang compiler
#if defined(KOKKOS_COMPILER_CRAY_LLVM) && defined(KOKKOS_ENABLE_OPENMPTARGET)
  GTEST_SKIP() << "known to fail with OpenMPTarget+Cray LLVM";
#endif
  const int start = -5;
  const int end   = 11;
  for (int i = start; i < end; ++i) {
    using T   = Kokkos::complex<float>;
    T old_val = static_cast<T>(i);
    T update  = static_cast<T>(end - i - start);
    ASSERT_TRUE(
        (atomic_op_test<AddAtomicTest, T, TEST_EXECSPACE>(old_val, update)));
    ASSERT_TRUE(
        (atomic_op_test<SubAtomicTest, T, TEST_EXECSPACE>(old_val, update)));
    ASSERT_TRUE(
        (atomic_op_test<MulAtomicTest, T, TEST_EXECSPACE>(old_val, update)));

    // FIXME_32BIT disable division test for 32bit where we have accuracy issues
    // with division atomics still compile it though
    if (sizeof(void*) == 8) {
      ASSERT_TRUE((update != 0
                       ? atomic_op_test<DivAtomicTest, T, TEST_EXECSPACE>(
                             old_val, update)
                       : true));
    }
    ASSERT_TRUE((atomic_op_test<LoadStoreAtomicTest, T, TEST_EXECSPACE>(
        old_val, update)));
  }
}
}  // namespace Test
