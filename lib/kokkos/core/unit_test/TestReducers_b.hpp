// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <TestReducers.hpp>

// FIXME_OPENMPTARGET - Fails at runtime post clang/16
#if defined(KOKKOS_ENABLE_OPENMPTARGER) && defined(KOKKOS_COMPILER_CLANG) && \
    (KOKKOS_COMPILER_CLANG >= 1600)
namespace Test {
TEST(TEST_CATEGORY, reducers_size_t) {
  TestReducers<size_t, TEST_EXECSPACE>::execute_integer();
}
}  // namespace Test
#endif
