//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

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
