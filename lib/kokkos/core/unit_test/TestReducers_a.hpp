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

namespace Test {
TEST(TEST_CATEGORY, reducers_int) {
#ifdef KOKKOS_ENABLE_OPENMPTARGET  // FIXME_OPENMPTARGET Failing with clang 17
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "Known to fail with OpenMPTarget and clang 17";
#endif
  TestReducers<int, TEST_EXECSPACE>::execute_integer();
}

}  // namespace Test
