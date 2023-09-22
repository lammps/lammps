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

#include <TestViewAPI.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_api_c) {
  TestViewAPI<double, TEST_EXECSPACE>::run_test_deep_copy_empty();
  TestViewAPI<double, TEST_EXECSPACE>::run_test_view_operator_b();
}

}  // namespace Test
