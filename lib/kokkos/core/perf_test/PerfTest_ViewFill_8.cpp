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

#include <PerfTest_ViewFill.hpp>

namespace Test {
TEST(default_exec, ViewFill_Rank8) {
  printf("ViewFill Performance for LayoutLeft:\n");
  run_fillview_tests8<Kokkos::LayoutLeft>(10, 1);
  printf("ViewFill Performance for LayoutRight:\n");
  run_fillview_tests8<Kokkos::LayoutRight>(10, 1);
}
}  // namespace Test
