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

#include <PerfTest_ViewResize.hpp>

namespace Test {

TEST(default_exec, ViewResize_Rank8) {
// FIXME_SYCL Avoid running out of resources on the CUDA GPU used in the CI
#ifdef KOKKOS_ENABLE_SYCL
  printf("Resize View Performance for LayoutLeft:\n");
  run_resizeview_tests8<Kokkos::LayoutLeft>(9, 1);
  printf("Resize View Performance for LayoutRight:\n");
  run_resizeview_tests8<Kokkos::LayoutRight>(9, 1);
#else
  printf("Resize View Performance for LayoutLeft:\n");
  run_resizeview_tests8<Kokkos::LayoutLeft>(10, 1);
  printf("Resize View Performance for LayoutRight:\n");
  run_resizeview_tests8<Kokkos::LayoutRight>(10, 1);
#endif
}

}  // namespace Test
