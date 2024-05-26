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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

template <class ExecutionSpace>
void test_kokkos_printf() {
  ::testing::internal::CaptureStdout();
  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(0, 1),
      KOKKOS_LAMBDA(int) { Kokkos::printf("Print an integer: %d", 2); });
  Kokkos::fence();
  auto const captured = ::testing::internal::GetCapturedStdout();
  std::string expected_string("Print an integer: 2");
  ASSERT_EQ(captured, expected_string);
}

// FIXME_OPENMPTARGET non-string-literal argument used in printf is not
// supported for spir64
#if !(defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(KOKKOS_ARCH_INTEL_GPU))
TEST(TEST_CATEGORY, kokkos_printf) { test_kokkos_printf<TEST_EXECSPACE>(); }
#endif
