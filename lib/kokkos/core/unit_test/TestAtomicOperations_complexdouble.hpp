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

#include <TestAtomicOperations.hpp>

namespace Test {
TEST(TEST_CATEGORY, atomic_operations_complexdouble) {
#if defined(KOKKOS_ENABLE_SYCL) && \
    !defined(KOKKOS_IMPL_SYCL_DEVICE_GLOBAL_SUPPORTED)
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::SYCL>)
    GTEST_SKIP() << "skipping since device_global variables are not available";
#endif
  const int start = 1;  // Avoid zero for division.
  const int end   = 11;
  for (int i = start; i < end; ++i) {
    ASSERT_TRUE(
        (TestAtomicOperations::MulAtomicTest<Kokkos::complex<double>,
                                             TEST_EXECSPACE>(start, end - i)));
    ASSERT_TRUE(
        (TestAtomicOperations::DivAtomicTest<Kokkos::complex<double>,
                                             TEST_EXECSPACE>(start, end - i)));
  }
}
}  // namespace Test
