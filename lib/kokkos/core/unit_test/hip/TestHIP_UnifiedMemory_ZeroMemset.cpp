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

#include <TestHIP_Category.hpp>
#include <Kokkos_Core.hpp>

namespace Test {

// On MI300a with ROCM <= 6.2.0, hipMemsetAsync was failing with an error when
// called on host-allocated buffers. The fix was in PR 7380 to use a
// parallel_for to zero memory
TEST(hip, unified_memory_zero_memset) {
#if !defined(KOKKOS_IMPL_HIP_UNIFIED_MEMORY)
  GTEST_SKIP()
      << "this test should only be run with HIP unified memory enabled";
#endif

  constexpr size_t N = 1024 * 1024;  // size doesn't matter
  std::vector<int> v(N, 1);          // initialize to non-zero
  Kokkos::View<int*, Kokkos::HIPSpace> a(v.data(), N);

  // zero with deep_copy (this is where the error occurs)
  Kokkos::deep_copy(a, 0);

  // see if it was zeroed
  int err;
  Kokkos::parallel_reduce(
      N, KOKKOS_LAMBDA(int i, int& lerr) { lerr += (a[i] != 0); }, err);
  EXPECT_EQ(err, 0);
}
}  // namespace Test
