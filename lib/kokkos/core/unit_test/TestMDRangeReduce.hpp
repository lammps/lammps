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

namespace {

template <typename T>
void MDRangeReduceTester([[maybe_unused]] int bound, int k) {
  const auto policy_MD = Kokkos::MDRangePolicy<Kokkos::Rank<2>, TEST_EXECSPACE>(
      {0, 0}, {bound, 2});

  // No explicit fence() calls needed because result is in HostSpace
  {
    T lor_MD = 0;
    Kokkos::parallel_reduce(
        policy_MD,
        KOKKOS_LAMBDA(const int i, const int, T& res) { res = res || i == k; },
        Kokkos::LOr<T>(lor_MD));
    EXPECT_EQ(lor_MD, 1);
  }
  {
    // Stick just a few true values in the Logical-OR reduction space,
    // to try to make sure every value is being captured
    T land_MD = 0;
    Kokkos::parallel_reduce(
        policy_MD, KOKKOS_LAMBDA(const int, const int, T& res) { res = 1; },
        Kokkos::LAnd<T>(land_MD));
    EXPECT_EQ(land_MD, 1);
  }
}

TEST(TEST_CATEGORY, mdrange_parallel_reduce_primitive_types) {
#if defined(KOKKOS_ENABLE_OPENMPTARGET)
  GTEST_SKIP() << "FIXME OPENMPTARGET Tests of MDRange reduce over values "
                  "smaller than int would fail";
#elif defined(KOKKOS_ENABLE_CUDA) && !defined(KOKKOS_ENABLE_CUDA_LAMBDA)
  GTEST_SKIP() << "Skipped ENABLE_CUDA_LAMBDA";
#else
  for (int bound : {0, 1, 7, 32, 65, 7000}) {
    for (int k = 0; k < bound; ++k) {
      MDRangeReduceTester<bool>(bound, k);
      MDRangeReduceTester<signed char>(bound, k);
      MDRangeReduceTester<int8_t>(bound, k);
      MDRangeReduceTester<int16_t>(bound, k);
      MDRangeReduceTester<int32_t>(bound, k);
      MDRangeReduceTester<int64_t>(bound, k);
    }
  }
#endif
}

}  // namespace
