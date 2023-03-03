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

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

/// @Kokkos_Feature_Level_Required:16
// Incremental test for parallel_scan.
// perform scan on a 1D view of double's and check for correctness.

namespace Test {

using value_type = double;
const int N      = 10;

template <class ExecSpace>
struct TestScan {
  // 1D  View of double
  using View_1D = typename Kokkos::View<value_type *, ExecSpace>;

  void parallel_scan() {
    View_1D d_data("data", N);

    // Initialize data.
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecSpace>(0, N),
        KOKKOS_LAMBDA(const int i) { d_data(i) = i * 0.5; });

    // Exclusive parallel_scan call.
    Kokkos::parallel_scan(
        Kokkos::RangePolicy<ExecSpace>(0, N),
        KOKKOS_LAMBDA(const int i, value_type &update_value, const bool final) {
          const value_type val_i = d_data(i);
          if (final) d_data(i) = update_value;

          update_value += val_i;
        });

    // Copy back the data.
    auto h_data =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), d_data);

    // Check Correctness
    ASSERT_EQ(h_data(0), 0.0);
    value_type upd = h_data(0);
    for (int i = 1; i < N; ++i) {
      upd += (i - 1) * 0.5;
      ASSERT_EQ(h_data(i), upd);
    }
  }
};

TEST(TEST_CATEGORY, IncrTest_16_parallelscan) {
  TestScan<TEST_EXECSPACE> test;
  test.parallel_scan();
}

}  // namespace Test
