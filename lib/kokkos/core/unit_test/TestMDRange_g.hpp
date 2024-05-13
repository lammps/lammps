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

//#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace Test {

template <typename View>
struct SumView {
  const View m_view;
  KOKKOS_FUNCTION void operator()(const int i, const int j, int& update) const {
    update += m_view(i, j);
  }

  SumView(View view) : m_view(view) {}

  int run() {
    int sum_view = 0;
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<typename View::execution_space, Kokkos::Rank<2>>(
            {0, 0}, {m_view.extent(0), m_view.extent(1)}),
        *this, sum_view);
    return sum_view;
  }
};

template <typename ExecutionSpace>
struct TestMDRangeLargeDeepCopy {
  static void run() {
    ExecutionSpace exec;
    using MemorySpace = typename ExecutionSpace::memory_space;
    // FIXME_SYCL
#ifdef KOKKOS_ENABLE_SYCL
    const int s = 13;
#else
    const int s = 45;
#endif
    const int step_sizes[2] = {1, 10000};
    Kokkos::View<int**, MemorySpace> view("v", s * step_sizes[0],
                                          (s + 1) * step_sizes[1]);
    Kokkos::deep_copy(exec, view, 1);
    for (int step = 2; step < view.extent_int(0); ++step) {
      auto subview =
          Kokkos::subview(view, std::make_pair(0, (step + 1) * step_sizes[0]),
                          std::make_pair(0, (step + 2) * step_sizes[1]));
      Kokkos::View<int**, MemorySpace> subview_copy(
          "subview_copy", subview.extent(0), subview.extent(1));
      Kokkos::deep_copy(TEST_EXECSPACE{}, subview_copy, subview);
      exec.fence();

      SumView<decltype(subview)> sum_subview(subview);
      int total_subview = sum_subview.run();
      SumView<decltype(subview_copy)> sum_subview_copy(subview_copy);
      int total_subview_copy = sum_subview_copy.run();

      ASSERT_EQ(total_subview, total_subview_copy);
    }
  }
};

// Check that deep_copy with a large range for a dimension different from the
// first one works successfully. There was a problem with this in the Cuda
// backend.
TEST(TEST_CATEGORY, mdrange_large_deep_copy) {
  TestMDRangeLargeDeepCopy<TEST_EXECSPACE>::run();
}

}  // namespace Test
