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

#include <TestDynViewAPI.hpp>

namespace Test {
TEST(TEST_CATEGORY, dyn_rank_view_api_operator_rank12345) {
  TestDynViewAPI<double, TEST_EXECSPACE>::run_operator_test_rank12345();
}

template <typename SharedMemorySpace>
void test_dyn_rank_view_resize() {
  int n = 1000;
  Kokkos::DynRankView<double, SharedMemorySpace> device_view("device view", n);
  // Make sure we don't deallocate memory in Kokkos::resize
  auto device_view_copy = device_view;

  Kokkos::resize(device_view, 2 * n);

  // Loop in reverse to increase likelihood of missing fence detection assuming
  // that resize copies values in order.
  for (int i = 2 * n - 1; i >= 0; --i) device_view(i) = i + 1;

  Kokkos::fence();

  // Check that Kokkos::resize completed before setting the values on the host
  // manually (possibly because of missing fences).
  for (int i = 0; i < 2 * n; ++i) ASSERT_EQ(device_view(i), i + 1);
}

template <typename SharedMemorySpace>
void test_dyn_rank_view_realloc() {
  int n = 1000;
  Kokkos::DynRankView<double, SharedMemorySpace> device_view("device view", n);
  // Make sure we don't deallocate memory in Kokkos::realloc
  auto device_view_copy = device_view;

  Kokkos::realloc(device_view, 2 * n);

  // Loop in reverse to increase likelihood of missing fence detection assuming
  // that realloc sets values in order.
  for (int i = 2 * n - 1; i >= 0; --i) device_view(i) = i + 1;

  Kokkos::fence();

  // Check that Kokkos::realloc completed before setting the values on the host
  // manually (possibly because of missing fences).
  for (int i = 0; i < 2 * n; ++i) ASSERT_EQ(device_view(i), i + 1);
}

#ifdef KOKKOS_HAS_SHARED_SPACE
TEST(TEST_CATEGORY, dyn_rank_view_check_fence_resize_realloc) {
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::DefaultExecutionSpace>) {
    test_dyn_rank_view_resize<Kokkos::SharedSpace>();
    test_dyn_rank_view_realloc<Kokkos::SharedSpace>();
  } else {
    GTEST_SKIP() << "skipping since not default execution space";
  }
}
#endif

}  // namespace Test
