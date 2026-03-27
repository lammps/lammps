// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

namespace {

TEST(TEST_CATEGORY_DEATH, view_subview_constructor_layout_compatibility) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  int N    = 10;
  using LR = Kokkos::LayoutRight;
  using LL = Kokkos::LayoutLeft;
  using LS = Kokkos::LayoutStride;

  Kokkos::View<int**, LL> a1("A1", N, N);
  {
    // Using subview dims (ALL, 1). For a LayoutLeft,
    // any subview layout should be appropriate.
    (void)Kokkos::View<int*, LL>(a1, Kokkos::ALL, 1);
    (void)Kokkos::View<int*, LS>(a1, Kokkos::ALL, 1);

    // FIXME: This doesn't compile for BasicView, but should
    //(void)Kokkos::View<int*, LR>(a, Kokkos::ALL, 1);
  }
  {
    // Using subview dims (1, ALL). For a LayoutLeft,
    // resulting subview must be strided.
    const std::string msg = "View assignment must have compatible layouts";
    ASSERT_DEATH(((void)Kokkos::View<int*, LL>(a1, 1, Kokkos::ALL)), msg);
    (void)Kokkos::View<int*, LS>(a1, 1, Kokkos::ALL);
    ASSERT_DEATH(((void)Kokkos::View<int*, LR>(a1, 1, Kokkos::ALL)), msg);
  }

  Kokkos::View<int**, LL> a2("A2", 1, N);
  {
    // Using subview dims (0, ALL). Any subview layout should be appropriate.
    (void)Kokkos::View<int*, LL>(a2, 0, Kokkos::ALL);
    (void)Kokkos::View<int*, LS>(a2, 0, Kokkos::ALL);
    (void)Kokkos::View<int*, LR>(a2, 0, Kokkos::ALL);
  }
}
}  // namespace
