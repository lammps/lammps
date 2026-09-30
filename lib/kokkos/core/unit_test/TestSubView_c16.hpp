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
// Test that the subview constructor works for constructing subviews from
// Views where one of them has the Kokkos::Unmanaged memory trait
TEST(TEST_CATEGORY, view_create_unmanaged_subview_from_managed) {
  int N    = 10;
  using LR = Kokkos::LayoutRight;
  using LL = Kokkos::LayoutLeft;
  using LS = Kokkos::LayoutStride;

  Kokkos::View<int**, LL> a1("A1", N, N);
  Kokkos::View<int**, LL, Kokkos::MemoryTraits<Kokkos::Unmanaged>> a1_unmanaged(
      a1.data(), N, N);
  {
    // Using subview dims (ALL, 1). For a LayoutLeft,
    // any subview layout should be appropriate.
    (void)Kokkos::View<int*, LL, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
        a1, Kokkos::ALL, 1);
    (void)Kokkos::View<int*, LS, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
        a1, Kokkos::ALL, 1);

    (void)Kokkos::View<int*, LL>(a1_unmanaged, Kokkos::ALL, 1);
    (void)Kokkos::View<int*, LS>(a1_unmanaged, Kokkos::ALL, 1);
    // FIXME: This doesn't compile for BasicView, but should
    // (void)Kokkos::View<int*, LR>(a1_unmanaged, Kokkos::ALL, 1);
  }

  Kokkos::View<int**, LL> a2("A2", 1, N);
  Kokkos::View<int**, LL, Kokkos::MemoryTraits<Kokkos::Unmanaged>> a2_unmanaged(
      a2.data(), 1, N);
  {
    // Using subview dims (0, ALL). Any subview layout should be appropriate.
    (void)Kokkos::View<int*, LL, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
        a2, 0, Kokkos::ALL);
    (void)Kokkos::View<int*, LS, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
        a2, 0, Kokkos::ALL);
    (void)Kokkos::View<int*, LR, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
        a2, 0, Kokkos::ALL);

    (void)Kokkos::View<int*, LL>(a2_unmanaged, 0, Kokkos::ALL);
    (void)Kokkos::View<int*, LS>(a2_unmanaged, 0, Kokkos::ALL);
    (void)Kokkos::View<int*, LR>(a2_unmanaged, 0, Kokkos::ALL);
  }
}
}  // namespace
