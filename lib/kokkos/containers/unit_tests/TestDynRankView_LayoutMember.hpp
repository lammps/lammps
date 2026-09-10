// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_DynRankView.hpp>

namespace {

template <class Layout>
void test_dyn_rank_view_layout_member() {
  constexpr bool is_ll = std::is_same_v<Layout, Kokkos::LayoutLeft>;
  {
    Kokkos::DynRankView<int, Layout> a(
        Kokkos::View<int***, Layout>("A", 11, 7, 5));
    auto l = a.layout();
    ASSERT_EQ(l.dimension[0], 11lu);
    ASSERT_EQ(l.dimension[1], 7lu);
    ASSERT_EQ(l.dimension[2], 5lu);
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
    ASSERT_TRUE((l.stride == is_ll ? 11lu : KOKKOS_INVALID_INDEX));
#else
    ASSERT_TRUE((l.stride == is_ll ? 11lu : 5lu));
#endif
  }
  {
    Kokkos::DynRankView<int, Layout> a(Kokkos::View<int**, Layout>("A", 7, 5));
    auto l = a.layout();
    ASSERT_EQ(l.dimension[0], 7lu);
    ASSERT_EQ(l.dimension[1], 5lu);
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
    ASSERT_TRUE((l.stride == is_ll ? 7lu : KOKKOS_INVALID_INDEX));
#else
    ASSERT_TRUE((l.stride == is_ll ? 7lu : 5lu));
#endif
  }
}

}  // namespace

TEST(TEST_CATEGORY, dyn_rank_view_layout_member) {
  test_dyn_rank_view_layout_member<Kokkos::LayoutRight>();
  test_dyn_rank_view_layout_member<Kokkos::LayoutLeft>();
}

TEST(TEST_CATEGORY, dyn_rank_view_layout_bug) {
  Kokkos::DynRankView<double, Kokkos::LayoutRight> a("A", 10, 11);
  // The following call failed with an error with legacy View off
  Kokkos::DynRankView<double, Kokkos::LayoutRight> b("B", a.layout());
#ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  ASSERT_EQ(a.mapping(), b.mapping());
#endif
}
