// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_DynRankView.hpp>

namespace {

void test_dyn_rank_view_ctor_from_members() {
  using drview_t = Kokkos::DynRankView<double, TEST_EXECSPACE>;
  using view_t   = typename drview_t::view_type;

  {
    view_t data("data", 2, 3, 4, 5, 6, 7, 8);
    drview_t drv(data, 7);

    ASSERT_EQ(drv.rank(), 7lu);
    for (size_t r = 0; r < drv.rank(); r++)
      ASSERT_EQ(static_cast<size_t>(drv.extent(r)), r + 2);
    ASSERT_EQ(data.use_count(), 2);
    ASSERT_EQ(data.data(), drv.data());
  }
  {
    view_t data("data", 2, 3, 4, 1, 1, 1, 1);
    drview_t drv(data, 3);

    ASSERT_EQ(drv.rank(), 3lu);
    for (size_t r = 0; r < drv.rank(); r++)
      ASSERT_EQ(static_cast<size_t>(drv.extent(r)), r + 2);
    ASSERT_EQ(data.use_count(), 2);
    ASSERT_EQ(data.data(), drv.data());
  }
  {
    view_t data("data", 1, 1, 1, 1, 1, 1, 1);
    drview_t drv(data, 0);

    ASSERT_EQ(drv.rank(), 0lu);
    ASSERT_EQ(data.use_count(), 2);
    ASSERT_EQ(data.data(), drv.data());
  }
}

TEST(TEST_CATEGORY, dyn_rank_view_ctor_from_members) {
  test_dyn_rank_view_ctor_from_members();
}

#ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
void test_dyn_rank_view_ctor_from_layout_stride() {
  using drview_t =
      Kokkos::DynRankView<double, Kokkos::LayoutStride, TEST_EXECSPACE>;
  using view_t = typename drview_t::view_type;
  {
    typename view_t::mapping_type map(
        Kokkos::dextents<int, 7>(2, 3, 4, 5, 6, 7, 8),
        std::array{1, 2, 6, 24, 120, 720, 5040});
    view_t data("data", map);
    drview_t drv(data, 7);

    ASSERT_EQ(drv.rank(), 7lu);
    for (size_t r = 0; r < drv.rank(); r++)
      ASSERT_EQ(static_cast<size_t>(drv.extent(r)), r + 2);
    ASSERT_EQ(data.use_count(), 2);
    ASSERT_EQ(data.data(), drv.data());
  }
  {
    typename view_t::mapping_type map(
        Kokkos::dextents<int, 7>(2, 3, 4, 1, 1, 1, 1),
        std::array{1, 2, 6, 1, 1, 1, 1});
    view_t data("data", map);
    drview_t drv(data, 3);

    ASSERT_EQ(drv.rank(), 3lu);
    for (size_t r = 0; r < drv.rank(); r++)
      ASSERT_EQ(static_cast<size_t>(drv.extent(r)), r + 2);
    ASSERT_EQ(data.use_count(), 2);
    ASSERT_EQ(data.data(), drv.data());
  }
  {
    typename view_t::mapping_type map(
        Kokkos::dextents<int, 7>(1, 1, 1, 1, 1, 1, 1),
        std::array{1, 1, 1, 1, 1, 1, 1});
    view_t data("data", map);
    drview_t drv(data, 0);

    ASSERT_EQ(drv.rank(), 0lu);
    ASSERT_EQ(data.use_count(), 2);
    ASSERT_EQ(data.data(), drv.data());
  }
}

TEST(TEST_CATEGORY, dyn_rank_view_ctor_from_layout_stride) {
  test_dyn_rank_view_ctor_from_layout_stride();
}
#endif

}  // namespace
