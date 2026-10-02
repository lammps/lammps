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
template <class DataType, class ExecSpace, class LayOut>
static void test_required_allocation_size() {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using DynRankType  = Kokkos::DynRankView<DataType, LayOut, ExecSpace>;
  const size_t bytes = sizeof(DataType);
  auto size_length_1 = DynRankType::required_allocation_size(10);
  ASSERT_EQ(size_length_1, bytes * 10);
  auto size_length_2 = DynRankType::required_allocation_size(10, 20);
  ASSERT_EQ(size_length_2, bytes * 10 * 20);
  auto size_length_3 = DynRankType::required_allocation_size(10, 20, 3);
  ASSERT_EQ(size_length_3, bytes * 10 * 20 * 3);
  auto size_length_7 =
      DynRankType::required_allocation_size(2, 3, 4, 5, 6, 7, 8);
  ASSERT_EQ(size_length_7, bytes * 2 * 3 * 4 * 5 * 6 * 7 * 8);
}

template <class DataType, class ExecSpace, class LayOut>
static void test_required_allocation_size_death() {
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_5
#ifndef KOKKOS_ENABLE_DEBUG
  GTEST_SKIP() << "only enforced when debug checks are enabled";
  KOKKOS_IMPL_UNREACHABLE();
#endif
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  using DynRankType = Kokkos::DynRankView<DataType, LayOut, ExecSpace>;
  ASSERT_DEATH(DynRankType::required_allocation_size(2, 3, 4, 5, 6, 7, 8, 9),
               "Cannot allocate 8 dimensions");
#endif
}

TEST(TEST_CATEGORY, dyn_rank_view_ctor_from_members) {
  test_dyn_rank_view_ctor_from_members();
}

TEST(TEST_CATEGORY, dyn_rank_view_required_allocation_size) {
  test_required_allocation_size<double, TEST_EXECSPACE, Kokkos::LayoutRight>();
  test_required_allocation_size<int, TEST_EXECSPACE, Kokkos::LayoutRight>();
}

TEST(TEST_CATEGORY_DEATH, dyn_rank_view_required_allocation_size_death) {
  test_required_allocation_size_death<double, TEST_EXECSPACE,
                                      Kokkos::LayoutRight>();
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
