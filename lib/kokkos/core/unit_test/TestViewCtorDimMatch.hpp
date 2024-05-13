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

namespace Test {

template <int rank, int dynrank, class RankType, std::size_t... Is>
void test_matching_arguments_rank_helper(std::index_sequence<Is...>) {
  constexpr int nargs = sizeof...(Is);
  using view_type     = Kokkos::View<RankType>;
  if (nargs == rank || nargs == dynrank) {
    EXPECT_NO_THROW({ view_type v("v", ((Is * 0) + 1)...); });
    EXPECT_NO_THROW({ view_type v(nullptr, ((Is * 0) + 1)...); });
  } else {
    ASSERT_DEATH(
        { view_type v("v", ((Is * 0) + 1)...); },
        "Constructor for Kokkos::View 'v' has mismatched number of arguments. "
        "The number of arguments = " +
            std::to_string(nargs) +
            " neither matches the dynamic rank = " + std::to_string(dynrank) +
            " nor the total rank = " + std::to_string(rank));
    ASSERT_DEATH(
        { view_type v(nullptr, ((Is * 0) + 1)...); },
        "Constructor for Kokkos::View 'UNMANAGED' has mismatched number of "
        "arguments. "
        "The number of arguments = " +
            std::to_string(nargs) +
            " neither matches the dynamic rank = " + std::to_string(dynrank) +
            " nor the total rank = " + std::to_string(rank));
  }
}

template <int rank, int dynrank, template <int> class RankType>
void test_matching_arguments_rank() {
  test_matching_arguments_rank_helper<rank, dynrank,
                                      typename RankType<rank>::type>(
      std::make_index_sequence<0>());
  test_matching_arguments_rank_helper<rank, dynrank,
                                      typename RankType<rank>::type>(
      std::make_index_sequence<1>());
  test_matching_arguments_rank_helper<rank, dynrank,
                                      typename RankType<rank>::type>(
      std::make_index_sequence<2>());
  test_matching_arguments_rank_helper<rank, dynrank,
                                      typename RankType<rank>::type>(
      std::make_index_sequence<3>());
  test_matching_arguments_rank_helper<rank, dynrank,
                                      typename RankType<rank>::type>(
      std::make_index_sequence<4>());
  test_matching_arguments_rank_helper<rank, dynrank,
                                      typename RankType<rank>::type>(
      std::make_index_sequence<5>());
  test_matching_arguments_rank_helper<rank, dynrank,
                                      typename RankType<rank>::type>(
      std::make_index_sequence<6>());
  test_matching_arguments_rank_helper<rank, dynrank,
                                      typename RankType<rank>::type>(
      std::make_index_sequence<7>());
  test_matching_arguments_rank_helper<rank, dynrank,
                                      typename RankType<rank>::type>(
      std::make_index_sequence<8>());
}

template <int rank>
struct DynamicRank {
  using type = typename DynamicRank<rank - 1>::type*;
};

template <>
struct DynamicRank<0> {
  using type = int;
};

// Skip test execution when KOKKOS_ENABLE_OPENMPTARGET is enabled until
// Kokkos::abort() aborts properly on that backend
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY_DEATH, view_construction_with_wrong_params_dyn) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECKS
  test_matching_arguments_rank<0, 0, DynamicRank>();  // dim = 0, dynamic = 0
  test_matching_arguments_rank<1, 1, DynamicRank>();  // dim = 1, dynamic = 1
  test_matching_arguments_rank<2, 2, DynamicRank>();  // dim = 2, dynamic = 2
  test_matching_arguments_rank<3, 3, DynamicRank>();  // dim = 3, dynamic = 3
  test_matching_arguments_rank<4, 4, DynamicRank>();  // dim = 4, dynamic = 4
  test_matching_arguments_rank<5, 5, DynamicRank>();  // dim = 5, dynamic = 5
  test_matching_arguments_rank<6, 6, DynamicRank>();  // dim = 6, dynamic = 6
  test_matching_arguments_rank<7, 7, DynamicRank>();  // dim = 7, dynamic = 7
  test_matching_arguments_rank<8, 8, DynamicRank>();  // dim = 8, dynamic = 8
#endif
}

template <int rank>
struct StaticRank {
  using type = typename StaticRank<rank - 1>::type[1];
};

template <>
struct StaticRank<0> {
  using type = int;
};

TEST(TEST_CATEGORY_DEATH, view_construction_with_wrong_params_stat) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECKS
  test_matching_arguments_rank<0, 0, StaticRank>();  // dim = 0, dynamic = 0
  test_matching_arguments_rank<1, 0, StaticRank>();  // dim = 1, dynamic = 0
  test_matching_arguments_rank<2, 0, StaticRank>();  // dim = 2, dynamic = 0
  test_matching_arguments_rank<3, 0, StaticRank>();  // dim = 3, dynamic = 0
  test_matching_arguments_rank<4, 0, StaticRank>();  // dim = 4, dynamic = 0
  test_matching_arguments_rank<5, 0, StaticRank>();  // dim = 5, dynamic = 0
  test_matching_arguments_rank<6, 0, StaticRank>();  // dim = 6, dynamic = 0
  test_matching_arguments_rank<7, 0, StaticRank>();  // dim = 7, dynamic = 0
  test_matching_arguments_rank<8, 0, StaticRank>();  // dim = 8, dynamic = 0
#endif
}

template <int rank>
struct MixedRank {
  using type = typename DynamicRank<rank - 1>::type[1];
};

template <>
struct MixedRank<0> {
  using type = int;
};

TEST(TEST_CATEGORY_DEATH, view_construction_with_wrong_params_mix) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECKS
  test_matching_arguments_rank<0, 0, MixedRank>();  // dim = 0, dynamic = 0
  test_matching_arguments_rank<1, 0, MixedRank>();  // dim = 1, dynamic = 0
  test_matching_arguments_rank<2, 1, MixedRank>();  // dim = 2, dynamic = 1
  test_matching_arguments_rank<3, 2, MixedRank>();  // dim = 3, dynamic = 2
  test_matching_arguments_rank<4, 3, MixedRank>();  // dim = 4, dynamic = 3
  test_matching_arguments_rank<5, 4, MixedRank>();  // dim = 5, dynamic = 4
  test_matching_arguments_rank<6, 5, MixedRank>();  // dim = 6, dynamic = 5
  test_matching_arguments_rank<7, 6, MixedRank>();  // dim = 7, dynamic = 6
  test_matching_arguments_rank<8, 7, MixedRank>();  // dim = 8, dynamic = 7
#endif
}

#define CHECK_DEATH(EXPR)                                                     \
  ASSERT_DEATH(EXPR,                                                          \
               "The specified run-time extent for Kokkos::View 'v' does not " \
               "match the compile-time extent in dimension 0. The given "     \
               "extent is 2 but should be 1.")

#define CHECK_DEATH_UNMANAGED(EXPR)                                          \
  ASSERT_DEATH(                                                              \
      EXPR,                                                                  \
      "The specified run-time extent for Kokkos::View 'UNMANAGED' does not " \
      "match the compile-time extent in dimension 0. The given "             \
      "extent is 2 but should be 1.")

TEST(TEST_CATEGORY_DEATH, view_construction_with_wrong_static_extents) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECKS
  // clang-format off
  CHECK_DEATH({ Kokkos::View<int[1]>                      v("v", 2); });
  CHECK_DEATH({ Kokkos::View<int[1][1]>                   v("v", 2, 1); });
  CHECK_DEATH({ Kokkos::View<int[1][1][1]>                v("v", 2, 1, 1); });
  CHECK_DEATH({ Kokkos::View<int[1][1][1][1]>             v("v", 2, 1, 1, 1); });
  CHECK_DEATH({ Kokkos::View<int[1][1][1][1][1]>          v("v", 2, 1, 1, 1, 1); });
  CHECK_DEATH({ Kokkos::View<int[1][1][1][1][1][1]>       v("v", 2, 1, 1, 1, 1, 1); });
  CHECK_DEATH({ Kokkos::View<int[1][1][1][1][1][1][1]>    v("v", 2, 1, 1, 1, 1, 1, 1); });
  CHECK_DEATH({ Kokkos::View<int[1][1][1][1][1][1][1][1]> v("v", 2, 1, 1, 1, 1, 1, 1, 1); });

  CHECK_DEATH_UNMANAGED({ Kokkos::View<int[1]>                      v(nullptr, 2); });
  CHECK_DEATH_UNMANAGED({ Kokkos::View<int[1][1]>                   v(nullptr, 2, 1); });
  CHECK_DEATH_UNMANAGED({ Kokkos::View<int[1][1][1]>                v(nullptr, 2, 1, 1); });
  CHECK_DEATH_UNMANAGED({ Kokkos::View<int[1][1][1][1]>             v(nullptr, 2, 1, 1, 1); });
  CHECK_DEATH_UNMANAGED({ Kokkos::View<int[1][1][1][1][1]>          v(nullptr, 2, 1, 1, 1, 1); });
  CHECK_DEATH_UNMANAGED({ Kokkos::View<int[1][1][1][1][1][1]>       v(nullptr, 2, 1, 1, 1, 1, 1); });
  CHECK_DEATH_UNMANAGED({ Kokkos::View<int[1][1][1][1][1][1][1]>    v(nullptr, 2, 1, 1, 1, 1, 1, 1); });
  CHECK_DEATH_UNMANAGED({ Kokkos::View<int[1][1][1][1][1][1][1][1]> v(nullptr, 2, 1, 1, 1, 1, 1, 1, 1); });
  // clang-format on
#endif
}

#undef CHECK_DEATH
#endif  // KOKKOS_ENABLE_OPENMPTARGET
}  // namespace Test
