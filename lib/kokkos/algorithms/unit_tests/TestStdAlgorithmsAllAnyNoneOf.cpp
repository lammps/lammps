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

#include <TestStdAlgorithmsCommon.hpp>
#include <std_algorithms/Kokkos_BeginEnd.hpp>
#include <std_algorithms/Kokkos_AllOf.hpp>
#include <std_algorithms/Kokkos_AnyOf.hpp>
#include <std_algorithms/Kokkos_NoneOf.hpp>
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace AllAnyNoneOf {

namespace KE = Kokkos::Experimental;

template <class ViewType>
void test_all_of(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;
  const auto equals_zero  = EqualsValFunctor<value_t>(0);

  view_host_space_t expected("all_of_expected", view.extent(0));
  compare_views(expected, view);

  // reference result
  EXPECT_TRUE(std::all_of(KE::begin(expected), KE::end(expected), equals_zero));

  // pass iterators
  EXPECT_TRUE(
      KE::all_of(exespace(), KE::begin(view), KE::end(view), equals_zero));
  // pass view
  EXPECT_TRUE(KE::all_of(exespace(), view, equals_zero));

  fill_views_inc(view, expected);

  if (view.extent(0) > 1) {
    // reference result
    EXPECT_FALSE(
        std::all_of(KE::begin(expected), KE::end(expected), equals_zero));

    // pass const iterators
    EXPECT_FALSE(
        KE::all_of(exespace(), KE::cbegin(view), KE::cend(view), equals_zero));
    // pass view
    EXPECT_FALSE(KE::all_of("label", exespace(), view, equals_zero));
  }
}

template <class ViewType>
void test_any_of(const ViewType view) {
  using value_t              = typename ViewType::value_type;
  using view_host_space_t    = Kokkos::View<value_t*, Kokkos::HostSpace>;
  const auto not_equals_zero = NotEqualsZeroFunctor<value_t>();

  view_host_space_t expected("any_of_expected", view.extent(0));
  compare_views(expected, view);

  // reference result
  EXPECT_FALSE(
      std::any_of(KE::begin(expected), KE::end(expected), not_equals_zero));

  // pass iterators
  EXPECT_FALSE(
      KE::any_of(exespace(), KE::begin(view), KE::end(view), not_equals_zero));
  // pass view
  EXPECT_FALSE(KE::any_of(exespace(), view, not_equals_zero));

  fill_views_inc(view, expected);

  if (view.extent(0) > 1) {
    // reference result
    EXPECT_TRUE(
        std::any_of(KE::begin(expected), KE::end(expected), not_equals_zero));

    // pass const iterators
    EXPECT_TRUE(KE::any_of(exespace(), KE::cbegin(view), KE::cend(view),
                           not_equals_zero));
    // pass view
    EXPECT_TRUE(KE::any_of("label", exespace(), view, not_equals_zero));
  }
}

template <class ViewType>
void test_none_of(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;
  const auto is_positive  = IsPositiveFunctor<value_t>();

  view_host_space_t expected("none_of_expected", view.extent(0));
  compare_views(expected, view);

  // reference result
  EXPECT_TRUE(
      std::none_of(KE::begin(expected), KE::end(expected), is_positive));

  // pass iterators
  EXPECT_TRUE(
      KE::none_of(exespace(), KE::begin(view), KE::end(view), is_positive));
  // pass view
  EXPECT_TRUE(KE::none_of(exespace(), view, is_positive));

  fill_views_inc(view, expected);

  if (view.extent(0) > 1) {
    // reference result
    EXPECT_FALSE(
        std::none_of(KE::begin(expected), KE::end(expected), is_positive));

    // pass const iterators
    EXPECT_FALSE(
        KE::none_of(exespace(), KE::cbegin(view), KE::cend(view), is_positive));
    // pass view
    EXPECT_FALSE(KE::none_of("label", exespace(), view, is_positive));
  }
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  for (const auto& scenario : default_scenarios) {
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "all_of");
      test_all_of(view);
    }
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "any_of");
      test_any_of(view);
    }
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "none_of");
      test_none_of(view);
    }
  }
}

TEST(std_algorithms_all_any_none_of_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoTag, int>();
  run_all_scenarios<StridedThreeTag, unsigned>();
}

}  // namespace AllAnyNoneOf
}  // namespace stdalgos
}  // namespace Test
