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
#include <iterator>
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace Find {

namespace KE = Kokkos::Experimental;

template <class ViewType>
void test_find(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;

  view_host_space_t expected("count_expected", view.extent(0));
  compare_views(expected, view);
  constexpr value_t find_value = 13;

  // value not found, return last
  ASSERT_EQ(KE::end(expected),
            std::find(KE::begin(expected), KE::end(expected), find_value));

  // pass const iterators, returns const iterator
  ASSERT_EQ(KE::cend(view),
            KE::find(exespace(), KE::cbegin(view), KE::cend(view), find_value));
  // pass view, returns iterator
  ASSERT_EQ(KE::end(view), KE::find(exespace(), view, find_value));

  fill_views_inc(view, expected);

  auto std_result =
      std::find(KE::begin(expected), KE::end(expected), find_value);
  auto distance = std::distance(KE::begin(expected), std_result);

  // pass iterators, returns iterator
  ASSERT_EQ(KE::begin(view) + distance,
            KE::find(exespace(), KE::begin(view), KE::end(view), find_value));
  // pass view, returns iterator
  ASSERT_EQ(KE::begin(view) + distance, KE::find(exespace(), view, find_value));
}

template <class ViewType>
void test_find_if(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;

  view_host_space_t expected("count_expected", view.extent(0));
  compare_views(expected, view);

  const auto not_equals_zero = NotEqualsZeroFunctor<value_type>();

  // value not found, return last
  ASSERT_EQ(
      KE::end(expected),
      std::find_if(KE::begin(expected), KE::end(expected), not_equals_zero));

  // pass iterators, returns iterator
  ASSERT_EQ(KE::end(view), KE::find_if(exespace(), KE::begin(view),
                                       KE::end(view), not_equals_zero));
  // pass view, returns iterator
  ASSERT_EQ(KE::end(view), KE::find_if(exespace(), view, not_equals_zero));

  fill_views_inc(view, expected);

  constexpr value_t find_value = 13;
  const auto equals_val        = EqualsValFunctor<value_type>(find_value);
  auto std_result =
      std::find_if(KE::begin(expected), KE::end(expected), equals_val);
  auto distance = std::distance(KE::begin(expected), std_result);

  // pass const iterators, returns const iterator
  ASSERT_EQ(
      KE::cbegin(view) + distance,
      KE::find_if(exespace(), KE::cbegin(view), KE::cend(view), equals_val));
  // pass view, returns iterator
  ASSERT_EQ(KE::begin(view) + distance,
            KE::find_if(exespace(), view, equals_val));
}

template <class ViewType>
void test_find_if_not(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;

  view_host_space_t expected("count_expected", view.extent(0));
  compare_views(expected, view);

  const auto not_equals_zero = NotEqualsZeroFunctor<value_type>();

  // first value matches
  ASSERT_EQ(KE::begin(expected),
            std::find_if_not(KE::begin(expected), KE::end(expected),
                             not_equals_zero));

  // pass iterators, returns iterator
  ASSERT_EQ(KE::begin(view), KE::find_if_not(exespace(), KE::begin(view),
                                             KE::end(view), not_equals_zero));
  // pass view, returns iterator
  ASSERT_EQ(KE::begin(view),
            KE::find_if_not(exespace(), view, not_equals_zero));

  fill_views_inc(view, expected);

  const auto equals_zero = EqualsValFunctor<value_type>(0);
  auto std_result =
      std::find_if_not(KE::begin(expected), KE::end(expected), equals_zero);
  auto distance = std::distance(KE::begin(expected), std_result);

  // pass const iterators, returns const iterator
  ASSERT_EQ(KE::cbegin(view) + distance,
            KE::find_if_not(exespace(), KE::cbegin(view), KE::cend(view),
                            equals_zero));
  // pass view, returns const iterator
  ASSERT_EQ(KE::begin(view) + distance,
            KE::find_if_not(exespace(), view, equals_zero));
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  for (const auto& scenario : default_scenarios) {
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "find");
      test_find(view);
    }
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "find_if");
      test_find_if(view);
    }
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "find_if_not");
      test_find_if_not(view);
    }
  }
}

TEST(std_algorithms_find_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoTag, int>();
  run_all_scenarios<StridedThreeTag, unsigned>();
}

}  // namespace Find
}  // namespace stdalgos
}  // namespace Test
