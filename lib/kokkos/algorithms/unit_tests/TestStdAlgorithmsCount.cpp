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
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace Count {

namespace KE = Kokkos::Experimental;

template <class ViewType>
void test_count(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;

  view_host_space_t expected("count_expected", view.extent(0));
  compare_views(expected, view);

  {
    const value_t count_value = 0;
    const auto std_result =
        std::count(KE::cbegin(expected), KE::cend(expected), count_value);
    ASSERT_EQ(view.extent(0), size_t(std_result));

    // pass const iterators
    ASSERT_EQ(std_result, KE::count(exespace(), KE::cbegin(view),
                                    KE::cend(view), count_value));
    // pass view
    ASSERT_EQ(std_result, KE::count(exespace(), view, count_value));
  }

  {
    const value_t count_value = 13;
    const auto std_result =
        std::count(KE::cbegin(expected), KE::cend(expected), count_value);

    // pass iterators
    ASSERT_EQ(std_result, KE::count("label", exespace(), KE::begin(view),
                                    KE::end(view), count_value));
    // pass view
    ASSERT_EQ(std_result, KE::count("label", exespace(), view, count_value));
  }
}

template <class ViewType>
void test_count_if(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;

  view_host_space_t expected("count_expected", view.extent(0));
  compare_views(expected, view);

  // no positive elements (all zeroes)
  const auto predicate = IsPositiveFunctor<value_type>();
  ASSERT_EQ(0,
            std::count_if(KE::begin(expected), KE::end(expected), predicate));

  // pass iterators
  ASSERT_EQ(
      0, KE::count_if(exespace(), KE::begin(view), KE::end(view), predicate));
  // pass view
  ASSERT_EQ(0, KE::count_if(exespace(), view, predicate));

  fill_views_inc(view, expected);

  const auto std_result =
      std::count_if(KE::begin(expected), KE::end(expected), predicate);
  // pass const iterators
  ASSERT_EQ(std_result, KE::count_if("label", exespace(), KE::cbegin(view),
                                     KE::cend(view), predicate));
  // pass view
  ASSERT_EQ(std_result, KE::count_if("label", exespace(), view, predicate));
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  for (const auto& scenario : default_scenarios) {
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "count");
      test_count(view);
    }
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "count");
      test_count_if(view);
    }
  }
}

TEST(std_algorithms_count_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoTag, int>();
  run_all_scenarios<StridedThreeTag, unsigned>();
}

}  // namespace Count
}  // namespace stdalgos
}  // namespace Test
