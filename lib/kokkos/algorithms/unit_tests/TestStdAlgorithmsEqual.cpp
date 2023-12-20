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
namespace Equal {

namespace KE = Kokkos::Experimental;

template <class ViewType>
void test_equal(const ViewType view) {
  auto copy = create_deep_copyable_compatible_clone(view);

  // pass iterators
  EXPECT_TRUE(
      KE::equal(exespace(), KE::begin(view), KE::end(view), KE::begin(copy)));
  // pass views
  EXPECT_TRUE(KE::equal(exespace(), view, copy));

  // modify copy - make the last element different
  const auto extent = view.extent(0);
  if (extent > 0) {
    KE::fill(exespace(), KE::end(copy) - 1, KE::end(copy), 1);

    // pass const iterators
    EXPECT_FALSE(KE::equal(exespace(), KE::cbegin(view), KE::cend(view),
                           KE::cbegin(copy)));
    // pass views
    EXPECT_FALSE(KE::equal("label", exespace(), view, copy));
  }
}

template <class ViewType>
void test_equal_custom_comparator(const ViewType view) {
  using value_t = typename ViewType::value_type;
  const auto p  = CustomEqualityComparator<value_t>();
  auto copy     = create_deep_copyable_compatible_clone(view);

  // pass iterators
  EXPECT_TRUE(KE::equal(exespace(), KE::begin(view), KE::end(view),
                        KE::begin(copy), p));
  // pass views
  EXPECT_TRUE(KE::equal(exespace(), view, copy, p));

  // modify copy - make the last element different
  const auto extent = view.extent(0);
  if (extent > 0) {
    KE::fill(exespace(), KE::end(copy) - 1, KE::end(copy), 1);

    // pass const iterators
    EXPECT_FALSE(KE::equal("label", exespace(), KE::cbegin(view),
                           KE::cend(view), KE::cbegin(copy), p));
    // pass views
    EXPECT_FALSE(KE::equal(exespace(), view, copy, p));
  }
}

template <class ViewType>
void test_equal_4_iterators(const ViewType view) {
  using value_t = typename ViewType::value_type;
  const auto p  = CustomEqualityComparator<value_t>();
  auto copy     = create_deep_copyable_compatible_clone(view);

  // pass iterators
  EXPECT_TRUE(KE::equal(exespace(), KE::begin(view), KE::end(view),
                        KE::begin(copy), KE::end(copy)));
  // pass const and non-const iterators, custom comparator
  EXPECT_TRUE(KE::equal("label", exespace(), KE::cbegin(view), KE::cend(view),
                        KE::begin(copy), KE::end(copy), p));

  const auto extent = view.extent(0);
  if (extent > 0) {
    // use different length ranges, pass const iterators
    EXPECT_FALSE(KE::equal(exespace(), KE::cbegin(view), KE::cend(view),
                           KE::cbegin(copy), KE::cend(copy) - 1));

    // modify copy - make the last element different
    KE::fill(exespace(), KE::end(copy) - 1, KE::end(copy), 1);
    // pass const iterators
    EXPECT_FALSE(KE::equal(exespace(), KE::cbegin(view), KE::cend(view),
                           KE::cbegin(copy), KE::cend(copy)));
  }
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  for (const auto& scenario : default_scenarios) {
    auto view = create_view<ValueType>(Tag{}, scenario.second, "equal");
    test_equal(view);
    test_equal_custom_comparator(view);
    test_equal_4_iterators(view);
  }
}

TEST(std_algorithms_equal_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoTag, int>();
  run_all_scenarios<StridedThreeTag, unsigned>();
}

}  // namespace Equal
}  // namespace stdalgos
}  // namespace Test
