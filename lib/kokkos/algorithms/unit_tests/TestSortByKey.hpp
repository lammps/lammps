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

#ifndef KOKKOS_ALGORITHMS_UNITTESTS_TEST_SORT_BY_KEY_HPP
#define KOKKOS_ALGORITHMS_UNITTESTS_TEST_SORT_BY_KEY_HPP

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Sort.hpp>

#include <utility>  // pair

namespace Test {
namespace SortImpl {

struct Less {
  template <class ValueType>
  KOKKOS_INLINE_FUNCTION bool operator()(const ValueType &lhs,
                                         const ValueType &rhs) const {
    return lhs < rhs;
  }
};

struct Greater {
  template <class ValueType>
  KOKKOS_INLINE_FUNCTION bool operator()(const ValueType &lhs,
                                         const ValueType &rhs) const {
    return lhs > rhs;
  }
};

template <class ExecutionSpace, class Keys, class Permute,
          class Comparator = Less>
struct is_sorted_by_key_struct {
  Keys keys;
  Keys keys_orig;
  Permute permute;
  Comparator comparator;

  is_sorted_by_key_struct(Keys keys_, Keys keys_orig_, Permute permute_,
                          Comparator comparator_ = Comparator{})
      : keys(keys_),
        keys_orig(keys_orig_),
        permute(permute_),
        comparator(comparator_) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(int i, unsigned int &count) const {
    if (i < keys.extent_int(0) - 1 && comparator(keys(i + 1), keys(i))) ++count;
    if (keys(i) != keys_orig(permute(i))) ++count;
  }
};

template <typename ExecutionSpace, typename ViewType>
void iota(ExecutionSpace const &space, ViewType const &v,
          typename ViewType::value_type value = 0) {
  using ValueType = typename ViewType::value_type;
  Kokkos::parallel_for(
      "Kokkos::Algorithms::iota",
      Kokkos::RangePolicy<ExecutionSpace>(space, 0, v.extent(0)),
      KOKKOS_LAMBDA(int i) { v(i) = value + (ValueType)i; });
}

}  // namespace SortImpl

TEST(TEST_CATEGORY, SortByKeyEmptyView) {
  using ExecutionSpace = TEST_EXECSPACE;

  // does not matter if we use int or something else
  Kokkos::View<int *, ExecutionSpace> keys("keys", 0);
  Kokkos::View<float *, ExecutionSpace> values("values", 0);

  ASSERT_NO_THROW(
      Kokkos::Experimental::sort_by_key(ExecutionSpace(), keys, values));
}

// Test #7036
TEST(TEST_CATEGORY, SortByKeyEmptyViewHost) {
  using ExecutionSpace = Kokkos::DefaultHostExecutionSpace;

  // does not matter if we use int or something else
  Kokkos::View<int *, ExecutionSpace> keys("keys", 0);
  Kokkos::View<float *, ExecutionSpace> values("values", 0);

  ASSERT_NO_THROW(
      Kokkos::Experimental::sort_by_key(ExecutionSpace(), keys, values));
}

TEST(TEST_CATEGORY, SortByKey) {
  using ExecutionSpace = TEST_EXECSPACE;
  using MemorySpace    = typename ExecutionSpace::memory_space;

  ExecutionSpace space{};

  for (auto keys_vector : {std::vector<int>{36, 19, 25, 17, 3, 7, 1, 2, 9},
                           std::vector<int>{36, 19, 25, 17, 3, 9, 1, 2, 7},
                           std::vector<int>{100, 19, 36, 17, 3, 25, 1, 2, 7},
                           std::vector<int>{15, 5, 11, 3, 4, 8}}) {
    auto const n = keys_vector.size();

    auto keys = Kokkos::create_mirror_view_and_copy(
        MemorySpace{},
        Kokkos::View<int *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
            keys_vector.data(), n));

    auto keys_orig = Kokkos::create_mirror(space, keys);
    Kokkos::deep_copy(space, keys_orig, keys);

    Kokkos::View<int *, ExecutionSpace> permute("permute", n);
    SortImpl::iota(space, permute);

    Kokkos::Experimental::sort_by_key(space, keys, permute);

    unsigned int sort_fails = 0;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<ExecutionSpace>(space, 0, n),
        SortImpl::is_sorted_by_key_struct<ExecutionSpace, decltype(keys),
                                          decltype(permute)>(keys, keys_orig,
                                                             permute),
        sort_fails);

    ASSERT_EQ(sort_fails, 0u);
  }
}

TEST(TEST_CATEGORY, SortByKeyWithComparator) {
  using ExecutionSpace = TEST_EXECSPACE;
  using MemorySpace    = typename ExecutionSpace::memory_space;

  ExecutionSpace space{};

  SortImpl::Greater comparator;

  for (auto keys_vector : {std::vector<int>{36, 19, 25, 17, 3, 7, 1, 2, 9},
                           std::vector<int>{36, 19, 25, 17, 3, 9, 1, 2, 7},
                           std::vector<int>{100, 19, 36, 17, 3, 25, 1, 2, 7},
                           std::vector<int>{15, 5, 11, 3, 4, 8}}) {
    auto const n = keys_vector.size();

    auto keys = Kokkos::create_mirror_view_and_copy(
        MemorySpace{},
        Kokkos::View<int *, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>(
            keys_vector.data(), n));

    auto keys_orig = Kokkos::create_mirror(space, keys);
    Kokkos::deep_copy(space, keys_orig, keys);

    Kokkos::View<int *, ExecutionSpace> permute("permute", n);
    SortImpl::iota(space, permute);

    Kokkos::Experimental::sort_by_key(space, keys, permute, comparator);

    unsigned int sort_fails = 0;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<ExecutionSpace>(space, 0, n),
        SortImpl::is_sorted_by_key_struct<ExecutionSpace, decltype(keys),
                                          decltype(permute), SortImpl::Greater>(
            keys, keys_orig, permute, comparator),
        sort_fails);

    ASSERT_EQ(sort_fails, 0u);
  }
}

TEST(TEST_CATEGORY, SortByKeyStaticExtents) {
  using ExecutionSpace = TEST_EXECSPACE;

  ExecutionSpace space{};

  Kokkos::View<int[10], ExecutionSpace> keys("keys");

  Kokkos::View<int[10], ExecutionSpace> values_static("values_static");
  ASSERT_NO_THROW(
      Kokkos::Experimental::sort_by_key(space, keys, values_static));

  Kokkos::View<int *, ExecutionSpace> values_dynamic("values_dynamic", 10);
  ASSERT_NO_THROW(
      Kokkos::Experimental::sort_by_key(space, keys, values_dynamic));
}

template <typename ExecutionSpace, typename Keys, typename Values>
void buildViewsForStrided(ExecutionSpace const &space, int n, Keys &keys,
                          Values &values) {
  Kokkos::parallel_for(
      "create_data",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>(space, {0, 0, 0},
                                                             {n, n, n}),
      KOKKOS_LAMBDA(int i, int j, int k) {
        keys(i, j, k)   = n - i;
        values(i, j, k) = j;
      });
}

TEST(TEST_CATEGORY, SortByKeyWithStrides) {
  using ExecutionSpace = TEST_EXECSPACE;

  ExecutionSpace space{};

  auto const n = 10;

  Kokkos::View<int ***, ExecutionSpace> keys("keys", n, n, n);
  Kokkos::View<int ***, ExecutionSpace> values("values", n, n, n);
  buildViewsForStrided(space, n, keys, values);

  auto keys_sub   = Kokkos::subview(keys, Kokkos::ALL(), 1, 2);
  auto values_sub = Kokkos::subview(values, 4, Kokkos::ALL(), 6);

  auto keys_orig = Kokkos::create_mirror(space, keys_sub);
  Kokkos::deep_copy(space, keys_orig, keys_sub);

  Kokkos::Experimental::sort_by_key(space, keys_sub, values_sub);

  unsigned int sort_fails = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecutionSpace>(space, 0, n),
      SortImpl::is_sorted_by_key_struct<ExecutionSpace, decltype(keys_sub),
                                        decltype(values_sub)>(
          keys_sub, keys_orig, values_sub),
      sort_fails);

  ASSERT_EQ(sort_fails, 0u);
}

TEST(TEST_CATEGORY, SortByKeyKeysLargerThanValues) {
  using ExecutionSpace = TEST_EXECSPACE;

  // does not matter if we use int or something else
  Kokkos::View<int *, ExecutionSpace> keys("keys", 3);
  Kokkos::View<float *, ExecutionSpace> values("values", 1);

  ASSERT_DEATH(
      Kokkos::Experimental::sort_by_key(ExecutionSpace(), keys, values),
      "values and keys extents must be the same");
  ASSERT_DEATH(Kokkos::Experimental::sort_by_key(ExecutionSpace(), keys, values,
                                                 SortImpl::Greater{}),
               "values and keys extents must be the same");
}

}  // namespace Test
#endif
