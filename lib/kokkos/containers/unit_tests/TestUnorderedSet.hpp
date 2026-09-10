// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_UNORDERED_SET_HPP
#define KOKKOS_TEST_UNORDERED_SET_HPP

#include <random>

#include <gtest/gtest.h>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.random;
import kokkos.sort;
import kokkos.std_algorithms;
import kokkos.unordered_map;
#else
#include <Kokkos_Core.hpp>
#include <Kokkos_Random.hpp>
#include <Kokkos_Sort.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <Kokkos_UnorderedMap.hpp>
#endif

namespace {

template <typename MapType>
struct Inserter {
  MapType map;

  template <std::integral T>
  KOKKOS_FUNCTION void operator()(const T idx) const {
    Kokkos::UnorderedMapInsertResult res = map.insert(idx);
    if (!res.success()) Kokkos::abort("Inserter: insert failed.");
  }
};

template <typename MapType, typename KeyViewType>
struct Eraser {
  MapType map;
  KeyViewType keys;

  template <std::integral T>
  KOKKOS_FUNCTION void operator()(const T idx) const {
    if (!map.erase(keys(idx))) Kokkos::abort("Eraser: erase failed.");
  }
};

// Ensure that inserting into, erasing from, and rehashing an unordered set
// works.
TEST(TEST_CATEGORY, UnorderedSet_insert_erase_and_rehash) {
  // First, size_all entries will be inserted into the set.
  // Then, size_erased will be erased from it before rehash.
  // Once rehashing is done, the capacity of the set is expected to
  // be smaller than the one of the initial set.
  constexpr size_t size_all    = 2579;
  constexpr size_t size_erased = 568;

  using uset_type = Kokkos::UnorderedMap<size_t, void, TEST_EXECSPACE>;
  static_assert(uset_type::is_set);

  using key_view_type =
      Kokkos::View<typename uset_type::key_type*, TEST_EXECSPACE>;

  const TEST_EXECSPACE exec{};

  // Initialize the set.
  uset_type uset(Kokkos::view_alloc(exec, "test uset"), size_all);
  const auto initial_capacity = uset.capacity();

  ASSERT_GE(uset.capacity(), size_all);
  ASSERT_EQ(uset.size(), 0u);

  // Insert size_all keys.
  Kokkos::parallel_for(
      "uset test - insert " + std::to_string(size_all) + " values",
      Kokkos::RangePolicy(exec, 0, size_all), Inserter<uset_type>{.map = uset});
  exec.fence();

  ASSERT_GE(uset.capacity(), size_all);
  ASSERT_EQ(uset.size(), size_all);

  // Always try to erase truly random indices by initializing the random pool
  // randomly.
  std::mt19937 rng(testing::UnitTest::GetInstance()->random_seed());
  const uint64_t seed = std::uniform_int_distribution<uint64_t>{}(rng);

  // Generate random indices between 0 and size_all. Those will be erased.
  Kokkos::Random_XorShift64_Pool<TEST_EXECSPACE> generator(seed);
  key_view_type keys_erased(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, exec,
                         "keys that will be erased"),
      size_erased);
  Kokkos::fill_random(exec, keys_erased, generator, size_all);
  Kokkos::sort(exec, keys_erased);

  // Filter out any duplicated index.
  namespace KE = Kokkos::Experimental;
  const auto keys_erased_end =
      KE::unique(exec, KE::begin(keys_erased), KE::end(keys_erased));
  const size_t keys_erased_size =
      KE::distance(KE::begin(keys_erased), keys_erased_end);
  const size_t expected_uset_size = size_all - keys_erased_size;

  // Start erasing at random keys.
  ASSERT_TRUE(uset.begin_erase());
  auto keys_erased_subview = Kokkos::subview(
      keys_erased, Kokkos::make_pair<size_t, size_t>(0, keys_erased_size));
  Kokkos::parallel_for(
      "uset test - erase " + std::to_string(keys_erased_size) + " values",
      Kokkos::RangePolicy(exec, 0, keys_erased_size),
      Eraser<uset_type, decltype(keys_erased_subview)>{
          .map = uset, .keys = std::move(keys_erased_subview)});
  exec.fence();
  ASSERT_TRUE(uset.end_erase());

  ASSERT_GE(uset.capacity(), size_all);
  ASSERT_EQ(uset.size(), expected_uset_size);

  // Rehashing should not increase the set capacity.
  uset.rehash(expected_uset_size);

  ASSERT_LE(uset.capacity(), initial_capacity);
  ASSERT_EQ(uset.size(), expected_uset_size);
}

}  // namespace

#endif  // KOKKOS_TEST_UNORDERED_SET_HPP
