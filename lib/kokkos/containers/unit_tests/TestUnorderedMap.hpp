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

#ifndef KOKKOS_TEST_UNORDERED_MAP_HPP
#define KOKKOS_TEST_UNORDERED_MAP_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <Kokkos_UnorderedMap.hpp>

namespace Test {

namespace Impl {

template <typename MapType,
          typename InsertOp = typename MapType::default_op_type,
          bool Near = false, bool CheckValues = false>
struct TestInsert {
  using map_type        = MapType;
  using execution_space = typename map_type::execution_space;
  using value_type      = uint32_t;

  struct ExpectedValues {
    unsigned map_idx;
    typename map_type::value_type v;
  };
  using expected_values_type = Kokkos::View<ExpectedValues *, execution_space>;
  expected_values_type expected_values;

  map_type map;
  uint32_t inserts;
  uint32_t collisions;
  InsertOp insert_op;

  TestInsert(map_type arg_map, uint32_t arg_inserts, uint32_t arg_collisions)
      : map(arg_map), inserts(arg_inserts), collisions(arg_collisions) {
    auto len = map.capacity() > arg_inserts ? map.capacity() : arg_inserts;
    expected_values = expected_values_type("ExpectedValues", len);
  }

  void testit(bool rehash_on_fail = true) {
    execution_space().fence();

    uint32_t failed_count = 0;
    do {
      failed_count = 0;
      Kokkos::parallel_reduce(inserts, *this, failed_count);

      if (rehash_on_fail && failed_count > 0u) {
        const uint32_t new_capacity = map.capacity() +
                                      ((map.capacity() * 3ull) / 20u) +
                                      failed_count / collisions;
        map.rehash(new_capacity);
      }
    } while (rehash_on_fail && failed_count > 0u);

    // Trigger the m_size mutable bug.
    auto map_h = create_mirror(map);
    execution_space().fence();
    Kokkos::deep_copy(map_h, map);
    execution_space().fence();
    ASSERT_EQ(map_h.size(), map.size());

    if (!rehash_on_fail && CheckValues) {
      typename expected_values_type::HostMirror expected_values_h =
          create_mirror_view(expected_values);
      Kokkos::deep_copy(expected_values_h, expected_values);
      for (unsigned i = 0; i < map_h.size(); i++) {
        auto map_idx = expected_values_h(i).map_idx;
        if (map_idx != static_cast<unsigned>(~0)) {
          ASSERT_EQ(expected_values_h(map_idx).v, map_h.value_at(map_idx));
        }
      }
    }

    const unsigned int old_size = map_h.size();
    map_h.clear();
    ASSERT_EQ(map.size(), old_size);
    ASSERT_EQ(map_h.size(), 0u);
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type &failed_count) const { failed_count = 0; }

  KOKKOS_INLINE_FUNCTION
  void join(value_type &failed_count, const value_type &count) const {
    failed_count += count;
  }

  template <typename UmapOpType = InsertOp>
  KOKKOS_FORCEINLINE_FUNCTION bool is_op_noop() const {
    using vt             = typename map_type::value_type;
    using Device         = typename map_type::device_type;
    using UmapOpTypeArg1 = Kokkos::View<
        std::remove_const_t<std::conditional_t<std::is_void_v<vt>, int, vt>> *,
        Device>;
    return std::is_base_of_v<
        InsertOp, typename Kokkos::UnorderedMapInsertOpTypes<UmapOpTypeArg1,
                                                             uint32_t>::NoOp>;
  }

  template <typename UmapOpType = InsertOp>
  KOKKOS_FORCEINLINE_FUNCTION bool is_op_atomic_add() const {
    using vt             = typename map_type::value_type;
    using Device         = typename map_type::device_type;
    using UmapOpTypeArg1 = Kokkos::View<
        std::remove_const_t<std::conditional_t<std::is_void_v<vt>, int, vt>> *,
        Device>;
    return std::is_base_of_v<UmapOpType,
                             typename Kokkos::UnorderedMapInsertOpTypes<
                                 UmapOpTypeArg1, uint32_t>::AtomicAdd>;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(uint32_t i, value_type &failed_count) const {
    const uint32_t key = Near ? i / collisions : i % (inserts / collisions);
    auto ret           = map.insert(key, i, insert_op);
    if (ret.failed()) {
      ++failed_count;
      expected_values(i).map_idx = static_cast<unsigned>(~0);
    } else if (CheckValues) {
      auto map_idx                     = map.find(key);
      expected_values(map_idx).map_idx = map_idx;
      auto ptr                         = expected_values.data();
      if (is_op_atomic_add()) {
        Kokkos::atomic_add(&((ptr + map_idx)[0].v), i);
      } else if (ret.success() && is_op_noop()) {
        Kokkos::atomic_store(&((ptr + map_idx)[0].v), i);
      }
    }
  }
};

template <typename MapType, bool Near>
struct TestErase {
  using self_type = TestErase<MapType, Near>;

  using map_type        = MapType;
  using execution_space = typename MapType::execution_space;

  map_type m_map;
  uint32_t m_num_erase;
  uint32_t m_num_duplicates;

  TestErase(map_type map, uint32_t num_erases, uint32_t num_duplicates)
      : m_map(map), m_num_erase(num_erases), m_num_duplicates(num_duplicates) {}

  void testit() {
    execution_space().fence();
    Kokkos::parallel_for(m_num_erase, *this);
    execution_space().fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type i) const {
    if (Near) {
      m_map.erase(i / m_num_duplicates);
    } else {
      m_map.erase(i % (m_num_erase / m_num_duplicates));
    }
  }
};

template <typename MapType>
struct TestFind {
  using map_type        = MapType;
  using execution_space = typename MapType::execution_space::execution_space;
  using value_type      = uint32_t;

  map_type m_map;
  uint32_t m_num_insert;
  uint32_t m_num_duplicates;
  uint32_t m_max_key;

  TestFind(map_type map, uint32_t num_inserts, uint32_t num_duplicates)
      : m_map(map),
        m_num_insert(num_inserts),
        m_num_duplicates(num_duplicates),
        m_max_key(((num_inserts + num_duplicates) - 1) / num_duplicates) {}

  void testit(value_type &errors) {
    execution_space().fence();
    Kokkos::parallel_reduce(m_map.capacity(), *this, errors);
    execution_space().fence();
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &dst) { dst = 0; }

  KOKKOS_INLINE_FUNCTION
  static void join(value_type &dst, const value_type &src) { dst += src; }

  KOKKOS_INLINE_FUNCTION
  void operator()(typename execution_space::size_type i,
                  value_type &errors) const {
    const bool expect_to_find_i =
        (i < typename execution_space::size_type(m_max_key));

    const bool exists = m_map.exists(i);

    if (expect_to_find_i && !exists) ++errors;
    if (!expect_to_find_i && exists) ++errors;
  }
};

}  // namespace Impl

// MSVC reports a syntax error for this test.
// WORKAROUND MSVC
#ifndef _WIN32
template <typename Device, class map_type, class const_map_type,
          class insert_op_type, bool check_values = false>
void test_insert(uint32_t num_nodes, uint32_t num_inserts,
                 uint32_t num_duplicates, bool near) {
  const uint32_t expected_inserts =
      (num_inserts + num_duplicates - 1u) / num_duplicates;
  typename map_type::size_type arg_capacity_hint = 0;
  typename map_type::hasher_type arg_hasher;
  typename map_type::equal_to_type arg_equal_to;

  map_type map(arg_capacity_hint, arg_hasher, arg_equal_to);
  map.rehash(num_nodes, false);

  if (near) {
    Impl::TestInsert<map_type, insert_op_type, true, check_values> test_insert(
        map, num_inserts, num_duplicates);
    test_insert.testit();
  } else {
    Impl::TestInsert<map_type, insert_op_type, false, check_values> test_insert(
        map, num_inserts, num_duplicates);
    test_insert.testit();
  }

  const bool print_list = false;
  if (print_list) {
    Kokkos::Impl::UnorderedMapPrint<map_type> f(map);
    f.apply();
  }

  const uint32_t map_size = map.size();

  ASSERT_FALSE(map.failed_insert());
  {
    EXPECT_EQ(expected_inserts, map_size);

    {
      uint32_t find_errors = 0;
      Impl::TestFind<map_type> test_find(map, num_inserts, num_duplicates);
      test_find.testit(find_errors);
      EXPECT_EQ(0u, find_errors);
    }

    map.begin_erase();
    Impl::TestErase<map_type, false> test_erase(map, num_inserts,
                                                num_duplicates);
    test_erase.testit();
    map.end_erase();
    EXPECT_EQ(0u, map.size());
  }

  // Check the values from the insert operation
  {
    Impl::TestInsert<map_type, insert_op_type, true> test_insert(
        map, num_inserts, num_duplicates);
    test_insert.testit(false);
  }
}

template <typename Device>
void test_inserts(uint32_t num_nodes, uint32_t num_inserts,
                  uint32_t num_duplicates, bool near) {
  using key_type        = uint32_t;
  using value_type      = uint32_t;
  using value_view_type = Kokkos::View<value_type *, Device>;
  using size_type       = uint32_t;
  using hasher_type     = typename Kokkos::pod_hash<key_type>;
  using equal_to_type   = typename Kokkos::pod_equal_to<key_type>;

  using map_op_type =
      Kokkos::UnorderedMapInsertOpTypes<value_view_type, size_type>;
  using noop_type = typename map_op_type::NoOp;

  using map_type = Kokkos::UnorderedMap<key_type, value_type, Device,
                                        hasher_type, equal_to_type>;
  using const_map_type =
      Kokkos::UnorderedMap<const key_type, const value_type, Device,
                           hasher_type, equal_to_type>;

  test_insert<Device, map_type, const_map_type, noop_type>(
      num_nodes, num_inserts, num_duplicates, near);
}

template <typename Device>
void test_all_insert_ops(uint32_t num_nodes, uint32_t num_inserts,
                         uint32_t num_duplicates, bool near) {
  using key_type        = uint32_t;
  using value_type      = uint32_t;
  using value_view_type = Kokkos::View<value_type *, Device>;
  using size_type       = uint32_t;
  using hasher_type     = typename Kokkos::pod_hash<key_type>;
  using equal_to_type   = typename Kokkos::pod_equal_to<key_type>;

  using map_op_type =
      Kokkos::UnorderedMapInsertOpTypes<value_view_type, size_type>;
  using noop_type       = typename map_op_type::NoOp;
  using atomic_add_type = typename map_op_type::AtomicAdd;

  using map_type = Kokkos::UnorderedMap<key_type, value_type, Device,
                                        hasher_type, equal_to_type>;
  using const_map_type =
      Kokkos::UnorderedMap<const key_type, const value_type, Device,
                           hasher_type, equal_to_type>;

  test_insert<Device, map_type, const_map_type, noop_type, true>(
      num_nodes, num_inserts, num_duplicates, near);
  test_insert<Device, map_type, const_map_type, atomic_add_type, true>(
      num_nodes, num_inserts, num_duplicates, near);
}
#endif

template <typename Device>
void test_failed_insert(uint32_t num_nodes) {
  using map_type = Kokkos::UnorderedMap<uint32_t, uint32_t, Device>;

  map_type map(num_nodes);
  Impl::TestInsert<map_type> test_insert(map, 2u * num_nodes, 1u);
  test_insert.testit(false /*don't rehash on fail*/);
  typename Device::execution_space().fence();

  EXPECT_TRUE(map.failed_insert());
}

template <typename Device>
void test_deep_copy(uint32_t num_nodes) {
  using map_type = Kokkos::UnorderedMap<uint32_t, uint32_t, Device>;
  using const_map_type =
      Kokkos::UnorderedMap<const uint32_t, const uint32_t, Device>;

  using host_map_type = typename map_type::HostMirror;

  map_type map;
  map.rehash(num_nodes, false);

  {
    Impl::TestInsert<map_type> test_insert(map, num_nodes, 1);
    test_insert.testit();
    ASSERT_EQ(map.size(), num_nodes);
    ASSERT_FALSE(map.failed_insert());
    {
      uint32_t find_errors = 0;
      Impl::TestFind<map_type> test_find(map, num_nodes, 1);
      test_find.testit(find_errors);
      EXPECT_EQ(find_errors, 0u);
    }
  }

  auto hmap = create_mirror(map);
  Kokkos::deep_copy(hmap, map);

  ASSERT_EQ(map.size(), hmap.size());
  ASSERT_EQ(map.capacity(), hmap.capacity());
  {
    uint32_t find_errors = 0;
    Impl::TestFind<host_map_type> test_find(hmap, num_nodes, 1);
    test_find.testit(find_errors);
    EXPECT_EQ(find_errors, 0u);
  }

  map_type mmap;
  mmap.allocate_view(hmap);
  Kokkos::deep_copy(mmap, hmap);

  const_map_type cmap = mmap;

  EXPECT_EQ(cmap.size(), num_nodes);

  {
    uint32_t find_errors = 0;
    Impl::TestFind<const_map_type> test_find(cmap, num_nodes, 1);
    test_find.testit(find_errors);
    EXPECT_EQ(find_errors, 0u);
  }
}

#if !defined(_WIN32)
TEST(TEST_CATEGORY, UnorderedMap_insert) {
  for (int i = 0; i < 500; ++i) {
    test_inserts<TEST_EXECSPACE>(100000, 90000, 100, true);
    test_inserts<TEST_EXECSPACE>(100000, 90000, 100, false);
  }
  for (int i = 0; i < 5; ++i) {
    test_all_insert_ops<TEST_EXECSPACE>(1000, 900, 10, true);
    test_all_insert_ops<TEST_EXECSPACE>(1000, 900, 10, false);
  }
}
#endif

TEST(TEST_CATEGORY, UnorderedMap_failed_insert) {
  for (int i = 0; i < 1000; ++i) test_failed_insert<TEST_EXECSPACE>(10000);
}

TEST(TEST_CATEGORY, UnorderedMap_deep_copy) {
  for (int i = 0; i < 2; ++i) test_deep_copy<TEST_EXECSPACE>(10000);
}

TEST(TEST_CATEGORY, UnorderedMap_valid_empty) {
  using Key   = int;
  using Value = int;
  using Map   = Kokkos::UnorderedMap<Key, Value, TEST_EXECSPACE>;

  Map m{};
  Map n{};
  n = Map{m.capacity()};
  n.rehash(m.capacity());
  n.create_copy_view(m);
  ASSERT_TRUE(m.is_allocated());
  ASSERT_TRUE(n.is_allocated());
}

/**
 * This helper is needed because NVCC does not like extended lambdas
 * in private member functions.
 * Google Test bodies are private member functions. So it is incompatible.
 * See also https://github.com/google/googletest/issues/4104.
 */
template <typename map_type>
struct UnorderedMapInsert {
  //! Type of range-for policy and its index type.
  using range_policy_t =
      Kokkos::RangePolicy<typename map_type::execution_space,
                          Kokkos::IndexType<unsigned short int>>;
  using index_t = typename range_policy_t::index_type;

  const map_type m_map;

  //! Ensure shared ownership of @ref m_map.
  UnorderedMapInsert(map_type map) : m_map(std::move(map)) {}

  //! Insert a single value.
  template <typename T>
  void insert_single(const T &arg) const {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename map_type::execution_space>(0, 1),
        // NOLINTNEXTLINE(kokkos-implicit-this-capture)
        KOKKOS_CLASS_LAMBDA(const index_t) { m_map.insert(arg); });
  }

  //! Insert multiple values.
  template <typename... Args>
  void insert(Args &&... args) const {
    static_assert(sizeof...(Args) > 1, "Prefer the single value version");
    constexpr size_t size = sizeof...(Args);
    Kokkos::Array<typename map_type::key_type, size> values{
        std::forward<Args>(args)...};
    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename map_type::execution_space>(0, size),
        // NOLINTNEXTLINE(kokkos-implicit-this-capture)
        KOKKOS_CLASS_LAMBDA(const index_t i) { m_map.insert(values[i]); });
  }
};

TEST(TEST_CATEGORY, UnorderedMap_clear_zero_size) {
  using map_type = Kokkos::UnorderedMap<int, void, TEST_EXECSPACE>;

  map_type m(11);

  ASSERT_EQ(0u, m.size());

  UnorderedMapInsert<map_type>(m).insert(2, 3, 5, 7);

  ASSERT_EQ(4u, m.size());
  m.rehash(0);
  ASSERT_EQ(128u, m.capacity());
  ASSERT_EQ(4u, m.size());

  m.clear();
  ASSERT_EQ(0u, m.size());
}

TEST(TEST_CATEGORY, UnorderedMap_consistent_size) {
  using map_type = Kokkos::UnorderedMap<int, void, TEST_EXECSPACE>;

  map_type m(11);
  UnorderedMapInsert<map_type> inserter(m);

  inserter.insert_single(7);

  ASSERT_EQ(1u, m.size());

  {
    auto m_copy = m;
    UnorderedMapInsert<decltype(m_copy)> inserter_copy(m_copy);
    inserter_copy.insert_single(2);
    // This line triggers modified flags to be cleared in both m and m2
    const auto sz = m_copy.size();
    ASSERT_EQ(2u, sz);
  }

  ASSERT_EQ(2u, m.size());
}

struct TestMapCopy {
  using map_type = Kokkos::UnorderedMap<int, void, TEST_EXECSPACE>;
  map_type m_map;

  KOKKOS_FUNCTION
  void test_insert_to_map_copy(map_type const &input_map, const int i) const {
    auto map = input_map;
    map.insert(i);
  }

  KOKKOS_FUNCTION
  void operator()(const int i) const { test_insert_to_map_copy(m_map, i); }
};

TEST(TEST_CATEGORY, UnorderedMap_shallow_copyable_on_device) {
  TestMapCopy test_map_copy;

  Kokkos::parallel_for(Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1),
                       test_map_copy);
  ASSERT_EQ(1u, test_map_copy.m_map.size());
}

#if !defined(KOKKOS_ENABLE_CUDA) || \
    (defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_LAMBDA))
void test_unordered_map_device_capture() {
  TestMapCopy::map_type map;

  Kokkos::parallel_for(
      Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1),
      KOKKOS_LAMBDA(int const i) { map.insert(i); });

  ASSERT_EQ(1u, map.size());
}

TEST(TEST_CATEGORY, UnorderedMap_lambda_capturable) {
  test_unordered_map_device_capture();
}
#endif

/**
 * @test This test ensures that an @ref UnorderedMap can be built
 *       with an execution space instance (using @ref view_alloc).
 */
TEST(TEST_CATEGORY, UnorderedMap_constructor_view_alloc) {
  using map_type = Kokkos::UnorderedMap<size_t, void, TEST_EXECSPACE>;
  map_type map(Kokkos::view_alloc(TEST_EXECSPACE{}, "test umap"), 150);
  ASSERT_EQ(map.size(), 0u);
  ASSERT_GE(map.capacity(), 150u);
  ASSERT_TRUE(map.is_allocated());
}

}  // namespace Test

#endif  // KOKKOS_TEST_UNORDERED_MAP_HPP
