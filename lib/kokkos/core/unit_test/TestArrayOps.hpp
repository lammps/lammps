// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <numeric>

namespace {

TEST(TEST_CATEGORY, array_capacity) {
  using A = Kokkos::Array<int, 2>;
  A a{{3, 5}};

  static_assert(noexcept(a.empty()));
  ASSERT_FALSE(a.empty());

  static_assert(noexcept(a.size()));
  ASSERT_EQ(a.size(), 2u);

  static_assert(noexcept(a.max_size()));
  ASSERT_EQ(a.max_size(), 2u);
}

enum Enum { EZero, EOne };
enum EnumShort : short { ESZero, ESOne };

TEST(TEST_CATEGORY, array_element_access) {
  using A = Kokkos::Array<int, 2>;
  A a{{3, 5}};
  A const& ca = a;

  size_t index = 1;
  ASSERT_EQ(a[index], 5);

  auto sc = static_cast<signed char>(index);
  ASSERT_EQ(a[sc], a[index]);
  ASSERT_EQ(ca[sc], a[index]);

  auto uc = static_cast<unsigned char>(index);
  ASSERT_EQ(a[uc], a[index]);
  ASSERT_EQ(ca[uc], a[index]);

  auto s = static_cast<short>(index);
  ASSERT_EQ(a[s], a[index]);
  ASSERT_EQ(ca[s], a[index]);

  auto us = static_cast<unsigned short>(index);
  ASSERT_EQ(a[us], a[index]);
  ASSERT_EQ(ca[us], a[index]);

  auto i = static_cast<int>(index);
  ASSERT_EQ(a[i], a[index]);
  ASSERT_EQ(ca[i], a[index]);

  auto ui = static_cast<unsigned int>(index);
  ASSERT_EQ(a[ui], a[index]);
  ASSERT_EQ(ca[ui], a[index]);

  auto l = static_cast<long>(index);
  ASSERT_EQ(a[l], a[index]);
  ASSERT_EQ(ca[l], a[index]);

  auto ul = static_cast<unsigned long>(index);
  ASSERT_EQ(a[ul], a[index]);
  ASSERT_EQ(ca[ul], a[index]);

  auto ll = static_cast<long long>(index);
  ASSERT_EQ(a[ll], a[index]);
  ASSERT_EQ(ca[ll], a[index]);

  auto ull = static_cast<unsigned long long>(index);
  ASSERT_EQ(a[ull], a[index]);
  ASSERT_EQ(ca[ull], a[index]);

  auto e = static_cast<Enum>(index);
  ASSERT_EQ(a[e], a[index]);
  ASSERT_EQ(ca[e], a[index]);

  auto es = static_cast<EnumShort>(index);
  ASSERT_EQ(a[es], a[index]);
  ASSERT_EQ(ca[es], a[index]);

  static_assert(noexcept(a.data()));
  ASSERT_EQ(a.data()[index], a[index]);

  static_assert(noexcept(ca.data()));
  ASSERT_EQ(ca.data()[index], a[index]);
}

TEST(TEST_CATEGORY, array_operator_equal) {
  using A = Kokkos::Array<int, 2>;
  constexpr A a{{3, 5}};
  constexpr A b{{3, 5}};
  constexpr A c{{5, 3}};

  static_assert(a == b);
  static_assert(!(a == c));
  static_assert(a != c);

  ASSERT_TRUE(a == b);
  ASSERT_FALSE(a == c);
  ASSERT_TRUE(a != c);

  using E = Kokkos::Array<int, 0>;
  constexpr E e;
  constexpr E f;

  static_assert(e == f);
  static_assert(!(e != f));

  ASSERT_TRUE(e == f);
  ASSERT_FALSE(e != f);
}

TEST(TEST_CATEGORY, array_zero_capacity) {
  using A = Kokkos::Array<int, 0>;
  A e;

  static_assert(noexcept(e.empty()));
  ASSERT_TRUE(e.empty());

  static_assert(noexcept(e.size()));
  ASSERT_EQ(e.size(), 0u);

  static_assert(noexcept(e.max_size()));
  ASSERT_EQ(e.max_size(), 0u);
}

TEST(TEST_CATEGORY, array_zero_data_nullptr) {
  using A = Kokkos::Array<int, 0>;

  A e;
  static_assert(noexcept(e.data()));
  ASSERT_EQ(e.data(), nullptr);

  const A& ce = e;
  static_assert(noexcept(ce.data()));
  ASSERT_EQ(ce.data(), nullptr);
}

TEST(TEST_CATEGORY, array_stl_compatibility) {
  using A = Kokkos::Array<int, 3>;
  A a{1, 2, 3};

  int sum  = std::reduce(a.begin(), a.end(), int{0}, std::plus<int>{});
  int csum = std::reduce(a.cbegin(), a.cend(), int{0}, std::plus<int>{});
  int expected_sum = 1 + 2 + 3;
  EXPECT_EQ(sum, expected_sum);
  EXPECT_EQ(csum, expected_sum);
}

}  // namespace
