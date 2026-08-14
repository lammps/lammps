// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif
#include <cstddef>

// Test View equality operators
// View equivalence defined as matching value_type,
// layout_type, memory_space, rank, span, data pointers and extents

namespace {

template <class T, class Enable = void>
struct Another;

template <class T>
struct Another<T, std::enable_if_t<Kokkos::is_array_layout_v<T>>> {
  using type = std::conditional_t<std::is_same_v<T, Kokkos::LayoutLeft>,
                                  Kokkos::LayoutRight, Kokkos::LayoutLeft>;
};

template <class T>
struct Another<T, std::enable_if_t<Kokkos::is_memory_traits_v<T>>> {
  using type =
      std::conditional_t<std::is_same_v<T, Kokkos::MemoryRandomAccess>,
                         Kokkos::MemoryUnmanaged, Kokkos::MemoryRandomAccess>;
};

template <class T>
struct Another<T, std::enable_if_t<Kokkos::is_memory_space_v<T>>> {
  using type = Kokkos::DefaultHostExecutionSpace;
};

template <class Left, class Right>
bool check_equal(Left l, Right r) {
  return l == r && !(l != r);
}

void test_view_equality_operator() {
  using T = double;
  using V = Kokkos::View<T, TEST_EXECSPACE>;

  {
    // Check for static properties
    using V_data_type_0 = Kokkos::View<V::data_type*, TEST_EXECSPACE>;
    using V_data_type_1 = Kokkos::View<V::value_type*, TEST_EXECSPACE>;
    using V_data_type_2 = Kokkos::View<V::value_type[1], TEST_EXECSPACE>;
    using V_data_type_3 = Kokkos::View<const V::value_type, TEST_EXECSPACE>;
    using V_data_type_4 = Kokkos::View<const V::data_type, TEST_EXECSPACE>;

    ASSERT_FALSE(check_equal(V(), V_data_type_0()));
    ASSERT_FALSE(check_equal(V(), V_data_type_1()));
    ASSERT_FALSE(check_equal(V(), V_data_type_2()));
    // Note: We do not enforce matching const'ness of data_type or value_type
    ASSERT_TRUE(check_equal(V(), V_data_type_3()));
    ASSERT_TRUE(check_equal(V(), V_data_type_4()));
  }

  {
    // Check for additional static properties
    using V_memory_traits =
        Kokkos::View<T, V::memory_space, Another<V::memory_traits>::type>;
    using V_layout_type =
        Kokkos::View<T, Another<V::array_layout>::type, TEST_EXECSPACE>;
    using V_memory_space_type = Kokkos::View<T, Another<V::memory_space>::type>;

    // Note: We do not enforce Traits::memory_traits equality
    ASSERT_TRUE(check_equal(V(), V_memory_traits()));
    ASSERT_FALSE(check_equal(V(), V_layout_type()));

    auto v_mem_space = V_memory_space_type();
    constexpr bool is_default_host_space =
        std::is_same_v<TEST_EXECSPACE::memory_space,
                       Kokkos::DefaultHostExecutionSpace::memory_space>;
    ASSERT_EQ(check_equal(V(), v_mem_space), is_default_host_space);
  }

  {
    // Check for pointer equality
    ASSERT_TRUE(check_equal(V(), V()));  // nullptr is equal
    using V_t = Kokkos::View<T*, TEST_EXECSPACE>;
    using V_unmanaged_t =
        Kokkos::View<V_t::data_type, TEST_EXECSPACE, Kokkos::MemoryUnmanaged>;
    auto v_0 = V_t("v_1D_0", 3);
    auto v_1 = V_t("v_1D_1", 3);
    ASSERT_FALSE(check_equal(v_0, v_1));
    ASSERT_TRUE(check_equal(v_0, V_unmanaged_t(v_0.data(), 3)));
  }

  {
    // Check for matching static extents (same span, same ptr)
    auto v_0 = Kokkos::View<T[1][3], TEST_EXECSPACE>("v_2D_0");
    auto v_1 = Kokkos::View<T[3][1], TEST_EXECSPACE, Kokkos::MemoryUnmanaged>(
        v_0.data());
    ASSERT_FALSE(check_equal(v_0, v_1));
  }

  {
    // Check for matching dynamic extents (same span, same ptr)
    using v_t = Kokkos::View<T**, TEST_EXECSPACE>;
    auto v_0  = v_t("v_2D_2", 3, 3);
    auto v_1  = Kokkos::View<T[3][3], TEST_EXECSPACE, Kokkos::MemoryUnmanaged>(
        v_0.data());
    ASSERT_TRUE(check_equal(v_0, v_1));
  }
}

TEST(TEST_CATEGORY, view_equality_operator) { test_view_equality_operator(); }

}  // namespace
