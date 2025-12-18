// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <impl/Kokkos_Utilities.hpp>

using type_list_empty_t = Kokkos::Impl::type_list<>;

using TypeList2 = Kokkos::Impl::type_list<void, bool>;
using TypeList3 = Kokkos::Impl::type_list<char, short, int>;
using TypeList223 =
    Kokkos::Impl::type_list<void, bool, void, bool, char, short, int>;
using TypeList223Void   = Kokkos::Impl::type_list<void, void>;
using TypeList223NoVoid = Kokkos::Impl::type_list<bool, bool, char, short, int>;

// concat_type_list
using ConcatTypeList2 = Kokkos::Impl::concat_type_list_t<TypeList2>;
static_assert(std::is_same_v<TypeList2, ConcatTypeList2>,
              "concat_type_list of a single type_list failed");

using ConcatTypeList223 =
    Kokkos::Impl::concat_type_list_t<TypeList2, TypeList2, TypeList3>;
static_assert(std::is_same_v<TypeList223, ConcatTypeList223>,
              "concat_type_list of three type_lists failed");

// filter_type_list
using FilterTypeList223Void =
    Kokkos::Impl::filter_type_list_t<std::is_void, TypeList223>;
static_assert(std::is_same_v<TypeList223Void, FilterTypeList223Void>,
              "filter_type_list with predicate value==true failed");

using FilterTypeList223NoVoid =
    Kokkos::Impl::filter_type_list_t<std::is_void, TypeList223, false>;
static_assert(std::is_same_v<TypeList223NoVoid, FilterTypeList223NoVoid>,
              "filter_type_list with predicate value==false failed");

constexpr bool test_type_list_any() {
  using Kokkos::Impl::type_list;
  using Kokkos::Impl::type_list_any_v;

  static_assert(!type_list_any_v<std::is_enum, type_list_empty_t>);
  static_assert(type_list_any_v<std::is_floating_point,
                                type_list<float, double, char, int>>);
  static_assert(type_list_any_v<std::is_integral,
                                type_list<float, char, int, std::size_t>>);
  static_assert(!type_list_any_v<std::is_enum, type_list<float, char, int>>);

  return true;
}
static_assert(test_type_list_any());

constexpr bool test_type_list_size() {
  using Kokkos::Impl::type_list_size_v;

  static_assert(type_list_size_v<type_list_empty_t> == 0);
  static_assert(type_list_size_v<TypeList2> == 2);
  static_assert(type_list_size_v<TypeList3> == 3);

  return true;
}
static_assert(test_type_list_size());

constexpr bool test_type_list_contains() {
  using Kokkos::Impl::type_list_contains_v;

  static_assert(!type_list_contains_v<int, type_list_empty_t>);

  static_assert(type_list_contains_v<char, TypeList3>);
  static_assert(type_list_contains_v<short, TypeList3>);
  static_assert(type_list_contains_v<int, TypeList3>);
  static_assert(!type_list_contains_v<double, TypeList3>);

  return true;
}
static_assert(test_type_list_contains());
