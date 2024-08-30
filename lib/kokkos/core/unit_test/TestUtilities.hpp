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

#include <gtest/gtest.h>

#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {

void test_is_specialization_of() {
  using Kokkos::Impl::is_specialization_of;
  static_assert(is_specialization_of<Kokkos::pair<float, int>, Kokkos::pair>{});
  static_assert(!is_specialization_of<Kokkos::View<int*>, Kokkos::pair>{});
  static_assert(is_specialization_of<Kokkos::View<int*>, Kokkos::View>{});
  // NOTE Not removing cv-qualifiers
  static_assert(
      !is_specialization_of<Kokkos::View<int*> const, Kokkos::View>{});
  // NOTE Would not compile because Kokkos::Array takes a non-type template
  // parameter
  // static_assert(is_specialization_of<Kokkos::Array<int, 4>,
  //               Kokkos::Array>{});
  // But this is fine of course
  static_assert(!is_specialization_of<Kokkos::Array<float, 2>, Kokkos::pair>{});
}

namespace {
enum Enum { EZero, EOne };
enum EnumBool : bool { EBFalse, EBTrue };
enum class ScopedEnum { SEZero, SEOne };
enum class ScopedEnumShort : short { SESZero, SESOne };
class Class {};

template <typename Base, typename Derived>
inline constexpr bool is_public_unambiguous_base_of_v =
    std::is_convertible_v<Derived*, Base*> && !std::is_same_v<Derived, Base>;
}  // namespace

void test_to_underlying() {
  using Kokkos::Impl::to_underlying;

  constexpr auto e0 = to_underlying(EZero);
  static_assert(e0 == 0);

  constexpr auto e1 = to_underlying(EOne);
  static_assert(e1 == 1);

  constexpr auto eb0 = to_underlying(EBFalse);
  constexpr bool b0  = false;
  static_assert(std::is_same_v<decltype(eb0), decltype(b0)>);
  static_assert(eb0 == b0);

  constexpr auto eb1 = to_underlying(EBTrue);
  constexpr bool b1  = true;
  static_assert(std::is_same_v<decltype(eb1), decltype(b1)>);
  static_assert(eb1 == b1);

  constexpr auto se0 = to_underlying(ScopedEnum::SEZero);
  static_assert(se0 == 0);

  constexpr auto se1 = to_underlying(ScopedEnum::SEOne);
  static_assert(se1 == 1);

  constexpr auto ses0 = to_underlying(ScopedEnumShort::SESZero);
  constexpr short s0  = 0;
  static_assert(std::is_same_v<decltype(ses0), decltype(s0)>);
  static_assert(ses0 == s0);

  constexpr auto ses1 = to_underlying(ScopedEnumShort::SESOne);
  constexpr short s1  = 1;
  static_assert(std::is_same_v<decltype(ses1), decltype(s1)>);
  static_assert(ses1 == s1);
}

void test_is_scoped_enum() {
  using Kokkos::Impl::is_scoped_enum;
  using Kokkos::Impl::is_scoped_enum_v;

  static_assert(!is_scoped_enum<int>{});
  static_assert(!is_scoped_enum<int>::value);
  static_assert(!is_scoped_enum_v<int>);
  static_assert(
      is_public_unambiguous_base_of_v<std::false_type, is_scoped_enum<int>>);

  static_assert(!is_scoped_enum<Class>{});
  static_assert(!is_scoped_enum<Class>::value);
  static_assert(!is_scoped_enum_v<Class>);
  static_assert(
      is_public_unambiguous_base_of_v<std::false_type, is_scoped_enum<Class>>);

  static_assert(!is_scoped_enum<Enum>{});
  static_assert(!is_scoped_enum<Enum>::value);
  static_assert(!is_scoped_enum_v<Enum>);
  static_assert(
      is_public_unambiguous_base_of_v<std::false_type, is_scoped_enum<Enum>>);

  static_assert(!is_scoped_enum<EnumBool>{});
  static_assert(!is_scoped_enum<EnumBool>::value);
  static_assert(!is_scoped_enum_v<EnumBool>);
  static_assert(is_public_unambiguous_base_of_v<std::false_type,
                                                is_scoped_enum<EnumBool>>);

  static_assert(is_scoped_enum<ScopedEnum>{});
  static_assert(is_scoped_enum<ScopedEnum>::value);
  static_assert(is_scoped_enum_v<ScopedEnum>);
  static_assert(is_public_unambiguous_base_of_v<std::true_type,
                                                is_scoped_enum<ScopedEnum>>);

  static_assert(is_scoped_enum<ScopedEnumShort>{});
  static_assert(is_scoped_enum<ScopedEnumShort>::value);
  static_assert(is_scoped_enum_v<ScopedEnumShort>);
  static_assert(
      is_public_unambiguous_base_of_v<std::true_type,
                                      is_scoped_enum<ScopedEnumShort>>);
}

}  // namespace Test
