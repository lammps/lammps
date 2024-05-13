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

#include <Kokkos_DetectionIdiom.hpp>

void test_nonesuch() {
  using Kokkos::nonesuch;
  static_assert(!std::is_constructible<nonesuch>::value);
  static_assert(!std::is_destructible<nonesuch>::value);
  static_assert(!std::is_copy_constructible<nonesuch>::value);
  static_assert(!std::is_move_constructible<nonesuch>::value);
  static_assert(!std::is_aggregate<nonesuch>::value);
}

namespace Example {
// Example from https://en.cppreference.com/w/cpp/experimental/is_detected
template <class T>
using copy_assign_t = decltype(std::declval<T&>() = std::declval<const T&>());

struct Meow {};
struct Purr {
  void operator=(const Purr&) = delete;
};

static_assert(Kokkos::is_detected<copy_assign_t, Meow>::value,
              "Meow should be copy assignable!");
static_assert(!Kokkos::is_detected<copy_assign_t, Purr>::value,
              "Purr should not be copy assignable!");
static_assert(Kokkos::is_detected_exact<Meow&, copy_assign_t, Meow>::value,
              "Copy assignment of Meow should return Meow&!");

template <class T>
using diff_t = typename T::difference_type;

template <class Ptr>
using difference_type = Kokkos::detected_or_t<std::ptrdiff_t, diff_t, Ptr>;

struct Woof {
  using difference_type = int;
};
struct Bark {};

static_assert(std::is_same<difference_type<Woof>, int>::value,
              "Woof's difference_type should be int!");
static_assert(std::is_same<difference_type<Bark>, std::ptrdiff_t>::value,
              "Bark's difference_type should be ptrdiff_t!");
}  // namespace Example
