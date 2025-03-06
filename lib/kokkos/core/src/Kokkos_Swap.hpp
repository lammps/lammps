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

#ifndef KOKKOS_SWAP_HPP
#define KOKKOS_SWAP_HPP

#include <Kokkos_Macros.hpp>

#include <cstddef>
#include <type_traits>
#include <utility>

namespace Kokkos {

template <class T>
KOKKOS_FUNCTION constexpr std::enable_if_t<std::is_move_constructible_v<T> &&
                                           std::is_move_assignable_v<T>>
kokkos_swap(T& a, T& b) noexcept(std::is_nothrow_move_constructible_v<T>&&
                                     std::is_nothrow_move_assignable_v<T>) {
  T t(std::move(a));
  a = std::move(b);
  b = std::move(t);
}

namespace Impl {

template <class T>
struct is_swappable {
  template <class U>
  static decltype(kokkos_swap(std::declval<T&>(), std::declval<T&>()))
  test_swap(int);
  struct Nope;
  template <class U>
  static Nope test_swap(long);
  static constexpr bool value =
      !std::is_same_v<decltype(test_swap<T>(0)), Nope>;
};

template <class T>
inline constexpr bool is_nothrow_swappable_v =
    noexcept(kokkos_swap(std::declval<T&>(), std::declval<T&>()));

}  // namespace Impl

template <class T, std::size_t N>
KOKKOS_FUNCTION constexpr std::enable_if_t<Impl::is_swappable<T>::value>
kokkos_swap(T (&a)[N], T (&b)[N]) noexcept(Impl::is_nothrow_swappable_v<T>) {
  for (std::size_t i = 0; i < N; ++i) {
    kokkos_swap(a[i], b[i]);
  }
}

}  // namespace Kokkos

#endif
