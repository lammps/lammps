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

#include <Kokkos_Core.hpp>

#include <type_traits>

namespace {

constexpr bool test_view_rank() {
  // clang-format off
  static_assert(Kokkos::View<int         >::rank == 0);
  static_assert(Kokkos::View<int[1]      >::rank == 1);
  static_assert(Kokkos::View<int *       >::rank == 1);
  static_assert(Kokkos::View<int[1][2]   >::rank == 2);
  static_assert(Kokkos::View<int * [1]   >::rank == 2);
  static_assert(Kokkos::View<int *  *    >::rank == 2);
  static_assert(Kokkos::View<int[1][2][3]>::rank == 3);
  static_assert(Kokkos::View<int * [1][2]>::rank == 3);
  static_assert(Kokkos::View<int *  * [1]>::rank == 3);
  static_assert(Kokkos::View<int *  *  * >::rank == 3);
  // clang-format on
  return true;
}
static_assert(test_view_rank());

constexpr bool test_is_view_type_trait() {
  struct NotView {};
  static_assert(Kokkos::is_view<Kokkos::View<int>>::value);
  static_assert(Kokkos::is_view_v<Kokkos::View<int>>);
  static_assert(!Kokkos::is_view_v<NotView>);
  static_assert(!Kokkos::is_view<NotView>::value);
  return true;
}
static_assert(test_is_view_type_trait());

}  // namespace
