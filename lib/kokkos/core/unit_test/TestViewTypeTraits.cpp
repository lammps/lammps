// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

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
