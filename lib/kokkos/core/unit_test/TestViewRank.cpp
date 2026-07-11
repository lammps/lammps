// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <cstddef>
#include <type_traits>

namespace {

template <class View, size_t Rank, size_t RankDynamic>
constexpr bool test_view_rank_and_dynamic_rank() {
  static_assert(View::rank == Rank);
  static_assert(View::rank() == Rank);
  static_assert(View::rank_dynamic == RankDynamic);
  static_assert(View::rank_dynamic() == RankDynamic);
  static_assert(std::is_convertible_v<decltype(View::rank), size_t>);
  static_assert(std::is_same_v<decltype(View::rank()), size_t>);
  static_assert(std::is_convertible_v<decltype(View::rank_dynamic), size_t>);
  static_assert(std::is_same_v<decltype(View::rank_dynamic()), size_t>);
  auto rank = View::rank;  // not an integral type in contrast to Kokkos version
  // less than 4.0.01
  static_assert(!std::is_integral_v<decltype(rank)>);
  auto rank_preferred = View::rank();  // since 4.0.01
  static_assert(std::is_same_v<decltype(rank_preferred), size_t>);
  return rank == rank_preferred;
}

// clang-format off
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<long long>, 0, 0>());

static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<unsigned[1]>, 1, 0>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<unsigned * >, 1, 1>());

static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<double[1][2]>, 2, 0>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<double * [2]>, 2, 1>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<double *  * >, 2, 2>());

static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<float[1][2][3]>, 3, 0>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<float * [2][3]>, 3, 1>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<float *  * [3]>, 3, 2>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<float *  *  * >, 3, 3>());

static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<int[1][2][3][4]>, 4, 0>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<int * [2][3][4]>, 4, 1>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<int *  * [3][4]>, 4, 2>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<int *  *  * [4]>, 4, 3>());
static_assert(test_view_rank_and_dynamic_rank<Kokkos::View<int *  *  *  * >, 4, 4>());
//clang-format on

}  // namespace
