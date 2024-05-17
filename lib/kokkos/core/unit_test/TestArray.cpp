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

#include <Kokkos_Array.hpp>

namespace {

KOKKOS_FUNCTION constexpr bool test_array() {
  constexpr Kokkos::Array<int, 3> a{{1, 2}};

  static_assert(!a.empty());
  static_assert(a.size() == 3);
  static_assert(a.max_size() == 3);

  static_assert(*a.data() == 1);
  static_assert(a[1] == 2);

  return true;
}

static_assert(test_array());

KOKKOS_FUNCTION constexpr bool test_array_structured_binding_support() {
  constexpr Kokkos::Array<float, 2> a{};
  auto& [xr, yr] = a;
  (void)xr;
  (void)yr;
  auto [x, y] = a;
  (void)x;
  (void)y;
  auto const& [xcr, ycr] = a;
  (void)xcr;
  (void)ycr;
  return true;
}

static_assert(test_array_structured_binding_support());

template <typename L, typename R>
KOKKOS_FUNCTION constexpr bool is_equal(L const& l, R const& r) {
  if (std::size(l) != std::size(r)) return false;

  for (size_t i = 0; i != std::size(l); ++i) {
    if (l[i] != r[i]) return false;
  }

  return true;
}

// Disable ctad test for intel versions < 2021, see issue #6702
#if !defined(KOKKOS_COMPILER_INTEL) || KOKKOS_COMPILER_INTEL >= 2021
KOKKOS_FUNCTION constexpr bool test_array_ctad() {
  constexpr int x = 10;
  constexpr Kokkos::Array a{1, 2, 3, 5, x};
  constexpr Kokkos::Array<int, 5> b{1, 2, 3, 5, x};

  return std::is_same_v<decltype(a), decltype(b)> && is_equal(a, b);
}

static_assert(test_array_ctad());
#endif

}  // namespace
