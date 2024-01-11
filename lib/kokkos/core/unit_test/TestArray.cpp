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

}  // namespace
