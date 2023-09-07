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
  static_assert(is_specialization_of<Kokkos::pair<float, int>, Kokkos::pair>{},
                "");
  static_assert(!is_specialization_of<Kokkos::View<int*>, Kokkos::pair>{}, "");
  static_assert(is_specialization_of<Kokkos::View<int*>, Kokkos::View>{}, "");
  // NOTE Not removing cv-qualifiers
  static_assert(!is_specialization_of<Kokkos::View<int*> const, Kokkos::View>{},
                "");
  // NOTE Would not compile because Kokkos::Array takes a non-type template
  // parameter
  // static_assert(is_specialization_of<Kokkos::Array<int, 4>, Kokkos::Array>{},
  // "");
  // But this is fine of course
  static_assert(!is_specialization_of<Kokkos::Array<float, 2>, Kokkos::pair>{},
                "");
}

}  // namespace Test
