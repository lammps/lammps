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

#include <Kokkos_Core.hpp>
#include <sstream>
#include <iostream>

namespace Test {

template <typename SpaceDst, typename SpaceSrc>
void view_space_assign() {
  Kokkos::View<double*, SpaceDst> a = Kokkos::View<double*, SpaceSrc>("a", 1);

  Kokkos::View<double*, Kokkos::LayoutLeft, SpaceDst> b =
      Kokkos::View<double*, Kokkos::LayoutLeft, SpaceSrc>("b", 1);

  Kokkos::View<double*, Kokkos::LayoutRight, SpaceDst> c =
      Kokkos::View<double*, Kokkos::LayoutRight, SpaceSrc>("c", 1);

  Kokkos::View<double*, SpaceDst, Kokkos::MemoryRandomAccess> d =
      Kokkos::View<double*, SpaceSrc>("d", 1);

  Kokkos::View<double*, Kokkos::LayoutLeft, SpaceDst,
               Kokkos::MemoryRandomAccess>
      e = Kokkos::View<double*, Kokkos::LayoutLeft, SpaceSrc>("e", 1);

  // Rank-one layout can assign:
  Kokkos::View<double*, Kokkos::LayoutRight, SpaceDst> f =
      Kokkos::View<double*, Kokkos::LayoutLeft, SpaceSrc>("f", 1);
}

}  // namespace Test
