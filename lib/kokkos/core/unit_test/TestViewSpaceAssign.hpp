// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
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
