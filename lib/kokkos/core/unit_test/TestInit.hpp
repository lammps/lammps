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

#include <cstdio>
#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {
TEST(TEST_CATEGORY, init) { ; }

#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA

template <class ExecSpace>
void test_dispatch() {
  const int repeat = 100;
  for (int i = 0; i < repeat; ++i) {
    for (int j = 0; j < repeat; ++j) {
      Kokkos::parallel_for(Kokkos::RangePolicy<TEST_EXECSPACE>(0, j),
                           KOKKOS_LAMBDA(int){});
    }
  }
}

TEST(TEST_CATEGORY, dispatch) { test_dispatch<TEST_EXECSPACE>(); }
#endif

}  // namespace Test
