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
#include <type_traits>
#include <utility>

namespace {

template <class ExecutionSpace>
struct TestSwap {
  KOKKOS_FUNCTION void operator()(int, int& err) const {
    {
      int a = 1;
      int b = 2;
      Kokkos::kokkos_swap(a, b);
      if (!(a == 2 && b == 1)) {
        Kokkos::printf("Failed Kokkos::kokkos_swap(int, int)\n");
        ++err;
      }
    }
    {
      float a = 1;
      float b = 2;
      Kokkos::kokkos_swap(a, b);
      if (!(a == 2 && b == 1)) {
        Kokkos::printf("Failed Kokkos::kokkos_swap(float, float)\n");
        ++err;
      }
    }
    {
      int a[3] = {1, 2, 3};
      int b[3] = {4, 5, 6};
      Kokkos::kokkos_swap(a, b);
      if (!(a[0] == 4 && a[1] == 5 && a[2] == 6 && b[0] == 1 && b[1] == 2 &&
            b[2] == 3)) {
        Kokkos::printf("Failed Kokkos::kokkos_swap(int[3], int[3])\n");
        ++err;
      }
    }
  }

  TestSwap() {
    int errors;
    Kokkos::parallel_reduce(
        "TestSwap", Kokkos::RangePolicy<ExecutionSpace>(0, 1), *this, errors);
    EXPECT_EQ(errors, 0);
  }
};

TEST(TEST_CATEGORY, kokkos_swap) { TestSwap<TEST_EXECSPACE>(); }

}  // namespace
