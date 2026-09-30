// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
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
