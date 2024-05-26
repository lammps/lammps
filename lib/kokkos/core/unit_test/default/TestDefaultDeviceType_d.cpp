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

#if !defined(KOKKOS_ENABLE_CUDA) || defined(__CUDACC__)

#include <TestDefaultDeviceType_Category.hpp>
#include <TestUtilities.hpp>

namespace Test {

TEST(defaultdevicetype, malloc) {
  int* data = static_cast<int*>(Kokkos::kokkos_malloc(100 * sizeof(int)));
  ASSERT_NO_THROW(data = static_cast<int*>(
                      Kokkos::kokkos_realloc(data, 120 * sizeof(int))));
  Kokkos::kokkos_free(data);

  int* data2 = static_cast<int*>(Kokkos::kokkos_malloc(0));
  ASSERT_EQ(data2, nullptr);
  Kokkos::kokkos_free(data2);
}

}  // namespace Test

#endif
