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
#include "TestResize.hpp"

namespace Test {

TEST(kokkosresize, host_space_access) {
  // Test with the default device type.
  using TestViewResize::testResize;
  using device_type = Kokkos::View<int *>::device_type;
  testResize<device_type>();
}

}  // namespace Test
