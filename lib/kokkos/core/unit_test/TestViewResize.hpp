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
#ifndef TESTVIEWRESIZE_HPP_
#define TESTVIEWRESIZE_HPP_

#include <gtest/gtest.h>
#include "TestResize.hpp"
#include "TestRealloc.hpp"

namespace Test {

TEST(TEST_CATEGORY, view_resize) {
  using ExecSpace = TEST_EXECSPACE;
  TestViewResize::testResize<ExecSpace>();
}

TEST(TEST_CATEGORY, view_realloc) {
  using ExecSpace = TEST_EXECSPACE;
  TestViewRealloc::testRealloc<ExecSpace>();
}

}  // namespace Test
#endif  // TESTVIEWRESIZE_HPP_
