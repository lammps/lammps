// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
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
