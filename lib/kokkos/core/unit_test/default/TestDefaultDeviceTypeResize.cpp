// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
