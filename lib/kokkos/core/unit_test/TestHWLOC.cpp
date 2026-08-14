// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <iostream>

#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_hwloc.hpp>

namespace Test {

TEST(hwloc, query) {
  std::cout << " NUMA[" << Kokkos::hwloc::get_available_numa_count() << "]"
            << " CORE[" << Kokkos::hwloc::get_available_cores_per_numa() << "]"
            << " PU[" << Kokkos::hwloc::get_available_threads_per_core() << "]"
            << std::endl;
}

}  // namespace Test
