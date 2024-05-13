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
