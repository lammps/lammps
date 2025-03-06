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
#include <TestDefaultDeviceType_Category.hpp>
#include <TestHalfConversion.hpp>
#include <TestHalfOperators.hpp>

#if !defined(KOKKOS_ENABLE_CUDA) || defined(__CUDACC__)

namespace Test {

TEST(TEST_CATEGORY, host_space_access) {
  using host_exec_space = Kokkos::HostSpace::execution_space;
  using device_space    = Kokkos::Device<host_exec_space, Kokkos::HostSpace>;
  using mirror_space =
      Kokkos::Impl::HostMirror<Kokkos::DefaultExecutionSpace>::Space;

  static_assert(Kokkos::SpaceAccessibility<host_exec_space,
                                           Kokkos::HostSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<device_space, Kokkos::HostSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<mirror_space, Kokkos::HostSpace>::accessible);
}

}  // namespace Test

#endif
