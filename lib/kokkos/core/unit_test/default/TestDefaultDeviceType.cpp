// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif
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
