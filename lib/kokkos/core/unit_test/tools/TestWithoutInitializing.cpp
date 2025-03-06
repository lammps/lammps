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

#include "include/ToolTestingUtilities.hpp"

TEST(kokkosp, create_mirror_no_init) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::View<int*, Kokkos::DefaultExecutionSpace> device_view("device view",
                                                                10);
  Kokkos::View<int*, Kokkos::HostSpace> host_view("host view", 10);

  auto success = validate_absence(
      [&]() {
        auto mirror_device =
            Kokkos::create_mirror(Kokkos::WithoutInitializing, device_view);
        auto mirror_host = Kokkos::create_mirror(
            Kokkos::WithoutInitializing, Kokkos::DefaultHostExecutionSpace{},
            host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, Kokkos::DefaultHostExecutionSpace{},
            host_view);
      },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found begin event"}};
      },
      [&](EndParallelForEvent) {
        return MatchDiagnostic{true, {"Found end event"}};
      });
  ASSERT_TRUE(success);
}

TEST(kokkosp, create_mirror_no_init_view_ctor) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::View<int*, Kokkos::DefaultExecutionSpace> device_view("device view",
                                                                10);
  Kokkos::View<int*, Kokkos::HostSpace> host_view("host view", 10);

  auto success = validate_absence(
      [&]() {
        auto mirror_device = Kokkos::create_mirror(
            Kokkos::view_alloc(Kokkos::HostSpace{},
                               Kokkos::WithoutInitializing),
            device_view);
        auto mirror_host = Kokkos::create_mirror(
            Kokkos::view_alloc(Kokkos::HostSpace{}, Kokkos::WithoutInitializing,
                               Kokkos::DefaultHostExecutionSpace{}),
            host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::HostSpace{},
                               Kokkos::WithoutInitializing),
            device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::HostSpace{}, Kokkos::WithoutInitializing,
                               Kokkos::DefaultHostExecutionSpace{}),
            host_view);
        mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::WithoutInitializing), host_view);
      },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found begin event"}};
      },
      [&](EndParallelForEvent) {
        return MatchDiagnostic{true, {"Found end event"}};
      });
  ASSERT_TRUE(success);
}

TEST(kokkosp, create_mirror_view_and_copy) {
#ifdef KOKKOS_ENABLE_OPENMPTARGET  // FIXME_OPENMPTARGET
  if (std::is_same<Kokkos::DefaultExecutionSpace,
                   Kokkos::Experimental::OpenMPTarget>::value)
    GTEST_SKIP() << "skipping since the OpenMPTarget has unexpected fences";
#endif

#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<Kokkos::DefaultExecutionSpace::memory_space,
                   Kokkos::CudaUVMSpace>::value)
    GTEST_SKIP()
        << "skipping since the CudaUVMSpace requires additional fences";
#endif

  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels(),
                     Config::EnableFences());
  Kokkos::View<int*, Kokkos::DefaultExecutionSpace> device_view;
  Kokkos::View<int*, Kokkos::HostSpace> host_view("host view", 10);

  auto success = validate_absence(
      [&]() {
        auto mirror_device = Kokkos::create_mirror_view_and_copy(
            Kokkos::view_alloc(
                Kokkos::DefaultExecutionSpace{},
                typename Kokkos::DefaultExecutionSpace::memory_space{}),
            host_view);
        // Avoid fences for deallocation when mirror_device goes out of scope.
        device_view = mirror_device;
      },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found parallel_for event"}};
      },
      [&](BeginFenceEvent) {
        return MatchDiagnostic{true, {"Found fence event"}};
      });
  ASSERT_TRUE(success);
}
