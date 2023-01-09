/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
        auto mirror_host =
            Kokkos::create_mirror(Kokkos::WithoutInitializing,
                                  Kokkos::DefaultExecutionSpace{}, host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, Kokkos::DefaultExecutionSpace{},
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
                               Kokkos::DefaultExecutionSpace{}),
            host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::HostSpace{},
                               Kokkos::WithoutInitializing),
            device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::HostSpace{}, Kokkos::WithoutInitializing,
                               Kokkos::DefaultExecutionSpace{}),
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
