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
#include <Kokkos_DualView.hpp>
#include <Kokkos_DynamicView.hpp>
#include <Kokkos_DynRankView.hpp>
#include <Kokkos_OffsetView.hpp>
#include <Kokkos_ScatterView.hpp>

#include <../../core/unit_test/tools/include/ToolTestingUtilities.hpp>

TEST(TEST_CATEGORY, resize_realloc_no_init_dualview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::DualView<int*** * [1][2][3][4], TEST_EXECSPACE> bla("bla", 5, 6, 7,
                                                              8);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(Kokkos::WithoutInitializing, bla, 5, 6, 7, 9);
        EXPECT_EQ(bla.template view<TEST_EXECSPACE>().label(), "bla");
        Kokkos::realloc(Kokkos::WithoutInitializing, bla, 8, 8, 8, 8);
        EXPECT_EQ(bla.template view<TEST_EXECSPACE>().label(), "bla");
        Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), bla, 5,
                        6, 7, 8);
        EXPECT_EQ(bla.template view<TEST_EXECSPACE>().label(), "bla");
      },
      [&](BeginParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found begin event"}};
        return MatchDiagnostic{false};
      },
      [&](EndParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found end event"}};
        return MatchDiagnostic{false};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, resize_realloc_no_alloc_dualview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels(),
                     Config::EnableAllocs());
  Kokkos::DualView<int*** * [1][2][3][4], TEST_EXECSPACE> bla("bla", 8, 7, 6,
                                                              5);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(bla, 8, 7, 6, 5);
        EXPECT_EQ(bla.template view<TEST_EXECSPACE>().label(), "bla");
        Kokkos::realloc(Kokkos::WithoutInitializing, bla, 8, 7, 6, 5);
        EXPECT_EQ(bla.template view<TEST_EXECSPACE>().label(), "bla");
      },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found begin event"}};
      },
      [&](EndParallelForEvent) {
        return MatchDiagnostic{true, {"Found end event"}};
      },
      [&](AllocateDataEvent) {
        return MatchDiagnostic{true, {"Found alloc event"}};
      },
      [&](DeallocateDataEvent) {
        return MatchDiagnostic{true, {"Found dealloc event"}};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, resize_exec_space_dualview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableFences(),
                     Config::EnableKernels());
  Kokkos::DualView<int*** * [1][2][3][4], TEST_EXECSPACE> bla("bla", 8, 7, 6,
                                                              5);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(
            Kokkos::view_alloc(TEST_EXECSPACE{}, Kokkos::WithoutInitializing),
            bla, 5, 6, 7, 8);
        EXPECT_EQ(bla.template view<TEST_EXECSPACE>().label(), "bla");
      },
      [&](BeginFenceEvent event) {
        if (event.descriptor().find("Kokkos::resize(View)") !=
            std::string::npos)
          return MatchDiagnostic{true, {"Found begin event"}};
        return MatchDiagnostic{false};
      },
      [&](EndFenceEvent event) {
        if (event.descriptor().find("Kokkos::resize(View)") !=
            std::string::npos)
          return MatchDiagnostic{true, {"Found end event"}};
        return MatchDiagnostic{false};
      },
      [&](BeginParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found begin event"}};
        return MatchDiagnostic{false};
      },
      [&](EndParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found end event"}};
        return MatchDiagnostic{false};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, realloc_exec_space_dualview) {
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<typename TEST_EXECSPACE::memory_space,
                   Kokkos::CudaUVMSpace>::value)
    GTEST_SKIP() << "skipping since CudaUVMSpace requires additional fences";
#endif

  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableFences());
  using view_type = Kokkos::DualView<int*, TEST_EXECSPACE>;
  view_type v(Kokkos::view_alloc(TEST_EXECSPACE{}, "bla"), 8);

  auto success = validate_absence(
      [&]() {
        Kokkos::realloc(Kokkos::view_alloc(TEST_EXECSPACE{}), v, 8);
        EXPECT_EQ(v.template view<TEST_EXECSPACE>().label(), "bla");
      },
      [&](BeginFenceEvent event) {
        if ((event.descriptor().find("Debug Only Check for Execution Error") !=
             std::string::npos) ||
            (event.descriptor().find("HostSpace fence") != std::string::npos))
          return MatchDiagnostic{false};
        return MatchDiagnostic{true, {"Found fence event!"}};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, resize_realloc_no_init_dynrankview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::DynRankView<int, TEST_EXECSPACE> bla("bla", 5, 6, 7, 8);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(Kokkos::WithoutInitializing, bla, 5, 6, 7, 9);
        EXPECT_EQ(bla.label(), "bla");
        Kokkos::realloc(Kokkos::WithoutInitializing, bla, 8, 8, 8, 8);
        EXPECT_EQ(bla.label(), "bla");
        Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), bla, 5,
                        6, 7, 8);
        EXPECT_EQ(bla.label(), "bla");
      },
      [&](BeginParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found begin event"}};
        return MatchDiagnostic{false};
      },
      [&](EndParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found end event"}};
        return MatchDiagnostic{false};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, resize_exec_space_dynrankview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableFences(),
                     Config::EnableKernels());
  Kokkos::DynRankView<int, TEST_EXECSPACE> bla("bla", 8, 7, 6, 5);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(
            Kokkos::view_alloc(TEST_EXECSPACE{}, Kokkos::WithoutInitializing),
            bla, 5, 6, 7, 8);
        EXPECT_EQ(bla.label(), "bla");
      },
      [&](BeginFenceEvent event) {
        if (event.descriptor().find("Kokkos::resize(View)") !=
            std::string::npos)
          return MatchDiagnostic{true, {"Found begin event"}};
        return MatchDiagnostic{false};
      },
      [&](EndFenceEvent event) {
        if (event.descriptor().find("Kokkos::resize(View)") !=
            std::string::npos)
          return MatchDiagnostic{true, {"Found end event"}};
        return MatchDiagnostic{false};
      },
      [&](BeginParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found begin event"}};
        return MatchDiagnostic{false};
      },
      [&](EndParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found end event"}};
        return MatchDiagnostic{false};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, realloc_exec_space_dynrankview) {
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<typename TEST_EXECSPACE::memory_space,
                   Kokkos::CudaUVMSpace>::value)
    GTEST_SKIP() << "skipping since CudaUVMSpace requires additional fences";
#endif
// FIXME_THREADS The Threads backend fences every parallel_for
#ifdef KOKKOS_ENABLE_THREADS
  if (std::is_same<TEST_EXECSPACE, Kokkos::Threads>::value)
    GTEST_SKIP() << "skipping since the Threads backend isn't asynchronous";
#endif

  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableFences());
  using view_type = Kokkos::DynRankView<int, TEST_EXECSPACE>;
  view_type outer_view, outer_view2;

  auto success = validate_absence(
      [&]() {
        view_type inner_view(Kokkos::view_alloc(TEST_EXECSPACE{}, "bla"), 8);
        // Avoid testing the destructor
        outer_view = inner_view;
        Kokkos::realloc(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, TEST_EXECSPACE{}),
            inner_view, 10);
        EXPECT_EQ(inner_view.label(), "bla");
        outer_view2 = inner_view;
      },
      [&](BeginFenceEvent event) {
        if ((event.descriptor().find("Debug Only Check for Execution Error") !=
             std::string::npos) ||
            (event.descriptor().find("HostSpace fence") != std::string::npos))
          return MatchDiagnostic{false};
        return MatchDiagnostic{true, {"Found fence event!"}};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, resize_realloc_no_init_scatterview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::Experimental::ScatterView<
      int*** * [1][2][3], typename TEST_EXECSPACE::array_layout, TEST_EXECSPACE>
      bla("bla", 4, 5, 6, 7);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(Kokkos::WithoutInitializing, bla, 4, 5, 6, 8);
        EXPECT_EQ(bla.subview().label(), "bla");
        Kokkos::realloc(Kokkos::WithoutInitializing, bla, 8, 8, 8, 8);
        EXPECT_EQ(bla.subview().label(), "bla");
        Kokkos::realloc(Kokkos::view_alloc(Kokkos::WithoutInitializing), bla, 5,
                        6, 7, 8);
        EXPECT_EQ(bla.subview().label(), "bla");
      },
      [&](BeginParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found begin event"}};
        return MatchDiagnostic{false};
      },
      [&](EndParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found end event"}};
        return MatchDiagnostic{false};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, resize_realloc_no_alloc_scatterview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels(),
                     Config::EnableAllocs());
  Kokkos::Experimental::ScatterView<
      int*** * [1][2][3], typename TEST_EXECSPACE::array_layout, TEST_EXECSPACE>
      bla("bla", 7, 6, 5, 4);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(bla, 7, 6, 5, 4);
        EXPECT_EQ(bla.subview().label(), "bla");
        Kokkos::realloc(Kokkos::WithoutInitializing, bla, 7, 6, 5, 4);
        EXPECT_EQ(bla.subview().label(), "bla");
      },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found begin event"}};
      },
      [&](EndParallelForEvent) {
        return MatchDiagnostic{true, {"Found end event"}};
      },
      [&](AllocateDataEvent) {
        return MatchDiagnostic{true, {"Found alloc event"}};
      },
      [&](DeallocateDataEvent) {
        return MatchDiagnostic{true, {"Found dealloc event"}};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, resize_exec_space_scatterview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableFences(),
                     Config::EnableKernels());
  Kokkos::Experimental::ScatterView<
      int*** * [1][2][3], typename TEST_EXECSPACE::array_layout, TEST_EXECSPACE>
      bla("bla", 7, 6, 5, 4);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(
            Kokkos::view_alloc(TEST_EXECSPACE{}, Kokkos::WithoutInitializing),
            bla, 5, 6, 7, 8);
        EXPECT_EQ(bla.subview().label(), "bla");
      },
      [&](BeginFenceEvent event) {
        if (event.descriptor().find("Kokkos::resize(View)") !=
            std::string::npos)
          return MatchDiagnostic{true, {"Found begin event"}};
        return MatchDiagnostic{false};
      },
      [&](EndFenceEvent event) {
        if (event.descriptor().find("Kokkos::resize(View)") !=
            std::string::npos)
          return MatchDiagnostic{true, {"Found end event"}};
        return MatchDiagnostic{false};
      },
      [&](BeginParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found begin event"}};
        return MatchDiagnostic{false};
      },
      [&](EndParallelForEvent event) {
        if (event.descriptor().find("initialization") != std::string::npos)
          return MatchDiagnostic{true, {"Found end event"}};
        return MatchDiagnostic{false};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, realloc_exec_space_scatterview) {
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<typename TEST_EXECSPACE::memory_space,
                   Kokkos::CudaUVMSpace>::value)
    GTEST_SKIP() << "skipping since CudaUVMSpace requires additional fences";
#endif
// FIXME_THREADS The Threads backend fences every parallel_for
#ifdef KOKKOS_ENABLE_THREADS
  if (std::is_same<typename TEST_EXECSPACE, Kokkos::Threads>::value)
    GTEST_SKIP() << "skipping since the Threads backend isn't asynchronous";
#endif

  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableFences());
  using view_type = Kokkos::Experimental::ScatterView<
      int*, typename TEST_EXECSPACE::array_layout, TEST_EXECSPACE>;
  view_type outer_view, outer_view2;

  auto success = validate_absence(
      [&]() {
        view_type inner_view(Kokkos::view_alloc(TEST_EXECSPACE{}, "bla"), 8);
        // Avoid testing the destructor
        outer_view = inner_view;
        Kokkos::realloc(
            Kokkos::view_alloc(Kokkos::WithoutInitializing, TEST_EXECSPACE{}),
            inner_view, 10);
        EXPECT_EQ(inner_view.subview().label(), "bla");
        outer_view2 = inner_view;
        Kokkos::realloc(Kokkos::view_alloc(TEST_EXECSPACE{}), inner_view, 10);
        EXPECT_EQ(inner_view.subview().label(), "bla");
      },
      [&](BeginFenceEvent event) {
        if ((event.descriptor().find("Debug Only Check for Execution Error") !=
             std::string::npos) ||
            (event.descriptor().find("HostSpace fence") != std::string::npos))
          return MatchDiagnostic{false};
        return MatchDiagnostic{true, {"Found fence event!"}};
      });
  ASSERT_TRUE(success);
  listen_tool_events(Config::DisableAll());
}

TEST(TEST_CATEGORY, create_mirror_no_init_dynrankview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::DynRankView<int, TEST_EXECSPACE> device_view("device view", 10);
  Kokkos::DynRankView<int, Kokkos::HostSpace> host_view("host view", 10);

  auto success = validate_absence(
      [&]() {
        auto mirror_device =
            Kokkos::create_mirror(Kokkos::WithoutInitializing, device_view);
        auto mirror_host = Kokkos::create_mirror(Kokkos::WithoutInitializing,
                                                 TEST_EXECSPACE{}, host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, TEST_EXECSPACE{}, host_view);
      },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found begin event"}};
      },
      [&](EndParallelForEvent) {
        return MatchDiagnostic{true, {"Found end event"}};
      });
  ASSERT_TRUE(success);
}

TEST(TEST_CATEGORY, create_mirror_no_init_dynrankview_viewctor) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::DynRankView<int, Kokkos::DefaultExecutionSpace> device_view(
      "device view", 10);
  Kokkos::DynRankView<int, Kokkos::HostSpace> host_view("host view", 10);

  auto success = validate_absence(
      [&]() {
        auto mirror_device = Kokkos::create_mirror(
            Kokkos::view_alloc(Kokkos::WithoutInitializing), device_view);
        auto mirror_host = Kokkos::create_mirror(
            Kokkos::view_alloc(Kokkos::WithoutInitializing,
                               Kokkos::DefaultExecutionSpace{}),
            host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::WithoutInitializing), device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::WithoutInitializing,
                               Kokkos::DefaultExecutionSpace{}),
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

TEST(TEST_CATEGORY, create_mirror_view_and_copy_dynrankview) {
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<typename TEST_EXECSPACE::memory_space,
                   Kokkos::CudaUVMSpace>::value)
    return;
#endif
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels(),
                     Config::EnableFences());

  Kokkos::DynRankView<int, Kokkos::HostSpace> host_view("host view", 10);
  decltype(Kokkos::create_mirror_view_and_copy(TEST_EXECSPACE{},
                                               host_view)) device_view;

  auto success = validate_absence(
      [&]() {
        auto mirror_device = Kokkos::create_mirror_view_and_copy(
            Kokkos::view_alloc(TEST_EXECSPACE{},
                               typename TEST_EXECSPACE::memory_space{}),
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

TEST(TEST_CATEGORY, create_mirror_no_init_offsetview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::Experimental::OffsetView<int*, TEST_EXECSPACE> device_view(
      "device view", {0, 10});
  Kokkos::Experimental::OffsetView<int*, Kokkos::HostSpace> host_view(
      "host view", {0, 10});

  auto success = validate_absence(
      [&]() {
        auto mirror_device =
            Kokkos::create_mirror(Kokkos::WithoutInitializing, device_view);
        auto mirror_host = Kokkos::create_mirror(Kokkos::WithoutInitializing,
                                                 TEST_EXECSPACE{}, host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, TEST_EXECSPACE{}, host_view);
      },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found begin event"}};
      },
      [&](EndParallelForEvent) {
        return MatchDiagnostic{true, {"Found end event"}};
      });
  ASSERT_TRUE(success);
}

TEST(TEST_CATEGORY, create_mirror_no_init_offsetview_view_ctor) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::Experimental::OffsetView<int*, Kokkos::DefaultExecutionSpace>
      device_view("device view", {0, 10});
  Kokkos::Experimental::OffsetView<int*, Kokkos::HostSpace> host_view(
      "host view", {0, 10});

  auto success = validate_absence(
      [&]() {
        auto mirror_device = Kokkos::create_mirror(
            Kokkos::view_alloc(Kokkos::WithoutInitializing), device_view);
        auto mirror_host = Kokkos::create_mirror(
            Kokkos::view_alloc(Kokkos::WithoutInitializing,
                               Kokkos::DefaultExecutionSpace{}),
            host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::WithoutInitializing), device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::WithoutInitializing,
                               Kokkos::DefaultExecutionSpace{}),
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

TEST(TEST_CATEGORY, create_mirror_view_and_copy_offsetview) {
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<typename TEST_EXECSPACE::memory_space,
                   Kokkos::CudaUVMSpace>::value)
    return;
#endif
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels(),
                     Config::EnableFences());

  Kokkos::Experimental::OffsetView<int*, Kokkos::HostSpace> host_view(
      "host view", {0, 10});
  decltype(Kokkos::create_mirror_view_and_copy(TEST_EXECSPACE{},
                                               host_view)) device_view;

  auto success = validate_absence(
      [&]() {
        auto mirror_device = Kokkos::create_mirror_view_and_copy(
            Kokkos::view_alloc(TEST_EXECSPACE{},
                               typename TEST_EXECSPACE::memory_space{}),
            host_view);
        // Avoid fences for deallocation when mirror_device goes out of scope.
        device_view               = mirror_device;
        auto mirror_device_mirror = Kokkos::create_mirror_view_and_copy(
            Kokkos::view_alloc(TEST_EXECSPACE{},
                               typename TEST_EXECSPACE::memory_space{}),
            mirror_device);
      },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found parallel_for event"}};
      },
      [&](BeginFenceEvent) {
        return MatchDiagnostic{true, {"Found fence event"}};
      });
  ASSERT_TRUE(success);
}

// FIXME OPENMPTARGET
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, create_mirror_no_init_dynamicview) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::Experimental::DynamicView<int*, TEST_EXECSPACE> device_view(
      "device view", 2, 10);
  Kokkos::Experimental::DynamicView<int*, Kokkos::HostSpace> host_view(
      "host view", 2, 10);

  auto success = validate_absence(
      [&]() {
        auto mirror_device =
            Kokkos::create_mirror(Kokkos::WithoutInitializing, device_view);
        auto mirror_host = Kokkos::create_mirror(Kokkos::WithoutInitializing,
                                                 TEST_EXECSPACE{}, host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::WithoutInitializing, TEST_EXECSPACE{}, host_view);
      },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found begin event"}};
      },
      [&](EndParallelForEvent) {
        return MatchDiagnostic{true, {"Found end event"}};
      });
  ASSERT_TRUE(success);
}

TEST(TEST_CATEGORY, create_mirror_view_and_copy_dynamicview) {
#ifdef KOKKOS_ENABLE_CUDA
  if (std::is_same<typename TEST_EXECSPACE::memory_space,
                   Kokkos::CudaUVMSpace>::value)
    return;
#endif
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels(),
                     Config::EnableFences());

  Kokkos::Experimental::DynamicView<int*, Kokkos::HostSpace> host_view(
      "host view", 2, 10);
  decltype(Kokkos::create_mirror_view_and_copy(TEST_EXECSPACE{},
                                               host_view)) device_view;

  auto success = validate_absence(
      [&]() {
        auto mirror_device = Kokkos::create_mirror_view_and_copy(
            Kokkos::view_alloc(TEST_EXECSPACE{},
                               typename TEST_EXECSPACE::memory_space{}),
            host_view);
        // Avoid fences for deallocation when mirror_device goes out of scope.
        device_view               = mirror_device;
        auto mirror_device_mirror = Kokkos::create_mirror_view_and_copy(
            Kokkos::view_alloc(TEST_EXECSPACE{},
                               typename TEST_EXECSPACE::memory_space{}),
            mirror_device);
      },
      [&](BeginFenceEvent event) {
        if (event.descriptor().find("DynamicView::resize_serial: Fence after "
                                    "copying chunks to the device") !=
            std::string::npos)
          return MatchDiagnostic{false};
        return MatchDiagnostic{true, {"Found fence event"}};
      },
      [&](EndFenceEvent) { return MatchDiagnostic{false}; },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found parallel_for event"}};
      });
  ASSERT_TRUE(success);
}
#endif

// FIXME OPENMPTARGET
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(TEST_CATEGORY, create_mirror_no_init_dynamicview_view_ctor) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::Experimental::DynamicView<int*, Kokkos::DefaultExecutionSpace>
      device_view("device view", 2, 10);
  Kokkos::Experimental::DynamicView<int*, Kokkos::HostSpace> host_view(
      "host view", 2, 10);

  auto success = validate_absence(
      [&]() {
        auto mirror_device = Kokkos::create_mirror(
            Kokkos::view_alloc(Kokkos::WithoutInitializing), device_view);
        auto mirror_host = Kokkos::create_mirror(
            Kokkos::view_alloc(Kokkos::WithoutInitializing,
                               Kokkos::DefaultExecutionSpace{}),
            host_view);
        auto mirror_device_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::WithoutInitializing), device_view);
        auto mirror_host_view = Kokkos::create_mirror_view(
            Kokkos::view_alloc(Kokkos::WithoutInitializing,
                               Kokkos::DefaultExecutionSpace{}),
            host_view);
      },
      [&](BeginFenceEvent event) {
        if (event.descriptor().find("DynamicView::resize_serial: Fence after "
                                    "copying chunks to the device") !=
            std::string::npos)
          return MatchDiagnostic{false};
        return MatchDiagnostic{true, {"Found fence event"}};
      },
      [&](EndFenceEvent) { return MatchDiagnostic{false}; },
      [&](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found begin event"}};
      },
      [&](EndParallelForEvent) {
        return MatchDiagnostic{true, {"Found end event"}};
      });
  ASSERT_TRUE(success);
}
#endif
