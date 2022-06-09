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

#include <Kokkos_Core.hpp>

#include "tools/include/ToolTestingUtilities.hpp"

TEST(TEST_CATEGORY, resize_realloc_no_init) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  Kokkos::View<int*** * [1][2][3][4], TEST_EXECSPACE> bla("bla", 5, 6, 7, 8);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(Kokkos::WithoutInitializing, bla, 5, 6, 7, 9);
        Kokkos::realloc(Kokkos::WithoutInitializing, bla, 8, 8, 8, 8);
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

TEST(TEST_CATEGORY, resize_realloc_no_alloc) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels(),
                     Config::EnableAllocs());
  Kokkos::View<int*** * [1][2][3][4], TEST_EXECSPACE> bla("bla", 8, 7, 6, 5);

  auto success = validate_absence(
      [&]() {
        Kokkos::resize(bla, 8, 7, 6, 5);
        Kokkos::realloc(Kokkos::WithoutInitializing, bla, 8, 7, 6, 5);
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
