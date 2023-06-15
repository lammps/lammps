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

#include <regex>
#include <Kokkos_Core.hpp>

#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC
TEST(TEST_CATEGORY_DEATH, abort_from_host) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  char msg[] = "Goodbye cruel world";
  EXPECT_DEATH({ Kokkos::abort(msg); }, msg);
}

template <class ExecutionSpace>
struct TestAbortPrintingToStdout {
  TestAbortPrintingToStdout() {
    ::testing::internal::CaptureStdout();
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0, 1), *this);
    Kokkos::fence();
    auto const captured = ::testing::internal::GetCapturedStdout();
    EXPECT_TRUE(std::regex_search(captured,
                                  std::regex("move along nothing to see here")))
        << "here is what was printed to stdout \"" << captured << "\"";
  }
  KOKKOS_FUNCTION void operator()(int) const {
    Kokkos::abort("move along nothing to see here");
  }
};

template <class ExecutionSpace>
struct TestAbortCausingAbnormalProgramTerminationButIgnoringErrorMessage {
  TestAbortCausingAbnormalProgramTerminationButIgnoringErrorMessage() {
    EXPECT_DEATH(
        {
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0, 1),
                               *this);
          Kokkos::fence();
        },
        ".*");
  }
  KOKKOS_FUNCTION void operator()(int) const { Kokkos::abort("ignored"); }
};

template <class ExecutionSpace>
struct TestAbortCausingAbnormalProgramTerminationAndPrinting {
  TestAbortCausingAbnormalProgramTerminationAndPrinting() {
    EXPECT_DEATH(
        {
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(0, 1),
                               *this);
          Kokkos::fence();
        },
        "Meurs, pourriture communiste !");
  }
  KOKKOS_FUNCTION void operator()(int) const {
    Kokkos::abort("Meurs, pourriture communiste !");
  }
};

template <class ExecutionSpace>
void test_abort_from_device() {
#if defined(KOKKOS_ENABLE_OPENMPTARGET)  // FIXME_OPENMPTARGET
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::OpenMPTarget>::value) {
    TestAbortPrintingToStdout<ExecutionSpace>();
  } else {
    TestAbortCausingAbnormalProgramTerminationAndPrinting<ExecutionSpace>();
  }
#elif defined(KOKKOS_ENABLE_SYCL)  // FIXME_SYCL
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value) {
#ifdef NDEBUG
    TestAbortPrintingToStdout<ExecutionSpace>();
#else
    TestAbortCausingAbnormalProgramTerminationAndPrinting<ExecutionSpace>();
#endif
  } else {
    TestAbortCausingAbnormalProgramTerminationAndPrinting<ExecutionSpace>();
  }
#elif defined(KOKKOS_IMPL_HIP_ABORT_DOES_NOT_PRINT_MESSAGE)
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::HIP>::value) {
    TestAbortCausingAbnormalProgramTerminationButIgnoringErrorMessage<
        ExecutionSpace>();
  } else {
    TestAbortCausingAbnormalProgramTerminationAndPrinting<ExecutionSpace>();
  }
#else
  TestAbortCausingAbnormalProgramTerminationAndPrinting<ExecutionSpace>();
#endif
}

TEST(TEST_CATEGORY_DEATH, abort_from_device) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  test_abort_from_device<TEST_EXECSPACE>();
}
#endif
