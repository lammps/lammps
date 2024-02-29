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

#include <regex>
#include <Kokkos_Core.hpp>

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
#elif defined(KOKKOS_ENABLE_OPENACC)  // FIXME_OPENACC
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::OpenACC>::value) {
    TestAbortPrintingToStdout<ExecutionSpace>();
  } else {
    TestAbortCausingAbnormalProgramTerminationAndPrinting<ExecutionSpace>();
  }
#elif defined(KOKKOS_ENABLE_SYCL)     // FIXME_SYCL
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value) {
#ifdef NDEBUG
    TestAbortPrintingToStdout<ExecutionSpace>();
#else
    TestAbortCausingAbnormalProgramTerminationAndPrinting<ExecutionSpace>();
#endif
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
