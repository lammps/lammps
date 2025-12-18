// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <regex>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

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
  if (std::is_same_v<ExecutionSpace, Kokkos::SYCL>) {
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
// FIXME_OPENACC FIXME_NVHPC: NVHPC fails when targetting CPUs.
#if defined(KOKKOS_ENABLE_OPENACC) && defined(KOKKOS_COMPILER_NVHPC) && \
    defined(KOKKOS_ENABLE_OPENACC_FORCE_HOST_AS_DEVICE)
  if (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenACC>) {
    GTEST_SKIP()
        << "skipping since the OpenACC backend compiled by NVHPC for CPU "
           "crashes at runtime.";
  }
#endif
  test_abort_from_device<TEST_EXECSPACE>();
}
