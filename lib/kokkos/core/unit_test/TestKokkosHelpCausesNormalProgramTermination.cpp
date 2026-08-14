// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include "KokkosExecutionEnvironmentNeverInitializedFixture.hpp"

namespace {

using KokkosHelpCausesNormalProgramTermination_DeathTest =
    KokkosExecutionEnvironmentNeverInitialized;

TEST_F(KokkosHelpCausesNormalProgramTermination_DeathTest,
       print_help_and_exit_early) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  int argc = 1;

  char const *argv[] = {
      "--kokkos-help",
      nullptr,
  };

  ::testing::internal::CaptureStdout();

  EXPECT_EXIT(
      {
        Kokkos::initialize(argc, const_cast<char **>(argv));
        Kokkos::abort("better exit before getting there");
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");

  (void)::testing::internal::GetCapturedStdout();
}

}  // namespace
