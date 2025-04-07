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

#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

#include <cstdlib>

#include "KokkosExecutionEnvironmentNeverInitializedFixture.hpp"

namespace {

using InitializeFinalize_DeathTest = KokkosExecutionEnvironmentNeverInitialized;

TEST_F(InitializeFinalize_DeathTest, initialize) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_EXIT(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");

  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::initialize();
      },
      "Error: Kokkos::initialize\\(\\) has already been called. Kokkos can be "
      "initialized at most once\\.");

  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        Kokkos::initialize();
      },
      "Error: Kokkos::initialize\\(\\) has already been called. Kokkos can be "
      "initialized at most once\\.");
}

TEST_F(InitializeFinalize_DeathTest, finalize) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_DEATH(
      { Kokkos::finalize(); },
      "Error: Kokkos::finalize\\(\\) may only be called after Kokkos has "
      "been initialized\\.");

  EXPECT_EXIT(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        std::exit(EXIT_SUCCESS);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");

  EXPECT_DEATH(
      {
        Kokkos::initialize();
        Kokkos::finalize();
        Kokkos::finalize();
      },
      "Error: Kokkos::finalize\\(\\) has already been called\\.");
}

TEST_F(InitializeFinalize_DeathTest, is_initialized) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_EXIT(
      {
        bool success = true;
        success &= !Kokkos::is_initialized();
        Kokkos::initialize();
        success &= Kokkos::is_initialized();
        Kokkos::finalize();
        success &= !Kokkos::is_initialized();
        std::exit(success ? EXIT_SUCCESS : EXIT_FAILURE);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
}

TEST_F(InitializeFinalize_DeathTest, is_finalized) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  EXPECT_EXIT(
      {
        bool success = true;
        success &= !Kokkos::is_finalized();
        Kokkos::initialize();
        success &= !Kokkos::is_finalized();
        Kokkos::finalize();
        success &= Kokkos::is_finalized();
        std::exit(success ? EXIT_SUCCESS : EXIT_FAILURE);
      },
      ::testing::ExitedWithCode(EXIT_SUCCESS), "");
}

}  // namespace
