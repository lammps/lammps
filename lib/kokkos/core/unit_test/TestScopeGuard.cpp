// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <cstdlib>
#include <gtest/gtest.h>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include "KokkosExecutionEnvironmentNeverInitializedFixture.hpp"

namespace {

using ScopeGuard_DeathTest = KokkosExecutionEnvironmentNeverInitialized;

/**
 * Test to create a scope guard normally.
 */
TEST_F(ScopeGuard_DeathTest, create) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  // run it in a different process so side effects are not kept
  EXPECT_EXIT(
      {
        {
          Kokkos::ScopeGuard guard{};

          if (!Kokkos::is_initialized()) std::exit(EXIT_FAILURE);
          if (Kokkos::is_finalized()) std::exit(EXIT_FAILURE);
        }

        if (Kokkos::is_initialized()) std::exit(EXIT_FAILURE);
        if (!Kokkos::is_finalized()) std::exit(EXIT_FAILURE);

        std::exit(EXIT_SUCCESS);
      },
      testing::ExitedWithCode(EXIT_SUCCESS), "");
}

/**
 * Test to create a scope guard with an argument.
 */
TEST_F(ScopeGuard_DeathTest, create_argument) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  // run it in a different process so side effects are not kept
  EXPECT_EXIT(
      {
        {
          Kokkos::InitializationSettings settings{};
          Kokkos::ScopeGuard guard{settings};
        }

        std::exit(EXIT_SUCCESS);
      },
      testing::ExitedWithCode(EXIT_SUCCESS), "");
}

/**
 * Test to create another scope guard when one has been created.
 */
TEST_F(ScopeGuard_DeathTest, create_while_initialize) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(
      {
        Kokkos::ScopeGuard guard1{};

        // create a second scope guard while there is one already existing
        Kokkos::ScopeGuard guard2{};
      },
      "Creating a ScopeGuard while Kokkos is initialized");
}

/**
 * Test to create a scope guard when initialization has been done manually.
 */
TEST_F(ScopeGuard_DeathTest, create_after_initialize) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(
      {
        Kokkos::initialize();

        // create a scope guard after manual initialization
        Kokkos::ScopeGuard guard{};
      },
      "Creating a ScopeGuard while Kokkos is initialized");
}

/**
 * Test to create another scope guard when one has been destroyed.
 */
TEST_F(ScopeGuard_DeathTest, create_after_finalize) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(
      {
        { Kokkos::ScopeGuard guard1{}; }

        // create a second scope guard while the first one has been destroyed
        // already
        Kokkos::ScopeGuard guard2{};
      },
      "Creating a ScopeGuard after Kokkos was finalized");
}

/**
 * Test to destroy a scope guard when finalization has been done manually.
 */
TEST_F(ScopeGuard_DeathTest, destroy_after_finalize) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH(
      {
        // create a scope guard and finalize it manually
        Kokkos::ScopeGuard guard{};
        Kokkos::finalize();
      },
      "Destroying a ScopeGuard after Kokkos was finalized");
}

/**
 * Static tests
 */

// Test scope guard is not copyable.
static_assert(!std::is_copy_assignable_v<Kokkos::ScopeGuard>);
static_assert(!std::is_copy_constructible_v<Kokkos::ScopeGuard>);

// Test scope guard is not movable.
static_assert(!std::is_move_assignable_v<Kokkos::ScopeGuard>);
static_assert(!std::is_move_constructible_v<Kokkos::ScopeGuard>);

}  // namespace
