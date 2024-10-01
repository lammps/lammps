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

#include <cstdlib>
#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace {

/**
 * Fixture that checks Kokkos is neither initialized nor finalized before and
 * after the test.
 */
class AssertEnvironmentTest : public ::testing::Test {
 protected:
  void SetUp() override {
    ASSERT_FALSE(Kokkos::is_initialized());
    ASSERT_FALSE(Kokkos::is_finalized());
  }

  void TearDown() override {
    ASSERT_FALSE(Kokkos::is_initialized());
    ASSERT_FALSE(Kokkos::is_finalized());
  }
};

using scope_guard_DeathTest = AssertEnvironmentTest;

/**
 * Test to create a scope guard normally.
 */
TEST_F(scope_guard_DeathTest, create) {
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
TEST_F(scope_guard_DeathTest, create_argument) {
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
TEST_F(scope_guard_DeathTest, create_while_initialize) {
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
TEST_F(scope_guard_DeathTest, create_after_initialize) {
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
TEST_F(scope_guard_DeathTest, create_after_finalize) {
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
TEST_F(scope_guard_DeathTest, destroy_after_finalize) {
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
static_assert(!std::is_copy_assignable<Kokkos::ScopeGuard>());
static_assert(!std::is_copy_constructible<Kokkos::ScopeGuard>());

// Test scope guard is not movable.
static_assert(!std::is_move_assignable<Kokkos::ScopeGuard>());
static_assert(!std::is_move_constructible<Kokkos::ScopeGuard>());

}  // namespace
