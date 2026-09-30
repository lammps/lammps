// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Timer.hpp>

#include <thread>
#include <type_traits>
#include <utility>

namespace {

TEST(TEST_CATEGORY, timer) {
  using namespace std::chrono_literals;

  Kokkos::Timer t;
  std::this_thread::sleep_for(5ms);
  auto elapsed = t.seconds();
  EXPECT_GE(elapsed, .005);
  EXPECT_LT(elapsed, 1.);

  std::this_thread::sleep_for(10ms);
  auto elapsed2 = std::as_const(t).seconds();
  EXPECT_GE(elapsed2, .015);
  EXPECT_GT(elapsed2, elapsed);

  t.reset();
  std::this_thread::sleep_for(5ms);
  auto elapsed3 = t.seconds();
  EXPECT_GE(elapsed3, .005);
  // Using the line below turned out to be problematic since there is no
  // guaranteed upper bound for the time std::this_thread::sleep_for is
  // blocking for.
  // EXPECT_LT(elapsed3, elapsed2);
}

static_assert(!std::is_copy_constructible_v<Kokkos::Timer>);
static_assert(!std::is_move_constructible_v<Kokkos::Timer>);
static_assert(!std::is_copy_assignable_v<Kokkos::Timer>);
static_assert(!std::is_move_assignable_v<Kokkos::Timer>);

}  // namespace
