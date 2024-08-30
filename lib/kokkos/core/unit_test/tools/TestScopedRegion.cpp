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

#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <gtest/gtest.h>

#include <string>
#include <stack>

namespace {

std::stack<std::string> test_region_stack;

// NOTE: cannot use lambdas because they can only be converted to function
// pointers if they don't capture anything
void test_push_region(char const *label) { test_region_stack.push(label); }

void test_pop_region() { test_region_stack.pop(); }

TEST(defaultdevicetype, scoped_profile_region) {
  Kokkos::Tools::Experimental::set_push_region_callback(test_push_region);
  Kokkos::Tools::Experimental::set_pop_region_callback(test_pop_region);

  ASSERT_TRUE(test_region_stack.empty());

  {
    std::string outer_identifier = "outer";
    Kokkos::Profiling::ScopedRegion guard_outer(outer_identifier);

    ASSERT_EQ(test_region_stack.size(), 1u);
    ASSERT_EQ(test_region_stack.top(), outer_identifier);

    {
      std::string inner_identifier = "inner";
      Kokkos::Profiling::ScopedRegion guard_inner(inner_identifier);
      ASSERT_EQ(test_region_stack.size(), 2u);
      ASSERT_EQ(test_region_stack.top(), inner_identifier);
    }

    ASSERT_EQ(test_region_stack.size(), 1u);
    ASSERT_EQ(test_region_stack.top(), outer_identifier);
  }

  ASSERT_TRUE(test_region_stack.empty());

  // Unset callbacks
  Kokkos::Tools::Experimental::set_push_region_callback(nullptr);
  Kokkos::Tools::Experimental::set_pop_region_callback(nullptr);
}

using Kokkos::Profiling::ScopedRegion;
static_assert(!std::is_default_constructible<ScopedRegion>::value);
static_assert(!std::is_copy_constructible<ScopedRegion>::value);
static_assert(!std::is_move_constructible<ScopedRegion>::value);
static_assert(!std::is_copy_assignable<ScopedRegion>::value);
static_assert(!std::is_move_assignable<ScopedRegion>::value);

}  // namespace
