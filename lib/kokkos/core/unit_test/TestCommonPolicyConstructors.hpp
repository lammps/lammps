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

namespace {

// Dummy policy for testing base class.
template <class... Args>
struct DummyPolicy : Kokkos::Impl::PolicyTraits<Args...> {
  using execution_policy = DummyPolicy;
};

// Asserts that a policy constructor is semiregular.
// Semiregular is copyable and default initializable
// (regular requires equality comparable).
template <class Policy>
constexpr bool check_semiregular() {
  static_assert(std::is_default_constructible_v<Policy>);
  static_assert(std::is_copy_constructible_v<Policy>);
  static_assert(std::is_move_constructible_v<Policy>);
  static_assert(std::is_copy_assignable_v<Policy>);
  static_assert(std::is_move_assignable_v<Policy>);
  static_assert(std::is_destructible_v<Policy>);

  return true;
}

static_assert(check_semiregular<DummyPolicy<>>());
static_assert(check_semiregular<Kokkos::RangePolicy<>>());
static_assert(check_semiregular<Kokkos::TeamPolicy<>>());
static_assert(check_semiregular<Kokkos::MDRangePolicy<Kokkos::Rank<2>>>());

// Assert that occupancy conversion and hints work properly.
template <class Policy>
void test_prefer_desired_occupancy() {
  Policy policy;

  using Kokkos::Experimental::DesiredOccupancy;
  using Kokkos::Experimental::MaximizeOccupancy;
  using Kokkos::Experimental::prefer;
  using Kokkos::Experimental::WorkItemProperty;

  static_assert(!Policy::experimental_contains_desired_occupancy);

  // MaximizeOccupancy -> MaximizeOccupancy
  auto const policy_still_no_occ = prefer(policy, MaximizeOccupancy{});
  static_assert(
      !decltype(policy_still_no_occ)::experimental_contains_desired_occupancy);

  // MaximizeOccupancy -> DesiredOccupancy
  auto const policy_with_occ =
      prefer(policy_still_no_occ, DesiredOccupancy{33});
  static_assert(
      decltype(policy_with_occ)::experimental_contains_desired_occupancy);
  EXPECT_EQ(policy_with_occ.impl_get_desired_occupancy().value(), 33);

  // DesiredOccupancy -> DesiredOccupancy
  auto const policy_change_occ = prefer(policy_with_occ, DesiredOccupancy{24});
  static_assert(
      decltype(policy_change_occ)::experimental_contains_desired_occupancy);
  EXPECT_EQ(policy_change_occ.impl_get_desired_occupancy().value(), 24);

  // DesiredOccupancy -> DesiredOccupancy w/ hint
  auto policy_with_occ_and_hint = Kokkos::Experimental::require(
      policy_change_occ,
      Kokkos::Experimental::WorkItemProperty::HintLightWeight);
  EXPECT_EQ(policy_with_occ_and_hint.impl_get_desired_occupancy().value(), 24);

  // DesiredOccupancy -> MaximizeOccupancy
  auto const policy_drop_occ =
      prefer(policy_with_occ_and_hint, MaximizeOccupancy{});
  static_assert(
      !decltype(policy_drop_occ)::experimental_contains_desired_occupancy);
}

TEST(TEST_CATEGORY, execution_policy_occupancy_and_hint) {
  test_prefer_desired_occupancy<DummyPolicy<>>();
  test_prefer_desired_occupancy<Kokkos::RangePolicy<>>();
  test_prefer_desired_occupancy<Kokkos::TeamPolicy<>>();
  test_prefer_desired_occupancy<Kokkos::MDRangePolicy<Kokkos::Rank<2>>>();
}

// Check that the policy size does not increase if the user does not specify the
// occupancy (only pay for what you use).
// Disabling since EBO was not working with VS 16.11.3 and CUDA 11.4.2
#if !(defined(_WIN32) && defined(KOKKOS_ENABLE_CUDA))
constexpr bool test_empty_base_optimization() {
  DummyPolicy<> policy;
  static_assert(sizeof(decltype(policy)) == 1);
  using Kokkos::Experimental::DesiredOccupancy;
  using Kokkos::Experimental::prefer;
  static_assert(sizeof(decltype(prefer(policy, DesiredOccupancy{33}))) ==
                sizeof(DesiredOccupancy));
  return true;
}
static_assert(test_empty_base_optimization());
#endif

}  // namespace
