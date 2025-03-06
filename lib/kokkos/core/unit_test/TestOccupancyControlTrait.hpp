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

template <class... Properties>
void test_policy_execution(const Kokkos::RangePolicy<Properties...>& policy) {
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int){});
}
template <class... Properties>
void test_policy_execution(const Kokkos::TeamPolicy<Properties...>& policy) {
  Kokkos::parallel_for(
      policy,
      KOKKOS_LAMBDA(
          const typename Kokkos::TeamPolicy<Properties...>::member_type&){});
}
template <class... Properties>
void test_policy_execution(const Kokkos::MDRangePolicy<Properties...>& policy) {
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(int, int){});
}

template <class Policy>
void test_prefer_desired_occupancy(Policy policy) {
  using Kokkos::Experimental::DesiredOccupancy;
  using Kokkos::Experimental::MaximizeOccupancy;
  using Kokkos::Experimental::prefer;
  using Kokkos::Experimental::WorkItemProperty;

  // MaximizeOccupancy -> MaximizeOccupancy
  auto const policy_still_no_occ = prefer(policy, MaximizeOccupancy{});
  test_policy_execution(policy_still_no_occ);

  // MaximizeOccupancy -> DesiredOccupancy
  auto const policy_with_occ =
      prefer(policy_still_no_occ, DesiredOccupancy{33});
  test_policy_execution(policy_with_occ);

  // DesiredOccupancy -> DesiredOccupancy
  auto const policy_change_occ = prefer(policy_with_occ, DesiredOccupancy{24});
  test_policy_execution(policy_change_occ);

  // DesiredOccupancy -> DesiredOccupancy w/ hint
  auto policy_with_occ_and_hint = Kokkos::Experimental::require(
      policy_change_occ,
      Kokkos::Experimental::WorkItemProperty::HintLightWeight);
  test_policy_execution(policy_with_occ_and_hint);

  // DesiredOccupancy -> MaximizeOccupancy
  auto const policy_drop_occ =
      prefer(policy_with_occ_and_hint, MaximizeOccupancy{});
  test_policy_execution(policy_drop_occ);
}

// FIXME_MSVC_WITH_CUDA
// This test doesn't compile with CUDA on Windows
#if !(defined(_WIN32) && defined(KOKKOS_ENABLE_CUDA))
TEST(TEST_CATEGORY, occupancy_control) {
  test_prefer_desired_occupancy(Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1));
  test_prefer_desired_occupancy(
      Kokkos::TeamPolicy<TEST_EXECSPACE>{1, Kokkos::AUTO});
  test_prefer_desired_occupancy(
      Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>{{0, 0}, {1, 1}});
}
#endif
}  // namespace
