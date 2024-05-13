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
#include <TestHPX_Category.hpp>

// These tests specifically check that work dispatched to independent instances
// is synchronized correctly on fences. A previous bug that this protects
// against is work being mistakenly dispatched to the default instance, but the
// fence fencing the independent instance. In that case these tests will fail.

namespace {
inline constexpr int n = 1 << 10;

TEST(hpx, independent_instances_synchronization_parallel_for_range_policy) {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a("a", n);

  Kokkos::Experimental::HPX instance{
      Kokkos::Experimental::HPX::instance_mode::independent};
  Kokkos::RangePolicy<Kokkos::Experimental::HPX> policy(instance, 0, n);
  Kokkos::parallel_for(
      "parallel_for_range_policy", policy,
      KOKKOS_LAMBDA(const auto i) { a[i] = i; });

  instance.fence();

  for (int i = 0; i < n; ++i) {
    ASSERT_EQ(a[i], i);
  }
}

TEST(hpx, independent_instances_synchronization_parallel_for_mdrange_policy) {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a("a", n);

  Kokkos::Experimental::HPX instance{
      Kokkos::Experimental::HPX::instance_mode::independent};
  Kokkos::MDRangePolicy<Kokkos::Experimental::HPX, Kokkos::Rank<2>> policy(
      instance, {{0, 0}}, {{n, 1}});
  Kokkos::parallel_for(
      "parallel_for_mdrange_policy", policy,
      KOKKOS_LAMBDA(const auto i, const auto) { a[i] = i; });

  instance.fence();

  for (int i = 0; i < n; ++i) {
    ASSERT_EQ(a[i], i);
  }
}

TEST(hpx, independent_instances_synchronization_parallel_for_team_policy) {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a("a", n);

  Kokkos::Experimental::HPX instance{
      Kokkos::Experimental::HPX::instance_mode::independent};
  Kokkos::TeamPolicy<Kokkos::Experimental::HPX> policy(instance, n, 1);
  Kokkos::parallel_for(
      "parallel_for_team_policy", policy, KOKKOS_LAMBDA(const auto &handle) {
        a[handle.league_rank()] = handle.league_rank();
      });

  instance.fence();

  for (int i = 0; i < n; ++i) {
    ASSERT_EQ(a[i], i);
  }
}

TEST(hpx, independent_instances_synchronization_parallel_reduce_range_policy) {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a("a", n);
  Kokkos::View<int, Kokkos::Experimental::HPX> b("b");

  Kokkos::Experimental::HPX instance{
      Kokkos::Experimental::HPX::instance_mode::independent};
  Kokkos::RangePolicy<Kokkos::Experimental::HPX> policy(instance, 0, n);
  Kokkos::parallel_reduce(
      "parallel_reduce_range_policy", policy,
      KOKKOS_LAMBDA(const int i, int &) { a[i] = i; }, b);

  instance.fence();

  for (int i = 0; i < n; ++i) {
    ASSERT_EQ(a[i], i);
  }
}

TEST(hpx,
     independent_instances_synchronization_parallel_reduce_mdrange_policy) {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a("a", n);
  Kokkos::View<int, Kokkos::Experimental::HPX> b("b");

  Kokkos::Experimental::HPX instance{
      Kokkos::Experimental::HPX::instance_mode::independent};
  Kokkos::MDRangePolicy<Kokkos::Experimental::HPX, Kokkos::Rank<2>> policy(
      instance, {{0, 0}}, {{n, 1}});
  Kokkos::parallel_reduce(
      "parallel_reduce_mdrange_policy", policy,
      KOKKOS_LAMBDA(const int i, const int, int &) { a[i] = i; }, b);

  instance.fence();

  for (int i = 0; i < n; ++i) {
    ASSERT_EQ(a[i], i);
  }
}

TEST(hpx, independent_instances_synchronization_parallel_reduce_team_policy) {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a("a", n);
  Kokkos::View<int, Kokkos::Experimental::HPX> b("b");

  Kokkos::Experimental::HPX instance{
      Kokkos::Experimental::HPX::instance_mode::independent};
  Kokkos::TeamPolicy<Kokkos::Experimental::HPX> policy(instance, n, 1);
  Kokkos::parallel_reduce(
      "parallel_reduce_team_policy", policy,
      KOKKOS_LAMBDA(const decltype(policy)::member_type &handle, int &) {
        a[handle.league_rank()] = handle.league_rank();
      },
      b);

  instance.fence();

  for (int i = 0; i < n; ++i) {
    ASSERT_EQ(a[i], i);
  }
}

TEST(hpx, independent_instances_synchronization_parallel_scan_range_policy) {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a("a", n);
  Kokkos::View<int *, Kokkos::Experimental::HPX> b("b", n);

  Kokkos::Experimental::HPX instance{
      Kokkos::Experimental::HPX::instance_mode::independent};
  Kokkos::RangePolicy<Kokkos::Experimental::HPX> policy(instance, 0, n);
  Kokkos::parallel_scan(
      "parallel_scan_range_policy", policy,
      KOKKOS_LAMBDA(const int i, int &, bool final) {
        if (!final) {
          a[i] = i;
        }
      },
      b);

  instance.fence();

  for (int i = 0; i < n; ++i) {
    ASSERT_EQ(a[i], i);
  }
}
}  // namespace
