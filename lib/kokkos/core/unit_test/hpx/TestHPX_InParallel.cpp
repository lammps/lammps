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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

// These tests specifically check that work dispatched to independent instances
// is synchronized correctly on fences. A previous bug that this protects
// against is work being mistakenly dispatched to the default instance, but the
// fence fencing the independent instance. In that case these tests will fail.

namespace {
inline constexpr int n = 1 << 10;

TEST(hpx, in_parallel_for_range_policy) {
  Kokkos::View<bool *, Kokkos::Experimental::HPX> a("a", n);

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  Kokkos::RangePolicy<Kokkos::Experimental::HPX> policy(0, n);
  Kokkos::parallel_for(
      "parallel_for_range_policy", policy, KOKKOS_LAMBDA(const int i) {
        a(i) = Kokkos::Experimental::HPX::in_parallel();
      });

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());
  Kokkos::fence();
  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  for (int i = 0; i < n; ++i) {
    ASSERT_TRUE(a(i));
  }
}

TEST(hpx, in_parallel_for_mdrange_policy) {
  Kokkos::View<bool *, Kokkos::Experimental::HPX> a("a", n);

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  Kokkos::MDRangePolicy<Kokkos::Experimental::HPX, Kokkos::Rank<2>> policy(
      {0, 0}, {n, 1});
  Kokkos::parallel_for(
      "parallel_for_mdrange_policy", policy,
      KOKKOS_LAMBDA(const int i, const int) {
        a(i) = Kokkos::Experimental::HPX::in_parallel();
      });

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());
  Kokkos::fence();
  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  for (int i = 0; i < n; ++i) {
    ASSERT_TRUE(a(i));
  }
}

TEST(hpx, in_parallel_for_team_policy) {
  Kokkos::View<bool *, Kokkos::Experimental::HPX> a("a", n);

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  Kokkos::TeamPolicy<Kokkos::Experimental::HPX> policy(n, 1);
  using member_type = decltype(policy)::member_type;
  Kokkos::parallel_for(
      "parallel_for_team_policy", policy,
      KOKKOS_LAMBDA(const member_type &handle) {
        a(handle.league_rank()) = Kokkos::Experimental::HPX::in_parallel();
      });

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());
  Kokkos::fence();
  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  for (int i = 0; i < n; ++i) {
    ASSERT_TRUE(a(i));
  }
}

TEST(hpx, in_parallel_reduce_range_policy) {
  Kokkos::View<int, Kokkos::Experimental::HPX> a("a");

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  Kokkos::RangePolicy<Kokkos::Experimental::HPX> policy(0, n);
  Kokkos::parallel_reduce(
      "parallel_reduce_range_policy", policy,
      KOKKOS_LAMBDA(const int, int &x) {
        if (!Kokkos::Experimental::HPX::in_parallel()) {
          ++x;
        }
      },
      a);

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());
  Kokkos::fence();
  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  ASSERT_EQ(a(), 0);
}

TEST(hpx, in_parallel_reduce_mdrange_policy) {
  Kokkos::View<int, Kokkos::Experimental::HPX> a("a");

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  Kokkos::MDRangePolicy<Kokkos::Experimental::HPX, Kokkos::Rank<2>> policy(
      {0, 0}, {n, 1});
  Kokkos::parallel_reduce(
      "parallel_reduce_mdrange_policy", policy,
      KOKKOS_LAMBDA(const int, const int, int &x) {
        if (!Kokkos::Experimental::HPX::in_parallel()) {
          ++x;
        }
      },
      a);

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());
  Kokkos::fence();
  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  ASSERT_EQ(a(), 0);
}

TEST(hpx, in_parallel_reduce_team_policy) {
  Kokkos::View<int, Kokkos::Experimental::HPX> a("a");

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  Kokkos::TeamPolicy<Kokkos::Experimental::HPX> policy(n, 1);
  using member_type = decltype(policy)::member_type;
  Kokkos::parallel_reduce(
      "parallel_reduce_team_policy", policy,
      KOKKOS_LAMBDA(const member_type &, int &x) {
        if (!Kokkos::Experimental::HPX::in_parallel()) {
          ++x;
        }
      },
      a);

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());
  Kokkos::fence();
  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  ASSERT_EQ(a(), 0);
}

TEST(hpx, in_parallel_scan_range_policy) {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a("a", n);

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  Kokkos::RangePolicy<Kokkos::Experimental::HPX> policy(0, n);
  Kokkos::parallel_scan(
      "parallel_scan_range_policy", policy,
      KOKKOS_LAMBDA(const int, int &x, bool) {
        if (!Kokkos::Experimental::HPX::in_parallel()) {
          ++x;
        }
      },
      a);

  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());
  Kokkos::fence();
  ASSERT_FALSE(Kokkos::Experimental::HPX::in_parallel());

  for (int i = 0; i < n; ++i) {
    ASSERT_EQ(a(i), 0);
  }
}
}  // namespace
