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
#include <impl/Kokkos_Stacktrace.hpp>
#include <cstdio>
#include <cstdint>
#include <sstream>
#include <type_traits>

namespace Test {

template <class ExecutionSpace, class DataType>
struct TestTeamScan {
  using execution_space = ExecutionSpace;
  using value_type      = DataType;
  using policy_type     = Kokkos::TeamPolicy<execution_space>;
  using member_type     = typename policy_type::member_type;
  using view_type       = Kokkos::View<value_type**, execution_space>;

  view_type a_d;
  view_type a_r;
  int32_t M = 0;
  int32_t N = 0;

  KOKKOS_FUNCTION
  void operator()(const member_type& team) const {
    auto leagueRank = team.league_rank();

    auto beg = 0;
    auto end = N;

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, beg, end),
        [&](const int i) { a_d(leagueRank, i) = leagueRank * N + i; });

    Kokkos::parallel_scan(Kokkos::TeamThreadRange(team, beg, end),
                          [&](int i, DataType& val, const bool final) {
                            val += a_d(leagueRank, i);
                            if (final) a_r(leagueRank, i) = val;
                          });
  }

  auto operator()(int32_t _M, int32_t _N) {
    std::stringstream ss;
    ss << Kokkos::Impl::demangle(typeid(*this).name());
    ss << "(/*M=*/" << _M << ", /*N=*/" << _N << ")";
    std::string const test_id = ss.str();

    M   = _M;
    N   = _N;
    a_d = view_type("a_d", M, N);
    a_r = view_type("a_r", M, N);

    // Set team size explicitly to check whether non-power-of-two team sizes can
    // be used.
    if (ExecutionSpace().concurrency() > 10000)
      Kokkos::parallel_for(policy_type(M, 127), *this);
    else if (ExecutionSpace().concurrency() > 2)
      Kokkos::parallel_for(policy_type(M, 3), *this);
    else
      Kokkos::parallel_for(policy_type(M, 1), *this);

    auto a_i = Kokkos::create_mirror_view(a_d);
    auto a_o = Kokkos::create_mirror_view(a_r);
    Kokkos::deep_copy(a_i, a_d);
    Kokkos::deep_copy(a_o, a_r);

    for (int32_t i = 0; i < M; ++i) {
      value_type scan_ref = 0;
      value_type scan_calc;
      value_type abs_err = 0;
      // each fp addition is subject to small loses in precision and these
      // compound as loop so we set the base error to be the machine epsilon and
      // then add in another epsilon each iteration. For example, with CUDA
      // backend + 32-bit float + large N values (e.g. 1,000) + high
      // thread-counts (e.g. 1024), this test will fail w/o epsilon
      // accommodation
      constexpr value_type epsilon = std::numeric_limits<value_type>::epsilon();
      for (int32_t j = 0; j < N; ++j) {
        scan_ref += a_i(i, j);
        scan_calc = a_o(i, j);
        if (std::is_integral<value_type>::value) {
          ASSERT_EQ(scan_ref, scan_calc)
              << test_id
              << " calculated scan output value differs from reference at "
                 "indices i="
              << i << " and j=" << j;
        } else {
          abs_err += epsilon;
          ASSERT_NEAR(scan_ref, scan_calc, abs_err)
              << test_id
              << " calculated scan output value differs from reference at "
                 "indices i="
              << i << " and j=" << j;
        }
      }
    }
  }
};

TEST(TEST_CATEGORY, team_scan) {
  TestTeamScan<TEST_EXECSPACE, int32_t>{}(0, 0);
  TestTeamScan<TEST_EXECSPACE, int32_t>{}(0, 1);
  TestTeamScan<TEST_EXECSPACE, int32_t>{}(1, 0);
  TestTeamScan<TEST_EXECSPACE, uint32_t>{}(99, 32);
  TestTeamScan<TEST_EXECSPACE, uint32_t>{}(139, 64);
  TestTeamScan<TEST_EXECSPACE, uint32_t>{}(163, 128);
  TestTeamScan<TEST_EXECSPACE, int64_t>{}(433, 256);
  TestTeamScan<TEST_EXECSPACE, uint64_t>{}(976, 512);
  TestTeamScan<TEST_EXECSPACE, uint64_t>{}(1234, 1024);
  TestTeamScan<TEST_EXECSPACE, float>{}(2596, 34);
  TestTeamScan<TEST_EXECSPACE, double>{}(2596, 59);
  TestTeamScan<TEST_EXECSPACE, float>{}(2596, 65);
  TestTeamScan<TEST_EXECSPACE, double>{}(2596, 371);
  TestTeamScan<TEST_EXECSPACE, int64_t>{}(2596, 987);
  TestTeamScan<TEST_EXECSPACE, double>{}(2596, 1311);
}

// Temporary: This condition will progressively be reduced when parallel_scan
// with return value will be implemented for more backends.
#if !defined(KOKKOS_ENABLE_OPENACC)
template <class ExecutionSpace, class DataType>
struct TestTeamScanRetVal {
  using execution_space = ExecutionSpace;
  using value_type      = DataType;
  using policy_type     = Kokkos::TeamPolicy<execution_space>;
  using member_type     = typename policy_type::member_type;
  using view_1d_type    = Kokkos::View<value_type*, execution_space>;
  using view_2d_type    = Kokkos::View<value_type**, execution_space>;

  view_2d_type a_d;
  view_2d_type a_r;
  view_1d_type a_s;
  int32_t M = 0;
  int32_t N = 0;

  KOKKOS_FUNCTION
  void operator()(const member_type& team) const {
    auto leagueRank = team.league_rank();

    auto beg = 0;
    auto end = N;

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, beg, end),
        [&](const int i) { a_d(leagueRank, i) = leagueRank * N + i; });

    DataType accum;
    Kokkos::parallel_scan(
        Kokkos::TeamThreadRange(team, beg, end),
        [&](int i, DataType& val, const bool final) {
          val += a_d(leagueRank, i);
          if (final) a_r(leagueRank, i) = val;
        },
        accum);

    // Save return value from parallel_scan
    Kokkos::single(Kokkos::PerTeam(team), [&]() { a_s(leagueRank) = accum; });
  }

  auto operator()(int32_t _M, int32_t _N) {
    std::stringstream ss;
    ss << Kokkos::Impl::demangle(typeid(*this).name());
    ss << "(/*M=*/" << _M << ", /*N=*/" << _N << ")";
    std::string const test_id = ss.str();

    M   = _M;
    N   = _N;
    a_d = view_2d_type("a_d", M, N);
    a_r = view_2d_type("a_r", M, N);
    a_s = view_1d_type("a_s", M);

    // Set team size explicitly to check whether non-power-of-two team sizes can
    // be used.
    if (ExecutionSpace().concurrency() > 10000)
      Kokkos::parallel_for(policy_type(M, 127), *this);
    else if (ExecutionSpace().concurrency() > 2)
      Kokkos::parallel_for(policy_type(M, 3), *this);
    else
      Kokkos::parallel_for(policy_type(M, 1), *this);

    Kokkos::fence();
    auto a_i  = Kokkos::create_mirror_view(a_d);
    auto a_o  = Kokkos::create_mirror_view(a_r);
    auto a_os = Kokkos::create_mirror_view(a_s);
    Kokkos::deep_copy(a_i, a_d);
    Kokkos::deep_copy(a_o, a_r);
    Kokkos::deep_copy(a_os, a_s);

    for (int32_t i = 0; i < M; ++i) {
      value_type scan_ref = 0;
      value_type scan_calc;
      value_type abs_err = 0;
      // each fp addition is subject to small loses in precision and these
      // compound as loop so we set the base error to be the machine epsilon and
      // then add in another epsilon each iteration. For example, with CUDA
      // backend + 32-bit float + large N values (e.g. 1,000) + high
      // thread-counts (e.g. 1024), this test will fail w/o epsilon
      // accommodation
      constexpr value_type epsilon = std::numeric_limits<value_type>::epsilon();
      for (int32_t j = 0; j < N; ++j) {
        scan_ref += a_i(i, j);
        scan_calc = a_o(i, j);
        if (std::is_integral<value_type>::value) {
          ASSERT_EQ(scan_ref, scan_calc)
              << test_id
              << " calculated scan output value differs from reference at "
                 "indices i="
              << i << " and j=" << j;
        } else {
          abs_err += epsilon;
          ASSERT_NEAR(scan_ref, scan_calc, abs_err)
              << test_id
              << " calculated scan output value differs from reference at "
                 "indices i="
              << i << " and j=" << j;
        }
      }
      // Validate return value from parallel_scan
      if (std::is_integral<value_type>::value) {
        ASSERT_EQ(scan_ref, a_os(i));
      } else {
        ASSERT_NEAR(scan_ref, a_os(i), abs_err);
      }
    }
  }
};

TEST(TEST_CATEGORY, team_scan_ret_val) {
  TestTeamScanRetVal<TEST_EXECSPACE, int32_t>{}(0, 0);
  TestTeamScanRetVal<TEST_EXECSPACE, int32_t>{}(0, 1);
  TestTeamScanRetVal<TEST_EXECSPACE, int32_t>{}(1, 0);
  TestTeamScanRetVal<TEST_EXECSPACE, uint32_t>{}(99, 32);
  TestTeamScanRetVal<TEST_EXECSPACE, uint32_t>{}(139, 64);
  TestTeamScanRetVal<TEST_EXECSPACE, uint32_t>{}(163, 128);
  TestTeamScanRetVal<TEST_EXECSPACE, int64_t>{}(433, 256);
  TestTeamScanRetVal<TEST_EXECSPACE, uint64_t>{}(976, 512);
  TestTeamScanRetVal<TEST_EXECSPACE, uint64_t>{}(1234, 1024);
  TestTeamScanRetVal<TEST_EXECSPACE, float>{}(2596, 34);
  TestTeamScanRetVal<TEST_EXECSPACE, double>{}(2596, 59);
  TestTeamScanRetVal<TEST_EXECSPACE, float>{}(2596, 65);
  TestTeamScanRetVal<TEST_EXECSPACE, double>{}(2596, 371);
  TestTeamScanRetVal<TEST_EXECSPACE, int64_t>{}(2596, 987);
  TestTeamScanRetVal<TEST_EXECSPACE, double>{}(2596, 1311);
}
#endif

}  // namespace Test
