/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
    // Set team size explicitly to
    // a) check whether this works in CPU backends with team_size > 1 and
    // b) make sure we have a power of 2 and for GPU backends due to limitation
    // of the scan algorithm implemented in CUDA etc.
    int team_size = 1;
    if (ExecutionSpace().concurrency() > 2) {
      if (ExecutionSpace().concurrency() > 10000)
        team_size = 128;
      else
        team_size = 3;
    }
    Kokkos::parallel_for(policy_type(M, team_size), *this);

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

}  // namespace Test
