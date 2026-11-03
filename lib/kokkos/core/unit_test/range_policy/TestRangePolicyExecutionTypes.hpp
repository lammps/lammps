// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>

namespace Test {

template <class Policy>
KOKKOS_INLINE_FUNCTION int check_runtime_inputs(
    Policy& p, const typename Policy::index_type expected_begin,
    const typename Policy::index_type expected_end,
    const typename Policy::index_type chunk_size = 1) {
  int nerrs = 0;

  if (p.begin() != expected_begin) ++nerrs;
  if (p.end() != expected_end) ++nerrs;

  auto p2 = p.set_chunk_size(chunk_size);
  if (p2.chunk_size() != chunk_size) ++nerrs;

  return nerrs;
}

void test_self_similar_range_policy_runtime() {
  using IndexType = typename Kokkos::DefaultExecutionSpace::size_type;

  IndexType beg        = 5;
  IndexType end        = 15;
  IndexType chunk_size = 10;

  auto p_execspace =
      Kokkos::RangePolicy(Kokkos::DefaultExecutionSpace(), beg, end);
  auto nerrs_exec_space =
      check_runtime_inputs(p_execspace, beg, end, chunk_size);
  ASSERT_EQ(nerrs_exec_space, 0);

  int nerrs_team_handle;
  using team_t = typename Kokkos::TeamPolicy<>::member_type;
  Kokkos::parallel_reduce(
      "check_runtime", Kokkos::TeamPolicy(1, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const team_t& team, int& nerrs) {
        auto p_teamhandle = Kokkos::RangePolicy(team, beg, end);
        auto tvr          = Kokkos::TeamVectorRange(team, beg, end);
        nerrs = check_runtime_inputs(p_teamhandle, tvr.start, tvr.end);
      },
      nerrs_team_handle);
  ASSERT_EQ(nerrs_team_handle, 0);
}

template <class Exec, class X, class Y>
KOKKOS_INLINE_FUNCTION void sum_views(const Exec& exec, const X& x,
                                      const Y& y) {
  auto policy = Kokkos::RangePolicy(exec, 0, x.extent(0));
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const int& i) { x(i) += y(i); });
}

void test_self_similar_range_policy_computation() {
  size_t N         = 7;
  size_t num_teams = 5;

  Kokkos::View<float*> v_x("v_x", N), v_y("v_y", N);
  Kokkos::View<float**> M_x("M_x", num_teams, N), M_y("M_y", num_teams, N);

  // Initialize v_x and v_y with values from 1 to N
  Kokkos::parallel_for(
      "init_v_x", Kokkos::RangePolicy<>(0, N),
      KOKKOS_LAMBDA(const int& i) { v_x(i) = static_cast<float>(i + 1); });
  Kokkos::parallel_for(
      "init_v_y", Kokkos::RangePolicy<>(0, N),
      KOKKOS_LAMBDA(const int& i) { v_y(i) = static_cast<float>(i + 1); });

  // Initialize M_x and M_y with values from 1 to M (flattened index)
  Kokkos::parallel_for(
      "init_M_x", Kokkos::RangePolicy<>(0, num_teams),
      KOKKOS_LAMBDA(const int& i) {
        for (size_t j = 0; j < N; j++) {
          M_x(i, j) = static_cast<float>(i * N + j + 1);
        }
      });
  Kokkos::parallel_for(
      "init_M_y", Kokkos::RangePolicy<>(0, num_teams),
      KOKKOS_LAMBDA(const int& i) {
        for (size_t j = 0; j < N; j++) {
          M_y(i, j) = static_cast<float>(i * N + j + 1);
        }
      });

  // Call sum_views(ExecSpace):
  sum_views(Kokkos::DefaultExecutionSpace(), v_x, v_y);

  // Call sum_views(TeamHandle)
  using team_t = typename Kokkos::TeamPolicy<>::member_type;
  Kokkos::parallel_for(
      "apxyFromTeam", Kokkos::TeamPolicy(num_teams, Kokkos::AUTO()),
      KOKKOS_LAMBDA(const team_t& team) {
        sum_views(team, Kokkos::subview(M_x, team.league_rank(), Kokkos::ALL()),
                  Kokkos::subview(M_y, team.league_rank(), Kokkos::ALL()));
      });

  // Check v_x
  size_t result = 0;
  Kokkos::parallel_reduce(
      "Check1", v_x.extent(0),
      KOKKOS_LAMBDA(int i, size_t& val) { val += v_x(i); }, result);
  size_t expected_v_x = N * (N + 1);
  ASSERT_EQ(result, expected_v_x);

  // Check individual elements of v_x
  Kokkos::parallel_reduce(
      "Check1_elements", v_x.extent(0),
      KOKKOS_LAMBDA(int i, size_t& errors) {
        float expected = static_cast<float>(2 * (i + 1));
        if (v_x(i) != expected) ++errors;
      },
      result);
  ASSERT_EQ(result, size_t(0));

  // Check M_x
  result = 0;
  Kokkos::parallel_reduce(
      "Check2", M_x.extent(0),
      KOKKOS_LAMBDA(int i, size_t& val) {
        for (int j = 0; j < M_x.extent_int(1); j++) val += M_x(i, j);
      },
      result);
  size_t M_total      = num_teams * N;
  size_t expected_M_x = M_total * (M_total + 1);
  ASSERT_EQ(result, expected_M_x);

  // Check individual elements of M_x
  Kokkos::parallel_reduce(
      "Check2_elements", M_x.extent(0),
      KOKKOS_LAMBDA(int i, size_t& errors) {
        for (int j = 0; j < M_x.extent_int(1); j++) {
          float expected = static_cast<float>(2 * (i * N + j + 1));
          if (M_x(i, j) != expected) ++errors;
        }
      },
      result);
  ASSERT_EQ(result, size_t(0));
}

TEST(TEST_CATEGORY, self_similar_range_policy_runtime) {
  test_self_similar_range_policy_runtime();
}

TEST(TEST_CATEGORY, self_similar_range_policy_computation) {
  test_self_similar_range_policy_computation();
}

}  // namespace Test
