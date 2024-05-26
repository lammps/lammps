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
#include <gtest/gtest.h>

namespace {

// Extended lambdas in parallel_for and parallel_reduce will not compile if
// KOKKOS_ENABLE_CUDA_LAMBDA is off
#if !defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_CUDA_LAMBDA)

struct TeamTeamCombinedReducer {
 public:
  void test_team_thread_range_only_scalars(const int n) {
    auto policy = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO);
    using team_member_type = decltype(policy)::member_type;

    auto teamView = Kokkos::View<int[4], TEST_EXECSPACE::memory_space>("view");

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(team_member_type const& team) {
          auto teamThreadRange = Kokkos::TeamThreadRange(team, n);
          int teamResult0, teamResult1, teamResult2, teamResult3;

          Kokkos::parallel_reduce(
              teamThreadRange,
              [=](int const& i, int& localVal0, int& localVal1, int& localVal2,
                  int& localVal3) {
                localVal0 += 1;
                localVal1 += i + 1;
                localVal2 += (i + 1) * n;
                localVal3 += n;
              },
              teamResult0, teamResult1, teamResult2, teamResult3);

          Kokkos::single(Kokkos::PerTeam(team), [=]() {
            teamView(0) = teamResult0;
            teamView(1) = teamResult1;
            teamView(2) = teamResult2;
            teamView(3) = teamResult3;
          });
        });

    auto hostView = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), teamView);

    if (n == 0) {
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(0));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(1));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(2));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(3));
    } else {
      EXPECT_EQ(n, hostView(0));
      EXPECT_EQ((n * (n + 1) / 2), hostView(1));
      EXPECT_EQ(n * n * (n + 1) / 2, hostView(2));
      EXPECT_EQ(n * n, hostView(3));
    }
  }

  void test_team_thread_range_only_builtin(const int n) {
    auto policy = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO);
    using team_member_type = decltype(policy)::member_type;

    auto teamView = Kokkos::View<int[4], TEST_EXECSPACE::memory_space>("view");

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(team_member_type const& team) {
          auto teamThreadRange = Kokkos::TeamThreadRange(team, n);
          int teamResult0, teamResult1, teamResult2, teamResult3;

          Kokkos::parallel_reduce(
              teamThreadRange,
              [=](int const& i, int& localVal0, int& localVal1, int& localVal2,
                  int& localVal3) {
                localVal0 += i + 1;
                localVal1 *= n;
                localVal2 = (localVal2 > (i + 1)) ? (i + 1) : localVal2;
                localVal3 = (localVal3 < (i + 1)) ? (i + 1) : localVal3;
              },
              Kokkos::Sum<int>(teamResult0), Kokkos::Prod<int>(teamResult1),
              Kokkos::Min<int>(teamResult2), Kokkos::Max<int>(teamResult3));

          Kokkos::single(Kokkos::PerTeam(team), [=]() {
            teamView(0) = teamResult0;
            teamView(1) = teamResult1;
            teamView(2) = teamResult2;
            teamView(3) = teamResult3;
          });
        });

    auto hostView = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), teamView);

    if (n == 0) {
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(0));
      EXPECT_EQ(Kokkos::reduction_identity<int>::prod(), hostView(1));
      EXPECT_EQ(Kokkos::reduction_identity<int>::min(), hostView(2));
      EXPECT_EQ(Kokkos::reduction_identity<int>::max(), hostView(3));
    } else {
      EXPECT_EQ((n * (n + 1) / 2), hostView(0));
      EXPECT_EQ(std::pow(n, n), hostView(1));
      EXPECT_EQ(1, hostView(2));
      EXPECT_EQ(n, hostView(3));
    }
  }

  void test_team_thread_range_combined_reducers(const int n) {
    auto policy = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO);
    using team_member_type = decltype(policy)::member_type;

    auto teamView = Kokkos::View<int*, TEST_EXECSPACE::memory_space>("view", 4);

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(team_member_type const& team) {
          auto teamThreadRange = Kokkos::TeamThreadRange(team, n);
          int teamResult0, teamResult1, teamResult2, teamResult3;

          Kokkos::parallel_reduce(
              teamThreadRange,
              [=](int const& i, int& localVal0, int& localVal1, int& localVal2,
                  int& localVal3) {
                localVal0 += i + 1;
                localVal1 += i + 1;
                localVal2 = (localVal2 < (i + 1)) ? (i + 1) : localVal2;
                localVal3 += n;
              },
              teamResult0, Kokkos::Sum<int>(teamResult1),
              Kokkos::Max<int>(teamResult2), teamResult3);

          Kokkos::single(Kokkos::PerTeam(team), [=]() {
            teamView(0) = teamResult0;
            teamView(1) = teamResult1;
            teamView(2) = teamResult2;
            teamView(3) = teamResult3;
          });
        });

    auto hostView = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), teamView);

    if (n == 0) {
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(0));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(1));
      EXPECT_EQ(Kokkos::reduction_identity<int>::max(), hostView(2));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(3));
    } else {
      EXPECT_EQ((n * (n + 1) / 2), hostView(0));
      EXPECT_EQ((n * (n + 1) / 2), hostView(1));
      EXPECT_EQ(n, hostView(2));
      EXPECT_EQ(n * n, hostView(3));
    }
  }

  void test_thread_vector_range_only_scalars(const int n) {
    auto policy = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO);
    using team_member_type = decltype(policy)::member_type;

    auto teamView = Kokkos::View<int[4], TEST_EXECSPACE::memory_space>("view");

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(team_member_type const& team) {
          auto teamThreadRange   = Kokkos::TeamThreadRange(team, 1);
          auto threadVectorRange = Kokkos::ThreadVectorRange(team, n);
          int teamResult0, teamResult1, teamResult2, teamResult3;

          Kokkos::parallel_for(teamThreadRange, [&](int const&) {
            Kokkos::parallel_reduce(
                threadVectorRange,
                [=](int const& i, int& localVal0, int& localVal1,
                    int& localVal2, int& localVal3) {
                  localVal0 += 1;
                  localVal1 += i + 1;
                  localVal2 += (i + 1) * n;
                  localVal3 += n;
                },
                teamResult0, teamResult1, teamResult2, teamResult3);

            teamView(0) = teamResult0;
            teamView(1) = teamResult1;
            teamView(2) = teamResult2;
            teamView(3) = teamResult3;
          });
        });

    auto hostView = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), teamView);

    if (n == 0) {
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(0));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(1));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(2));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(3));
    } else {
      EXPECT_EQ(n, hostView(0));
      EXPECT_EQ((n * (n + 1) / 2), hostView(1));
      EXPECT_EQ(n * n * (n + 1) / 2, hostView(2));
      EXPECT_EQ(n * n, hostView(3));
    }
  }

  void test_thread_vector_range_only_builtin(const int n) {
    auto policy = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO);
    using team_member_type = decltype(policy)::member_type;

    auto teamView = Kokkos::View<int[4], TEST_EXECSPACE::memory_space>("view");

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(team_member_type const& team) {
          auto teamThreadRange   = Kokkos::TeamThreadRange(team, 1);
          auto threadVectorRange = Kokkos::ThreadVectorRange(team, n);
          int teamResult0, teamResult1, teamResult2, teamResult3;

          Kokkos::parallel_for(teamThreadRange, [&](int const&) {
            Kokkos::parallel_reduce(
                threadVectorRange,
                [=](int const& i, int& localVal0, int& localVal1,
                    int& localVal2, int& localVal3) {
                  localVal0 += i + 1;
                  localVal1 *= n;
                  localVal2 = (localVal2 > (i + 1)) ? (i + 1) : localVal2;
                  localVal3 = (localVal3 < (i + 1)) ? (i + 1) : localVal3;
                },
                Kokkos::Sum<int>(teamResult0), Kokkos::Prod<int>(teamResult1),
                Kokkos::Min<int>(teamResult2), Kokkos::Max<int>(teamResult3));

            teamView(0) = teamResult0;
            teamView(1) = teamResult1;
            teamView(2) = teamResult2;
            teamView(3) = teamResult3;
          });
        });

    auto hostView = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), teamView);

    if (n == 0) {
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(0));
      EXPECT_EQ(Kokkos::reduction_identity<int>::prod(), hostView(1));
      EXPECT_EQ(Kokkos::reduction_identity<int>::min(), hostView(2));
      EXPECT_EQ(Kokkos::reduction_identity<int>::max(), hostView(3));
    } else {
      EXPECT_EQ((n * (n + 1) / 2), hostView(0));
      EXPECT_EQ(std::pow(n, n), hostView(1));
      EXPECT_EQ(1, hostView(2));
      EXPECT_EQ(n, hostView(3));
    }
  }

  void test_thread_vector_range_combined_reducers(const int n) {
    auto policy = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO);
    using team_member_type = decltype(policy)::member_type;

    auto teamView = Kokkos::View<int[4], TEST_EXECSPACE::memory_space>("view");

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(team_member_type const& team) {
          auto teamThreadRange   = Kokkos::TeamThreadRange(team, 1);
          auto threadVectorRange = Kokkos::ThreadVectorRange(team, n);
          int teamResult0, teamResult1, teamResult2, teamResult3;

          Kokkos::parallel_for(teamThreadRange, [&](int const&) {
            Kokkos::parallel_reduce(
                threadVectorRange,
                [=](int const& i, int& localVal0, int& localVal1,
                    int& localVal2, int& localVal3) {
                  localVal0 *= n;
                  localVal1 += i + 1;
                  localVal2 = (localVal2 > (i + 1)) ? (i + 1) : localVal2;
                  localVal3 += n;
                },
                Kokkos::Prod<int>(teamResult0), teamResult1,
                Kokkos::Min<int>(teamResult2), teamResult3);

            teamView(0) = teamResult0;
            teamView(1) = teamResult1;
            teamView(2) = teamResult2;
            teamView(3) = teamResult3;
          });
        });

    auto hostView = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), teamView);

    if (n == 0) {
      EXPECT_EQ(Kokkos::reduction_identity<int>::prod(), hostView(0));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(1));
      EXPECT_EQ(Kokkos::reduction_identity<int>::min(), hostView(2));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(3));
    } else {
      EXPECT_EQ(std::pow(n, n), hostView(0));
      EXPECT_EQ((n * (n + 1) / 2), hostView(1));
      EXPECT_EQ(1, hostView(2));
      EXPECT_EQ(n * n, hostView(3));
    }
  }

  void test_team_vector_range_only_scalars(const int n) {
    auto policy = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO);
    using team_member_type = decltype(policy)::member_type;

    auto teamView = Kokkos::View<int[4], TEST_EXECSPACE::memory_space>("view");

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(team_member_type const& team) {
          auto teamVectorRange = Kokkos::TeamVectorRange(team, n);
          int teamResult0, teamResult1, teamResult2, teamResult3;

          Kokkos::parallel_reduce(
              teamVectorRange,
              [=](int const& i, int& localVal0, int& localVal1, int& localVal2,
                  int& localVal3) {
                localVal0 += 1;
                localVal1 += i + 1;
                localVal2 += (i + 1) * n;
                localVal3 += n;
              },
              teamResult0, teamResult1, teamResult2, teamResult3);

          Kokkos::single(Kokkos::PerTeam(team), [=]() {
            teamView(0) = teamResult0;
            teamView(1) = teamResult1;
            teamView(2) = teamResult2;
            teamView(3) = teamResult3;
          });
        });

    auto hostView = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), teamView);

    if (n == 0) {
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(0));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(1));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(2));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(3));
    } else {
      EXPECT_EQ(n, hostView(0));
      EXPECT_EQ((n * (n + 1) / 2), hostView(1));
      EXPECT_EQ(n * n * (n + 1) / 2, hostView(2));
      EXPECT_EQ(n * n, hostView(3));
    }
  }

  void test_team_vector_range_only_builtin(const int n) {
    auto policy = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO);
    using team_member_type = decltype(policy)::member_type;

    auto teamView = Kokkos::View<int[4], TEST_EXECSPACE::memory_space>("view");

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(team_member_type const& team) {
          auto teamVectorRange = Kokkos::TeamVectorRange(team, n);
          int teamResult0, teamResult1, teamResult2, teamResult3;

          Kokkos::parallel_reduce(
              teamVectorRange,
              [=](int const& i, int& localVal0, int& localVal1, int& localVal2,
                  int& localVal3) {
                localVal0 += i + 1;
                localVal1 *= n;
                localVal2 = (localVal2 > (i + 1)) ? (i + 1) : localVal2;
                localVal3 = (localVal3 < (i + 1)) ? (i + 1) : localVal3;
              },
              Kokkos::Sum<int>(teamResult0), Kokkos::Prod<int>(teamResult1),
              Kokkos::Min<int>(teamResult2), Kokkos::Max<int>(teamResult3));

          Kokkos::single(Kokkos::PerTeam(team), [=]() {
            teamView(0) = teamResult0;
            teamView(1) = teamResult1;
            teamView(2) = teamResult2;
            teamView(3) = teamResult3;
          });
        });

    auto hostView = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), teamView);

    if (n == 0) {
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(0));
      EXPECT_EQ(Kokkos::reduction_identity<int>::prod(), hostView(1));
      EXPECT_EQ(Kokkos::reduction_identity<int>::min(), hostView(2));
      EXPECT_EQ(Kokkos::reduction_identity<int>::max(), hostView(3));
    } else {
      EXPECT_EQ((n * (n + 1) / 2), hostView(0));
      EXPECT_EQ(std::pow(n, n), hostView(1));
      EXPECT_EQ(1, hostView(2));
      EXPECT_EQ(n, hostView(3));
    }
  }

  void test_team_vector_range_combined_reducers(const int n) {
    auto policy = Kokkos::TeamPolicy<TEST_EXECSPACE>(1, Kokkos::AUTO);
    using team_member_type = decltype(policy)::member_type;

    auto teamView = Kokkos::View<int[4], TEST_EXECSPACE::memory_space>("view");

    Kokkos::parallel_for(
        policy, KOKKOS_LAMBDA(team_member_type const& team) {
          auto teamVectorRange = Kokkos::TeamVectorRange(team, n);
          int teamResult0, teamResult1, teamResult2, teamResult3;

          Kokkos::parallel_reduce(
              teamVectorRange,
              [=](int const& i, int& localVal0, int& localVal1, int& localVal2,
                  int& localVal3) {
                localVal0 += i + 1;
                localVal1 += i + 1;
                localVal2 = (localVal2 < (i + 1)) ? (i + 1) : localVal2;
                localVal3 += n;
              },
              teamResult0, Kokkos::Sum<int>(teamResult1),
              Kokkos::Max<int>(teamResult2), teamResult3);

          Kokkos::single(Kokkos::PerTeam(team), [=]() {
            teamView(0) = teamResult0;
            teamView(1) = teamResult1;
            teamView(2) = teamResult2;
            teamView(3) = teamResult3;
          });
        });

    auto hostView = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), teamView);

    if (n == 0) {
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(0));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(1));
      EXPECT_EQ(Kokkos::reduction_identity<int>::max(), hostView(2));
      EXPECT_EQ(Kokkos::reduction_identity<int>::sum(), hostView(3));
    } else {
      EXPECT_EQ((n * (n + 1) / 2), hostView(0));
      EXPECT_EQ((n * (n + 1) / 2), hostView(1));
      EXPECT_EQ(n, hostView(2));
      EXPECT_EQ(n * n, hostView(3));
    }
  }
};

TEST(TEST_CATEGORY, team_thread_range_combined_reducers) {
#if defined(KOKKOS_ENABLE_OPENMPTARGET)
  if constexpr (std::is_same_v<TEST_EXECSPACE,
                               Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "team_reduce with a generic reducer is not implemented for "
                 << TEST_EXECSPACE::name();

#elif defined(KOKKOS_ENABLE_OPENACC)
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenACC>)
    GTEST_SKIP() << "team_reduce with a generic reducer is not implemented for "
                 << TEST_EXECSPACE::name();
#endif

  TeamTeamCombinedReducer tester;
  tester.test_team_thread_range_only_scalars(5);
  tester.test_team_thread_range_only_builtin(7);
  tester.test_team_thread_range_combined_reducers(0);
  tester.test_team_thread_range_combined_reducers(9);
}

TEST(TEST_CATEGORY, thread_vector_range_combined_reducers) {
#if defined(KOKKOS_ENABLE_OPENMPTARGET)
  if constexpr (std::is_same_v<TEST_EXECSPACE,
                               Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "team_reduce with a generic reducer is not implemented for "
                 << TEST_EXECSPACE::name();

#elif defined(KOKKOS_ENABLE_OPENACC)
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenACC>)
    GTEST_SKIP() << "team_reduce with a generic reducer is not implemented for "
                 << TEST_EXECSPACE::name();
#endif

  TeamTeamCombinedReducer tester;
  tester.test_thread_vector_range_only_scalars(5);
  tester.test_thread_vector_range_only_builtin(7);
  tester.test_thread_vector_range_combined_reducers(0);
  tester.test_thread_vector_range_combined_reducers(9);
}

TEST(TEST_CATEGORY, team_vector_range_combined_reducers) {
#ifdef KOKKOS_ENABLE_OPENMPTARGET  // FIXME_OPENMPTARGET
  if constexpr (std::is_same_v<TEST_EXECSPACE,
                               Kokkos::Experimental::OpenMPTarget>)
    GTEST_SKIP() << "team_reduce with a generic reducer is not implemented for "
                 << TEST_EXECSPACE::name();
#endif

#ifdef KOKKOS_ENABLE_OPENACC  // FIXME_OPENACC
  if constexpr (std::is_same_v<TEST_EXECSPACE, Kokkos::Experimental::OpenACC>)
    GTEST_SKIP() << "team_reduce with a generic reducer is not implemented for "
                 << TEST_EXECSPACE::name();
#endif

  TeamTeamCombinedReducer tester;
  tester.test_team_vector_range_only_scalars(5);
  tester.test_team_vector_range_only_builtin(7);
  tester.test_team_vector_range_combined_reducers(0);
  tester.test_team_vector_range_combined_reducers(9);
}

#endif

}  // namespace
