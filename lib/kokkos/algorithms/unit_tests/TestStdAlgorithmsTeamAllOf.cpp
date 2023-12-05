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

#include <TestStdAlgorithmsCommon.hpp>

namespace Test {
namespace stdalgos {
namespace TeamAllOf {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct GreaterThanValueFunctor {
  ValueType m_val;

  KOKKOS_INLINE_FUNCTION
  GreaterThanValueFunctor(ValueType val) : m_val(val) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(ValueType val) const { return (val > m_val); }
};

template <class DataViewType, class AlllOfResultsViewType,
          class IntraTeamSentinelView, class UnaryPredType>
struct TestFunctorA {
  DataViewType m_dataView;
  AlllOfResultsViewType m_allOfResultsView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;
  UnaryPredType m_unaryPred;

  TestFunctorA(const DataViewType dataView,
               const AlllOfResultsViewType allOfResultsView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick,
               UnaryPredType unaryPred)
      : m_dataView(dataView),
        m_allOfResultsView(allOfResultsView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick),
        m_unaryPred(unaryPred) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();

    auto myRowViewFrom = Kokkos::subview(m_dataView, myRowIndex, Kokkos::ALL());
    bool result        = false;

    switch (m_apiPick) {
      case 0: {
        result = KE::all_of(member, KE::cbegin(myRowViewFrom),
                            KE::cend(myRowViewFrom), m_unaryPred);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_allOfResultsView(myRowIndex) = result;
        });
        break;
      }

      case 1: {
        result = KE::all_of(member, myRowViewFrom, m_unaryPred);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_allOfResultsView(myRowIndex) = result;
        });
        break;
      }
    }

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck = team_members_have_matching_result(
        member, result, m_allOfResultsView(myRowIndex));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(myRowIndex) = intraTeamCheck;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level all_of
   */

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // create a view in the memory space associated with default exespace
  // with as many rows as the number of teams and fill it with random
  // values from an arbitrary range.
  constexpr ValueType lowerBound = 5;
  constexpr ValueType upperBound = 523;
  const auto bounds              = make_bounds(lowerBound, upperBound);

  auto [dataView, dataViewBeforeOp_h] = create_random_view_and_host_clone(
      LayoutTag{}, numTeams, numCols, bounds, "dataView");

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // to verify that things work, each team stores the result of its all_of call,
  // and then we check that these match what we expect
  Kokkos::View<bool*> allOfResultsView("allOfResultsView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  GreaterThanValueFunctor unaryPred{lowerBound - 1};

  // use CTAD for functor
  TestFunctorA fnc(dataView, allOfResultsView, intraTeamSentinelView, apiId,
                   unaryPred);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel and check
  // -----------------------------------------------
  auto allOfResultsView_h      = create_host_space_copy(allOfResultsView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);

  for (std::size_t i = 0; i < dataView.extent(0); ++i) {
    auto rowFrom = Kokkos::subview(dataViewBeforeOp_h, i, Kokkos::ALL());
    const bool result =
        std::all_of(KE::cbegin(rowFrom), KE::cend(rowFrom), unaryPred);
    ASSERT_EQ(result, allOfResultsView_h(i));
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios() {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_all_of_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamAllOf
}  // namespace stdalgos
}  // namespace Test
