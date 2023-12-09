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
namespace TeamCountIf {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct GreaterThanValueFunctor {
  ValueType m_val;

  KOKKOS_INLINE_FUNCTION
  GreaterThanValueFunctor(ValueType val) : m_val(val) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(ValueType val) const { return (val > m_val); }
};

template <class ViewType, class CountsViewType, class IntraTeamSentinelView,
          class ValueType>
struct TestFunctorA {
  ViewType m_view;
  CountsViewType m_countsView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  ValueType m_threshold;
  int m_apiPick;

  TestFunctorA(const ViewType view, const CountsViewType countsView,
               const IntraTeamSentinelView intraTeamSentinelView,
               ValueType threshold, int apiPick)
      : m_view(view),
        m_countsView(countsView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_threshold(threshold),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowView        = Kokkos::subview(m_view, myRowIndex, Kokkos::ALL());
    std::size_t myCount   = 0;

    GreaterThanValueFunctor predicate(m_threshold);
    if (m_apiPick == 0) {
      myCount = KE::count_if(member, KE::begin(myRowView), KE::end(myRowView),
                             predicate);

      Kokkos::single(Kokkos::PerTeam(member),
                     [=, *this]() { m_countsView(myRowIndex) = myCount; });
    } else if (m_apiPick == 1) {
      myCount = KE::count_if(member, myRowView, predicate);
      Kokkos::single(Kokkos::PerTeam(member),
                     [=, *this]() { m_countsView(myRowIndex) = myCount; });
    }

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck = team_members_have_matching_result(
        member, myCount, m_countsView(myRowIndex));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(myRowIndex) = intraTeamCheck;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level count_if where only the values
     strictly greater than a threshold are counted
   */

  const auto threshold = static_cast<ValueType>(151);

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // create a view in the memory space associated with default exespace
  // with as many rows as the number of teams and fill it with random
  // values from an arbitrary range.
  constexpr ValueType lowerBound = 5;
  constexpr ValueType upperBound = 523;
  const auto bounds              = make_bounds(lowerBound, upperBound);

  auto [dataView, cloneOfDataViewBeforeOp_h] =
      create_random_view_and_host_clone(LayoutTag{}, numTeams, numCols, bounds,
                                        "dataView");

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // to verify that things work, each team stores the result
  // of its count_if call, and then we check
  // that these match what we expect
  Kokkos::View<std::size_t*> countsView("countsView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(dataView, countsView, intraTeamSentinelView, threshold,
                   apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  auto countsView_h            = create_host_space_copy(countsView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  for (std::size_t i = 0; i < cloneOfDataViewBeforeOp_h.extent(0); ++i) {
    std::size_t goldCountForRow = 0;
    for (std::size_t j = 0; j < cloneOfDataViewBeforeOp_h.extent(1); ++j) {
      if (cloneOfDataViewBeforeOp_h(i, j) > threshold) {
        goldCountForRow++;
      }
    }
    ASSERT_EQ(goldCountForRow, countsView_h(i));
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

TEST(std_algorithms_count_if_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamCountIf
}  // namespace stdalgos
}  // namespace Test
