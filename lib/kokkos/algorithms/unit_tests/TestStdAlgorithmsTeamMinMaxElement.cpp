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
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace TeamMinMaxElement {

namespace KE = Kokkos::Experimental;

template <class ViewType, class DistancesViewType, class IntraTeamSentinelView>
struct TestFunctorA {
  ViewType m_view;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const ViewType view, const DistancesViewType distancesView,
               IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_view(view),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowView        = Kokkos::subview(m_view, myRowIndex, Kokkos::ALL());
    ptrdiff_t resultDist1 = 0;
    ptrdiff_t resultDist2 = 0;

    if (m_apiPick == 0) {
      auto itPair = KE::minmax_element(member, KE::cbegin(myRowView),
                                       KE::cend(myRowView));
      resultDist1 = KE::distance(KE::cbegin(myRowView), itPair.first);
      resultDist2 = KE::distance(KE::cbegin(myRowView), itPair.second);

      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex, 0) = resultDist1;
        m_distancesView(myRowIndex, 1) = resultDist2;
      });
    }

    else if (m_apiPick == 1) {
      auto itPair = KE::minmax_element(member, myRowView);
      resultDist1 = KE::distance(KE::begin(myRowView), itPair.first);
      resultDist2 = KE::distance(KE::begin(myRowView), itPair.second);

      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex, 0) = resultDist1;
        m_distancesView(myRowIndex, 1) = resultDist2;
      });
    }
#if not defined KOKKOS_ENABLE_OPENMPTARGET
    else if (m_apiPick == 2) {
      using value_type = typename ViewType::value_type;
      auto itPair =
          KE::minmax_element(member, KE::cbegin(myRowView), KE::cend(myRowView),
                             CustomLessThanComparator<value_type>{});
      resultDist1 = KE::distance(KE::cbegin(myRowView), itPair.first);
      resultDist2 = KE::distance(KE::cbegin(myRowView), itPair.second);

      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex, 0) = resultDist1;
        m_distancesView(myRowIndex, 1) = resultDist2;
      });
    }

    else if (m_apiPick == 3) {
      using value_type = typename ViewType::value_type;
      auto itPair      = KE::minmax_element(member, myRowView,
                                       CustomLessThanComparator<value_type>{});
      resultDist1      = KE::distance(KE::begin(myRowView), itPair.first);
      resultDist2      = KE::distance(KE::begin(myRowView), itPair.second);

      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex, 0) = resultDist1;
        m_distancesView(myRowIndex, 1) = resultDist2;
      });
    }
#endif

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck1 = team_members_have_matching_result(
        member, resultDist1, m_distancesView(myRowIndex, 0));
    const bool intraTeamCheck2 = team_members_have_matching_result(
        member, resultDist2, m_distancesView(myRowIndex, 1));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(myRowIndex) = intraTeamCheck1 && intraTeamCheck2;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId) {
  /* description:
     team-level KE::minmax_element on a rank-2 view where
     data is filled randomly and we use one team per row.
   */

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // create a view in the memory space associated with default exespace
  // with as many rows as the number of teams and fill it with random
  auto [dataView, cloneOfDataViewBeforeOp_h] =
      create_random_view_and_host_clone(
          LayoutTag{}, numTeams, numCols,
          Kokkos::pair<ValueType, ValueType>{0, 1153}, "dataView");

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // each team stores the distance of the returned value from the
  // beginning of the interval that team operates on and then we check
  // that these distances match the expectation
  Kokkos::View<std::size_t**> distancesView("distancesView", numTeams, 2);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(dataView, distancesView, intraTeamSentinelView, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run std algo and check
  // -----------------------------------------------
  // here I can use cloneOfDataViewBeforeOp_h to run std algo on
  // since that contains a valid copy of the data
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  auto dataViewAfterOp_h       = create_host_space_copy(dataView);
  for (std::size_t i = 0; i < cloneOfDataViewBeforeOp_h.extent(0); ++i) {
    auto myRow = Kokkos::subview(cloneOfDataViewBeforeOp_h, i, Kokkos::ALL());

    std::size_t stdDistance[2];
    if (apiId <= 1) {
      auto itPair    = std::minmax_element(KE::cbegin(myRow), KE::cend(myRow));
      stdDistance[0] = KE::distance(KE::cbegin(myRow), itPair.first);
      stdDistance[1] = KE::distance(KE::cbegin(myRow), itPair.second);
    } else {
      auto itPair    = std::minmax_element(KE::cbegin(myRow), KE::cend(myRow),
                                        CustomLessThanComparator<value_type>{});
      stdDistance[0] = KE::distance(KE::cbegin(myRow), itPair.first);
      stdDistance[1] = KE::distance(KE::cbegin(myRow), itPair.second);
    }

    ASSERT_EQ(stdDistance[0], distancesView_h(i, 0));
    ASSERT_EQ(stdDistance[1], distancesView_h(i, 1));
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }

  // dataView should remain unchanged
  expect_equal_host_views(cloneOfDataViewBeforeOp_h, dataViewAfterOp_h);
}

template <class LayoutTag, class ValueType>
void run_all_scenarios() {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 5113}) {
      // for OpenMPTarget we need to avod api accepting a custom
      // comparator because it is not supported
      for (int apiId : {0, 1, 2, 3}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_minmax_element_team_test, test) {
#if not defined KOKKOS_ENABLE_OPENMPTARGET
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedTwoRowsTag, double>();
  run_all_scenarios<StridedThreeRowsTag, int>();
#endif
}

}  // namespace TeamMinMaxElement
}  // namespace stdalgos
}  // namespace Test
