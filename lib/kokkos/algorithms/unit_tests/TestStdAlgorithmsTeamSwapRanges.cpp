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
namespace TeamSwapRanges {

namespace KE = Kokkos::Experimental;

template <class View1Type, class View2Type, class DistancesViewType,
          class IntraTeamSentinelView>
struct TestFunctorA {
  View1Type m_view1;
  View2Type m_view2;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const View1Type view1, const View2Type view2,
               const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_view1(view1),
        m_view2(view2),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowView1       = Kokkos::subview(m_view1, myRowIndex, Kokkos::ALL());
    auto myRowView2       = Kokkos::subview(m_view2, myRowIndex, Kokkos::ALL());
    ptrdiff_t resultDist  = 0;

    if (m_apiPick == 0) {
      auto it    = KE::swap_ranges(member, KE::begin(myRowView1),
                                KE::end(myRowView1), KE::begin(myRowView2));
      resultDist = KE::distance(KE::begin(myRowView2), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    } else if (m_apiPick == 1) {
      auto it    = KE::swap_ranges(member, myRowView1, myRowView2);
      resultDist = KE::distance(KE::begin(myRowView2), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    }

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck = team_members_have_matching_result(
        member, resultDist, m_distancesView(myRowIndex));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(myRowIndex) = intraTeamCheck;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId) {
  /* description:
     randomly fill two views and do team level swap_ranges
   */

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  auto [dataView1, cloneOfDataView1BeforeOp_h] =
      create_random_view_and_host_clone(
          LayoutTag{}, numTeams, numCols,
          Kokkos::pair<ValueType, ValueType>{11, 523}, "dataView1");

  auto [dataView2, cloneOfDataView2BeforeOp_h] =
      create_random_view_and_host_clone(
          LayoutTag{}, numTeams, numCols,
          Kokkos::pair<ValueType, ValueType>{530, 1523}, "dataView2");

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // each team stores the distance of the returned iterator from the
  // beginning of the interval that team operates on and then we check
  // that these distances match the expectation
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(dataView1, dataView2, distancesView, intraTeamSentinelView,
                   apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto dataView1AfterOp_h      = create_host_space_copy(dataView1);
  auto dataView2AfterOp_h      = create_host_space_copy(dataView2);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);

  for (std::size_t i = 0; i < dataView1AfterOp_h.extent(0); ++i) {
    for (std::size_t j = 0; j < dataView1AfterOp_h.extent(1); ++j) {
      ASSERT_EQ(cloneOfDataView1BeforeOp_h(i, j), dataView2AfterOp_h(i, j));
      ASSERT_EQ(cloneOfDataView2BeforeOp_h(i, j), dataView1AfterOp_h(i, j));
    }
    // each team should return an iterator past the last column
    EXPECT_TRUE(distancesView_h(i) == numCols);
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios() {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 11113}) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_swap_ranges_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamSwapRanges
}  // namespace stdalgos
}  // namespace Test
