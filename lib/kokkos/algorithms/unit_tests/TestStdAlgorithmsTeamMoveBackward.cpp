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
namespace TeamMovebackward {

namespace KE = Kokkos::Experimental;

template <class SourceViewType, class DestViewType, class DistancesViewType,
          class IntraTeamSentinelView>
struct TestFunctorA {
  SourceViewType m_sourceView;
  DestViewType m_destView;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const SourceViewType sourceView, const DestViewType destView,
               const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_sourceView(sourceView),
        m_destView(destView),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowViewFrom =
        Kokkos::subview(m_sourceView, myRowIndex, Kokkos::ALL());
    auto myRowViewDest = Kokkos::subview(m_destView, myRowIndex, Kokkos::ALL());
    ptrdiff_t resultDist = 0;

    if (m_apiPick == 0) {
      auto it =
          KE::move_backward(member, KE::cbegin(myRowViewFrom),
                            KE::cend(myRowViewFrom), KE::end(myRowViewDest));
      resultDist = KE::distance(KE::begin(myRowViewDest), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    } else if (m_apiPick == 1) {
      auto it    = KE::move_backward(member, myRowViewFrom, myRowViewDest);
      resultDist = KE::distance(KE::begin(myRowViewDest), it);
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
     randomly fill a source view and do team KE::move_backward
     into a destination view. The operation is done via a
     team parfor with one row per team.
   */

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // create a view in the memory space associated with default exespace
  // with as many rows as the number of teams and fill it with random
  // values from an arbitrary range
  auto [sourceView, cloneOfSourceViewBeforeOp_h] =
      create_random_view_and_host_clone(
          LayoutTag{}, numTeams, numCols,
          Kokkos::pair<ValueType, ValueType>{11, 523}, "sourceView");

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());
  // create the destination view: for a meaningful test, the destination
  // view must have more columns that than the source view so that we
  // can check that the elements are copied into the right place
  constexpr std::size_t extra = 10;
  Kokkos::View<ValueType**> destView("destView", numTeams, numCols + extra);
  // make host copy of destView that should remain unchanged
  auto destViewBeforeOp_h = create_host_space_copy(destView);

  // each team stores the distance of the returned iterator from the
  // beginning of the interval that team operates on and then we check
  // that these distances match the expectation
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(sourceView, destView, distancesView, intraTeamSentinelView,
                   apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  // KE::move_backwar has been applied to sourceView, so we should
  // NOT use sourceView henceforth, becuase all its elements
  // have been moved from.

  auto distancesView_h         = create_host_space_copy(distancesView);
  auto destViewAfterOp_h       = create_host_space_copy(destView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  for (std::size_t i = 0; i < destViewAfterOp_h.extent(0); ++i) {
    // first extra num of columns should be unchanged
    for (std::size_t j = 0; j < extra; ++j) {
      EXPECT_TRUE(destViewAfterOp_h(i, j) == destViewBeforeOp_h(i, j));
    }

    // after extra # of column (inclusive) should match the source view
    for (std::size_t j = extra; j < destViewBeforeOp_h.extent(1); ++j) {
      EXPECT_TRUE(cloneOfSourceViewBeforeOp_h(i, j - extra) ==
                  destViewAfterOp_h(i, j));
    }

    // each team should have returned an interator whose distance
    // from the beginning of the row should satisfy this
    EXPECT_TRUE(distancesView_h(i) == extra);
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

TEST(std_algorithms_move_backward_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamMovebackward
}  // namespace stdalgos
}  // namespace Test
