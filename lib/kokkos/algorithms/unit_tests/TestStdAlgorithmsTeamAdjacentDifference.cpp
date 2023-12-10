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
namespace TeamAdjacentDifference {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct PlusFunctor {
  KOKKOS_INLINE_FUNCTION constexpr ValueType operator()(
      const ValueType& lhs, const ValueType& rhs) const {
    return lhs + rhs;
  }
};

template <class SourceViewType, class DestViewType, class DistancesViewType,
          class IntraTeamSentinelView, class BinaryOp>
struct TestFunctorA {
  SourceViewType m_sourceView;
  DestViewType m_destView;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;
  BinaryOp m_binaryOp;

  TestFunctorA(const SourceViewType sourceView, const DestViewType destView,
               const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick,
               BinaryOp binaryOp)
      : m_sourceView(sourceView),
        m_destView(destView),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick),
        m_binaryOp(std::move(binaryOp)) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowViewFrom =
        Kokkos::subview(m_sourceView, myRowIndex, Kokkos::ALL());
    auto myRowViewDest = Kokkos::subview(m_destView, myRowIndex, Kokkos::ALL());
    ptrdiff_t resultDist = 0;

    switch (m_apiPick) {
      case 0: {
        auto it    = KE::adjacent_difference(member, KE::cbegin(myRowViewFrom),
                                          KE::cend(myRowViewFrom),
                                          KE::begin(myRowViewDest));
        resultDist = KE::distance(KE::begin(myRowViewDest), it);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_distancesView(myRowIndex) = resultDist;
        });
        break;
      }

      case 1: {
        auto it    = KE::adjacent_difference(member, KE::cbegin(myRowViewFrom),
                                          KE::cend(myRowViewFrom),
                                          KE::begin(myRowViewDest), m_binaryOp);
        resultDist = KE::distance(KE::begin(myRowViewDest), it);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_distancesView(myRowIndex) = resultDist;
        });
        break;
      }

      case 2: {
        auto it = KE::adjacent_difference(member, myRowViewFrom, myRowViewDest);
        resultDist = KE::distance(KE::begin(myRowViewDest), it);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_distancesView(myRowIndex) = resultDist;
        });
        break;
      }

      case 3: {
        auto it = KE::adjacent_difference(member, myRowViewFrom, myRowViewDest,
                                          m_binaryOp);
        resultDist = KE::distance(KE::begin(myRowViewDest), it);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_distancesView(myRowIndex) = resultDist;
        });
        break;
      }
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
     use a rank-2 view randomly filled with values,
     and run a team-level adjacent_difference
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

  auto [sourceView, sourceViewBeforeOp_h] = create_random_view_and_host_clone(
      LayoutTag{}, numTeams, numCols, bounds, "sourceView");

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // to verify that things work, each team stores the result
  // of its adjacent_difference call, and then we check
  // that these match what we expect
  Kokkos::View<ValueType**> destView("destView", numTeams, numCols);

  // adjacent_difference returns an iterator so to verify that it is correct
  // each team stores the distance of the returned iterator from the beginning
  // of the interval that team operates on and then we check that these
  // distances match the std result
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(sourceView, destView, distancesView, intraTeamSentinelView,
                   apiId, PlusFunctor<ValueType>{});
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel and check
  // -----------------------------------------------
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  Kokkos::View<ValueType**, Kokkos::HostSpace> stdDestView("stdDestView",
                                                           numTeams, numCols);

  for (std::size_t i = 0; i < sourceView.extent(0); ++i) {
    auto rowFrom = Kokkos::subview(sourceViewBeforeOp_h, i, Kokkos::ALL());
    auto rowDest = Kokkos::subview(stdDestView, i, Kokkos::ALL());
    ASSERT_TRUE(intraTeamSentinelView_h(i));

    switch (apiId) {
      case 0:
      case 2: {
        auto it = std::adjacent_difference(KE::begin(rowFrom), KE::end(rowFrom),
                                           KE::begin(rowDest));

        const std::size_t stdDistance = KE::distance(KE::begin(rowDest), it);
        ASSERT_EQ(stdDistance, distancesView_h(i));
        break;
      }

      case 1:
      case 3: {
        auto it = std::adjacent_difference(KE::begin(rowFrom), KE::end(rowFrom),
                                           KE::begin(rowDest),
                                           PlusFunctor<ValueType>{});

        const std::size_t stdDistance = KE::distance(KE::begin(rowDest), it);
        ASSERT_EQ(stdDistance, distancesView_h(i));
        break;
      }
    }
  }

  auto dataViewAfterOp_h = create_host_space_copy(destView);
  expect_equal_host_views(stdDestView, dataViewAfterOp_h);
}

template <class LayoutTag, class ValueType>
void run_all_scenarios() {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1, 2, 3}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_adjacent_difference_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamAdjacentDifference
}  // namespace stdalgos
}  // namespace Test
