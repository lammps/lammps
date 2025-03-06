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
namespace TeamTransformBinaryOp {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct AddValuesBinaryOp {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& a, const ValueType& b) const {
    return a + b;
  }
};

template <class SourceView1Type, class SourceView2Type, class DestViewType,
          class DistancesViewType, class IntraTeamSentinelView>
struct TestFunctorA {
  SourceView1Type m_sourceView1;
  SourceView2Type m_sourceView2;
  DestViewType m_destView;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const SourceView1Type sourceView1,
               const SourceView2Type sourceView2, const DestViewType destView,
               const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_sourceView1(sourceView1),
        m_sourceView2(sourceView2),
        m_destView(destView),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowView1From =
        Kokkos::subview(m_sourceView1, myRowIndex, Kokkos::ALL());
    auto myRowView2From =
        Kokkos::subview(m_sourceView2, myRowIndex, Kokkos::ALL());
    auto myRowViewDest = Kokkos::subview(m_destView, myRowIndex, Kokkos::ALL());
    ptrdiff_t resultDist = 0;

    using value_type = typename SourceView1Type::value_type;
    if (m_apiPick == 0) {
      auto it = KE::transform(
          member, KE::cbegin(myRowView1From), KE::cend(myRowView1From),
          KE::cbegin(myRowView2From), KE::begin(myRowViewDest),
          AddValuesBinaryOp<value_type>());

      resultDist = KE::distance(KE::begin(myRowViewDest), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    } else if (m_apiPick == 1) {
      auto it = KE::transform(member, myRowView1From, myRowView2From,
                              myRowViewDest, AddValuesBinaryOp<value_type>());

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
     team level transform with each team handling a row of
     two rank-2 source views and applying a binary op that
     add each pair of element from those two views
   */

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // create a view in the memory space associated with default exespace
  // with as many rows as the number of teams and fill it with random
  // values from an arbitrary range
  auto [sourceView1, cloneOfSourceView1BeforeOp_h] =
      create_random_view_and_host_clone(
          LayoutTag{}, numTeams, numCols,
          Kokkos::pair<ValueType, ValueType>{0, 523}, "sourceView1",
          317539 /*random seed*/);
  auto [sourceView2, cloneOfSourceView2BeforeOp_h] =
      create_random_view_and_host_clone(
          LayoutTag{}, numTeams, numCols,
          Kokkos::pair<ValueType, ValueType>{0, 523}, "sourceView2",
          957313 /*random seed*/);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());
  // create the destination view
  Kokkos::View<ValueType**> destView("destView", numTeams, numCols);
  // make a host copy of the dest view that we can check below
  // to be all zeros
  auto destViewBeforeOp_h = create_host_space_copy(destView);

  // each team stores the distance of the returned iterator from the
  // beginning of the interval that team operates on and then we check
  // that these distances match the expectation
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(sourceView1, sourceView2, destView, distancesView,
                   intraTeamSentinelView, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto destViewAfterOp_h       = create_host_space_copy(destView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  for (std::size_t i = 0; i < destViewAfterOp_h.extent(0); ++i) {
    for (std::size_t j = 0; j < destViewAfterOp_h.extent(1); ++j) {
      // elements in dest view should be the sum of source elements
      ASSERT_DOUBLE_EQ(destViewAfterOp_h(i, j),
                       cloneOfSourceView1BeforeOp_h(i, j) +
                           cloneOfSourceView2BeforeOp_h(i, j));
    }

    // each team should return an iterator whose distance from the
    // beginning of the row equals the num of columns since
    // each team transforms all elements in each row
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

TEST(std_algorithms_transform_team_test, test_binary_op) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamTransformBinaryOp
}  // namespace stdalgos
}  // namespace Test
