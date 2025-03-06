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
namespace TeamShiftLeft {

namespace KE = Kokkos::Experimental;

template <class ViewType, class DistancesViewType, class IntraTeamSentinelView>
struct TestFunctorA {
  ViewType m_view;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  std::size_t m_shift;
  int m_apiPick;

  TestFunctorA(const ViewType view, const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView,
               std::size_t shift, int apiPick)
      : m_view(view),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_shift(shift),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowView        = Kokkos::subview(m_view, myRowIndex, Kokkos::ALL());
    ptrdiff_t resultDist  = 0;

    if (m_apiPick == 0) {
      auto it = KE::shift_left(member, KE::begin(myRowView), KE::end(myRowView),
                               m_shift);
      resultDist = KE::distance(KE::begin(myRowView), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    } else if (m_apiPick == 1) {
      auto it    = KE::shift_left(member, myRowView, m_shift);
      resultDist = KE::distance(KE::begin(myRowView), it);
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

// shift_left is only supported starting from C++20,
// so put here a working version of the std algo copied from
// https://github.com/llvm/llvm-project/blob/main/libcxx/include/__algorithm/shift_left.h
template <class ForwardIterator>
ForwardIterator my_std_shift_left(
    ForwardIterator first, ForwardIterator last,
    typename std::iterator_traits<ForwardIterator>::difference_type n) {
  if (n == 0) {
    return last;
  }

  ForwardIterator m = first;
  for (; n > 0; --n) {
    if (m == last) {
      return first;
    }
    ++m;
  }
  return std::move(m, last, first);
}

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, std::size_t shift,
            int apiId) {
  /* description:
     randomly fill a rank-2 view and do a team-level KE::shift_left
     using shift as the shift count.
   */

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // create a view in the memory space associated with default exespace
  // with as many rows as the number of teams and fill it with random
  // values from an arbitrary range
  auto [dataView, cloneOfDataViewBeforeOp_h] =
      create_random_view_and_host_clone(
          LayoutTag{}, numTeams, numCols,
          Kokkos::pair<ValueType, ValueType>{11, 523}, "dataView");

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
  TestFunctorA fnc(dataView, distancesView, intraTeamSentinelView, shift,
                   apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run std algo and check
  // -----------------------------------------------
  // here I can use cloneOfDataViewBeforeOp_h to run std algo on
  // since that contains a valid copy of the data
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  for (std::size_t i = 0; i < cloneOfDataViewBeforeOp_h.extent(0); ++i) {
    auto myRow = Kokkos::subview(cloneOfDataViewBeforeOp_h, i, Kokkos::ALL());
    auto it    = my_std_shift_left(KE::begin(myRow), KE::end(myRow), shift);
    const std::size_t stdDistance = KE::distance(KE::begin(myRow), it);
    ASSERT_EQ(stdDistance, distancesView_h(i));
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }

  auto dataViewAfterOp_h = create_host_space_copy(dataView);
  expect_equal_host_views(cloneOfDataViewBeforeOp_h, dataViewAfterOp_h);
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  // prepare a map where, for a given set of num cols
  // we provide a list of shifts to use for testing
  // key = num of columns
  // value = list of shifts
  // Note that the cornerCase number is here since the shiftLeft algo
  // should work even when the shift given is way larger than the range.
  constexpr std::size_t cornerCase                        = 110111;
  const std::map<int, std::vector<std::size_t>> scenarios = {
      {0, {0, cornerCase}},
      {2, {0, 1, 2, cornerCase}},
      {6, {0, 1, 2, 5, cornerCase}},
      {13, {0, 1, 2, 8, 11, cornerCase}},
      {56, {0, 1, 2, 8, 11, 33, 56, cornerCase}},
      {123, {0, 1, 11, 33, 56, 89, 112, cornerCase}},
      {3145, {0, 1, 11, 33, 56, 89, 112, 5677, cornerCase}}};

  for (int numTeams : teamSizesToTest) {
    for (const auto& scenario : scenarios) {
      const std::size_t numCols = scenario.first;
      for (int copyCount : scenario.second) {
        for (int apiId : {0, 1}) {
          test_A<Tag, ValueType>(numTeams, numCols, copyCount, apiId);
        }
      }
    }
  }
}

TEST(std_algorithms_shift_left_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamShiftLeft
}  // namespace stdalgos
}  // namespace Test
