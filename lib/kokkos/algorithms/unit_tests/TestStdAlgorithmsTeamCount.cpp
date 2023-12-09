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
namespace TeamCount {

namespace KE = Kokkos::Experimental;

template <class ViewType, class ValuesViewType, class CountsViewType,
          class IntraTeamSentinelView>
struct TestFunctorA {
  ViewType m_view;
  ValuesViewType m_valuesView;
  CountsViewType m_countsView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const ViewType view, const ValuesViewType valuesView,
               const CountsViewType countsView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_view(view),
        m_valuesView(valuesView),
        m_countsView(countsView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto rowIndex = member.league_rank();
    const auto value    = m_valuesView(rowIndex);
    auto rowView        = Kokkos::subview(m_view, rowIndex, Kokkos::ALL());
    std::size_t result  = 0;

    switch (m_apiPick) {
      case 0: {
        result =
            KE::count(member, KE::cbegin(rowView), KE::cend(rowView), value);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_countsView(rowIndex) = result; });

        break;
      }

      case 1: {
        result = KE::count(member, rowView, value);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_countsView(rowIndex) = result; });

        break;
      }
    }

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck = team_members_have_matching_result(
        member, result, m_countsView(rowIndex));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(rowIndex) = intraTeamCheck;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(const bool searched_value_exist, std::size_t numTeams,
            std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level count
   */

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // create a view in the memory space associated with default exespace
  // with as many rows as the number of teams and fill it with random
  // values from an arbitrary range.

  // Boundaries choosen so that every drawn number is at least once in the given
  // row
  const ValueType lowerBound = numCols / 4;
  const ValueType upperBound = 1 + numCols * 3 / 4;
  const auto bounds          = make_bounds(lowerBound, upperBound);

  auto [dataView, dataViewBeforeOp_h] = create_random_view_and_host_clone(
      LayoutTag{}, numTeams, numCols, bounds, "dataView");

  // If searched_value_exist == true, we want to ensure that count result is >
  // 0, so we randomly pick a value to look for from a given row.
  //
  // If searched_value_exist == false, we want to ensure that count returns 0,
  // so we pick a value that's outside of view boundaries.
  Kokkos::View<ValueType*> valuesView("valuesView", numTeams);
  auto valuesView_h = create_mirror_view(Kokkos::HostSpace(), valuesView);

  using rand_pool =
      Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>;
  rand_pool pool(lowerBound * upperBound);

  if (searched_value_exist) {
    Kokkos::View<std::size_t*, Kokkos::DefaultHostExecutionSpace> randomIndices(
        "randomIndices", numTeams);
    Kokkos::fill_random(randomIndices, pool, 0, numCols);

    for (std::size_t i = 0; i < numTeams; ++i) {
      const std::size_t j = randomIndices(i);
      valuesView_h(i)     = dataViewBeforeOp_h(i, j);
    }
  } else {
    Kokkos::fill_random(valuesView_h, pool, upperBound, upperBound * 2);
  }

  Kokkos::deep_copy(valuesView, valuesView_h);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // to verify that things work, each team stores the result of its count
  // call, and then we check that these match what we expect
  Kokkos::View<std::size_t*> countsView("countsView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(dataView, valuesView, countsView, intraTeamSentinelView,
                   apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  auto countsView_h            = create_host_space_copy(countsView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  if (searched_value_exist) {
    for (std::size_t i = 0; i < dataView.extent(0); ++i) {
      auto rowFrom = Kokkos::subview(dataViewBeforeOp_h, i, Kokkos::ALL());
      const auto rowFromBegin = KE::cbegin(rowFrom);
      const auto rowFromEnd   = KE::cend(rowFrom);
      const auto val          = valuesView_h(i);

      const std::size_t result = std::count(rowFromBegin, rowFromEnd, val);
      ASSERT_EQ(result, countsView_h(i));
      ASSERT_TRUE(intraTeamSentinelView_h(i));
    }
  } else {
    for (std::size_t i = 0; i < countsView.extent(0); ++i) {
      constexpr std::size_t zero = 0;
      ASSERT_EQ(countsView_h(i), zero);
      ASSERT_TRUE(intraTeamSentinelView_h(i));
    }
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios(const bool searchedValueExist) {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(searchedValueExist, numTeams, numCols,
                                     apiId);
      }
    }
  }
}

TEST(std_algorithms_count_team_test, count_returns_nonzero) {
  constexpr bool searchedValueExist = true;
  run_all_scenarios<DynamicTag, double>(searchedValueExist);
  run_all_scenarios<StridedTwoRowsTag, int>(searchedValueExist);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(searchedValueExist);
}

TEST(std_algorithms_count_team_test, count_returns_zero) {
  constexpr bool searchedValueExist = false;
  run_all_scenarios<DynamicTag, double>(searchedValueExist);
  run_all_scenarios<StridedTwoRowsTag, int>(searchedValueExist);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(searchedValueExist);
}

}  // namespace TeamCount
}  // namespace stdalgos
}  // namespace Test
