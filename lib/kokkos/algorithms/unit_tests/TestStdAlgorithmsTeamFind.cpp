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
namespace TeamFind {

namespace KE = Kokkos::Experimental;

template <class DataViewType, class SearchedValuesViewType,
          class DistancesViewType, class IntraTeamSentinelView>
struct TestFunctorA {
  DataViewType m_dataView;
  SearchedValuesViewType m_searchedValuesView;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const DataViewType dataView,
               const SearchedValuesViewType searchedValuesView,
               const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_dataView(dataView),
        m_searchedValuesView(searchedValuesView),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowViewFrom = Kokkos::subview(m_dataView, myRowIndex, Kokkos::ALL());
    const auto searchedValue = m_searchedValuesView(myRowIndex);
    ptrdiff_t resultDist     = 0;

    switch (m_apiPick) {
      case 0: {
        auto it    = KE::find(member, KE::cbegin(myRowViewFrom),
                           KE::cend(myRowViewFrom), searchedValue);
        resultDist = KE::distance(KE::cbegin(myRowViewFrom), it);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_distancesView(myRowIndex) = resultDist;
        });

        break;
      }

      case 1: {
        auto it    = KE::find(member, myRowViewFrom, searchedValue);
        resultDist = KE::distance(KE::begin(myRowViewFrom), it);
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
void test_A(const bool searchedValuesExist, std::size_t numTeams,
            std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level find
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

  // find returns an iterator so to verify that it is correct
  // each team stores the distance of the returned iterator from the
  // beginning of the interval that team operates on and then we check
  // that these distances match the std result
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // If searchedValuesExist == true we want to ensure that each value we're
  // looking for exists in dataView. To do that, for each numTeams, a random j
  // index from a range [0, numCols) is used to obtain a value from dataView.
  //
  // If searchedValuesExist == false we want to ensure the opposite, so every
  // value is less than a lower bound of dataView.
  Kokkos::View<ValueType*> searchedValuesView("searchValuesView", numTeams);
  auto searchedValuesView_h =
      create_mirror_view(Kokkos::HostSpace(), searchedValuesView);

  using rand_pool =
      Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>;
  rand_pool pool(lowerBound * upperBound);

  if (searchedValuesExist) {
    Kokkos::View<std::size_t*, Kokkos::DefaultHostExecutionSpace> randomIndices(
        "randomIndices", numTeams);

    Kokkos::fill_random(randomIndices, pool, 0, numCols);

    for (std::size_t i = 0; i < numTeams; ++i) {
      const std::size_t j     = randomIndices(i);
      searchedValuesView_h(i) = dataViewBeforeOp_h(i, j);
    }
  } else {
    Kokkos::fill_random(searchedValuesView_h, pool, 0, lowerBound);
  }

  Kokkos::deep_copy(searchedValuesView, searchedValuesView_h);

  // use CTAD for functor
  TestFunctorA fnc(dataView, searchedValuesView, distancesView,
                   intraTeamSentinelView, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel and check
  // -----------------------------------------------
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);

  for (std::size_t i = 0; i < dataView.extent(0); ++i) {
    auto rowFrom = Kokkos::subview(dataViewBeforeOp_h, i, Kokkos::ALL());
    const auto rowFromBegin = KE::cbegin(rowFrom);
    const auto rowFromEnd   = KE::cend(rowFrom);

    auto it = std::find(rowFromBegin, rowFromEnd, searchedValuesView_h(i));

    const std::size_t stdDistance      = KE::distance(rowFromBegin, it);
    const std::size_t beginEndDistance = KE::distance(rowFromBegin, rowFromEnd);

    if (searchedValuesExist) {
      EXPECT_LT(stdDistance, beginEndDistance);
    } else {
      ASSERT_EQ(stdDistance, beginEndDistance);
    }

    ASSERT_EQ(stdDistance, distancesView_h(i));
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios(const bool searchedValuesExist) {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(searchedValuesExist, numTeams, numCols,
                                     apiId);
      }
    }
  }
}

TEST(std_algorithms_find_team_test, searched_values_exist) {
  constexpr bool searchedValuesExist = true;

  run_all_scenarios<DynamicTag, double>(searchedValuesExist);
  run_all_scenarios<StridedTwoRowsTag, int>(searchedValuesExist);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(searchedValuesExist);
}

TEST(std_algorithms_find_team_test, searched_values_do_not_exist) {
  constexpr bool searchedValuesExist = false;

  run_all_scenarios<DynamicTag, double>(searchedValuesExist);
  run_all_scenarios<StridedTwoRowsTag, int>(searchedValuesExist);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(searchedValuesExist);
}

}  // namespace TeamFind
}  // namespace stdalgos
}  // namespace Test
