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
namespace TeamMismatch {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct EqualFunctor {
  KOKKOS_INLINE_FUNCTION bool operator()(const ValueType& lhs,
                                         const ValueType& rhs) const {
    return lhs == rhs;
  }
};

template <class DataViewType, class CompViewType, class ResultsViewType,
          class IntraTeamSentinelView, class BinaryOpType>
struct TestFunctorA {
  DataViewType m_dataView;
  CompViewType m_compView;
  ResultsViewType m_resultsView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;
  BinaryOpType m_binaryOp;

  TestFunctorA(const DataViewType dataView, const CompViewType compView,
               const ResultsViewType resultsView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick,
               BinaryOpType binaryOp)
      : m_dataView(dataView),
        m_compView(compView),
        m_resultsView(resultsView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick),
        m_binaryOp(binaryOp) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto rowIndex = member.league_rank();

    auto rowData   = Kokkos::subview(m_dataView, rowIndex, Kokkos::ALL());
    auto dataBegin = KE::begin(rowData);
    auto dataEnd   = KE::end(rowData);

    auto rowComp   = Kokkos::subview(m_compView, rowIndex, Kokkos::ALL());
    auto compBegin = KE::begin(rowComp);
    auto compEnd   = KE::end(rowComp);

    ptrdiff_t dataDist = 0;
    ptrdiff_t compDist = 0;

    switch (m_apiPick) {
      case 0: {
        auto [dataIt, compIt] =
            KE::mismatch(member, dataBegin, dataEnd, compBegin, compEnd);

        dataDist = KE::distance(dataBegin, dataIt);
        compDist = KE::distance(compBegin, compIt);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_resultsView(rowIndex) = Kokkos::make_pair(dataDist, compDist);
        });

        break;
      }

      case 1: {
        const auto [dataIt, compIt] = KE::mismatch(
            member, dataBegin, dataEnd, compBegin, compEnd, m_binaryOp);

        dataDist = KE::distance(dataBegin, dataIt);
        compDist = KE::distance(compBegin, compIt);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_resultsView(rowIndex) = Kokkos::make_pair(dataDist, compDist);
        });

        break;
      }

      case 2: {
        const auto [dataIt, compIt] = KE::mismatch(member, rowData, rowComp);

        dataDist = KE::distance(dataBegin, dataIt);
        compDist = KE::distance(compBegin, compIt);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_resultsView(rowIndex) = Kokkos::make_pair(dataDist, compDist);
        });

        break;
      }

      case 3: {
        const auto [dataIt, compIt] =
            KE::mismatch(member, rowData, rowComp, m_binaryOp);

        dataDist = KE::distance(dataBegin, dataIt);
        compDist = KE::distance(compBegin, compIt);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_resultsView(rowIndex) = Kokkos::make_pair(dataDist, compDist);
        });

        break;
      }
    }

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck1 = team_members_have_matching_result(
        member, dataDist, m_resultsView(rowIndex).first);
    const bool intraTeamCheck2 = team_members_have_matching_result(
        member, compDist, m_resultsView(rowIndex).second);
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(rowIndex) = intraTeamCheck1 && intraTeamCheck2;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(const bool viewsAreEqual, std::size_t numTeams, std::size_t numCols,
            int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level mismatch
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

  // create a view to compare it with dataView. If viewsAreEqual == true,
  // compView is a copy of dataView. If viewsAreEqual == false, compView is
  // randomly filled
  auto compView   = create_deep_copyable_compatible_clone(dataView);
  auto compView_h = create_mirror_view(Kokkos::HostSpace(), compView);
  if (viewsAreEqual) {
    Kokkos::deep_copy(compView_h, dataViewBeforeOp_h);
  } else {
    using rand_pool =
        Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>;
    rand_pool pool(lowerBound * upperBound);
    Kokkos::fill_random(compView_h, pool, lowerBound, upperBound);
  }

  Kokkos::deep_copy(compView, compView_h);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // create the view to store results of mismatch()
  Kokkos::View<Kokkos::pair<std::size_t, std::size_t>*> resultsView(
      "resultsView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  EqualFunctor<ValueType> binaryPred{};

  // use CTAD for functor
  TestFunctorA fnc(dataView, compView, resultsView, intraTeamSentinelView,
                   apiId, binaryPred);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel and check
  // -----------------------------------------------
  auto resultsView_h           = create_host_space_copy(resultsView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);

  for (std::size_t i = 0; i < dataView.extent(0); ++i) {
    ASSERT_TRUE(intraTeamSentinelView_h(i));

    auto rowData = Kokkos::subview(dataViewBeforeOp_h, i, Kokkos::ALL());

    const auto dataBegin = KE::cbegin(rowData);
    const auto dataEnd   = KE::cend(rowData);

    const std::size_t dataBeginEndDist = KE::distance(dataBegin, dataEnd);

    auto rowComp = Kokkos::subview(compView_h, i, Kokkos::ALL());

    const auto compBegin = KE::cbegin(rowComp);
    const auto compEnd   = KE::cend(rowComp);

    const std::size_t compBeginEndDist = KE::distance(compBegin, compEnd);

    switch (apiId) {
      case 0:
      case 2: {
        const auto [dataIt, compIt] =
            std::mismatch(dataBegin, dataEnd, compBegin, compEnd);

        const std::size_t dataDist = KE::distance(dataBegin, dataIt);
        const std::size_t compDist = KE::distance(compBegin, compIt);

        if (viewsAreEqual) {
          ASSERT_EQ(dataBeginEndDist, resultsView_h(i).first);
          ASSERT_EQ(compBeginEndDist, resultsView_h(i).second);
        } else {
          ASSERT_EQ(dataDist, resultsView_h(i).first);
          ASSERT_EQ(compDist, resultsView_h(i).second);
        }

        break;
      }

      case 1:
      case 3: {
        const auto [dataIt, compIt] =
            std::mismatch(dataBegin, dataEnd, compBegin, compEnd, binaryPred);

        const std::size_t dataDist = KE::distance(dataBegin, dataIt);
        const std::size_t compDist = KE::distance(compBegin, compIt);

        if (viewsAreEqual) {
          ASSERT_EQ(dataBeginEndDist, resultsView_h(i).first);
          ASSERT_EQ(compBeginEndDist, resultsView_h(i).second);
        } else {
          ASSERT_EQ(dataDist, resultsView_h(i).first);
          ASSERT_EQ(compDist, resultsView_h(i).second);
        }

        break;
      }
    }
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios(const bool viewsAreEqual) {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1, 2, 3}) {
        test_A<LayoutTag, ValueType>(viewsAreEqual, numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_mismatch_team_test, views_are_equal) {
  constexpr bool viewsAreEqual = true;
  run_all_scenarios<DynamicTag, double>(viewsAreEqual);
  run_all_scenarios<StridedTwoRowsTag, int>(viewsAreEqual);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(viewsAreEqual);
}

TEST(std_algorithms_mismatch_team_test, views_are_not_equal) {
  constexpr bool viewsAreEqual = false;
  run_all_scenarios<DynamicTag, double>(viewsAreEqual);
  run_all_scenarios<StridedTwoRowsTag, int>(viewsAreEqual);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(viewsAreEqual);
}

}  // namespace TeamMismatch
}  // namespace stdalgos
}  // namespace Test
