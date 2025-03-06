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
namespace TeamLexicographicalCompare {

namespace KE = Kokkos::Experimental;

enum class TestCaseType { ViewsAreEqual, FirstIsLess, FirstIsGreater };

template <class ValueType>
struct LessFunctor {
  KOKKOS_INLINE_FUNCTION bool operator()(const ValueType& lhs,
                                         const ValueType& rhs) const {
    return lhs < rhs;
  }
};

template <class DataViewType, class CompViewType, class ResultsViewType,
          class IntraTeamSentinelView, class BinaryCompType>
struct TestFunctorA {
  DataViewType m_dataView;
  CompViewType m_compView;
  ResultsViewType m_resultsView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;
  BinaryCompType m_binaryComp;

  TestFunctorA(const DataViewType dataView, const CompViewType compView,
               const ResultsViewType resultsView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick,
               BinaryCompType binaryComp)
      : m_dataView(dataView),
        m_compView(compView),
        m_resultsView(resultsView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick),
        m_binaryComp(binaryComp) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto rowIndex = member.league_rank();

    auto rowData         = Kokkos::subview(m_dataView, rowIndex, Kokkos::ALL());
    const auto dataBegin = KE::cbegin(rowData);
    const auto dataEnd   = KE::cend(rowData);

    auto rowComp         = Kokkos::subview(m_compView, rowIndex, Kokkos::ALL());
    const auto compBegin = KE::cbegin(rowComp);
    const auto compEnd   = KE::cend(rowComp);

    bool result = false;
    switch (m_apiPick) {
      case 0: {
        result = KE::lexicographical_compare(member, dataBegin, dataEnd,
                                             compBegin, compEnd);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 1: {
        result = KE::lexicographical_compare(member, rowData, rowComp);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 2: {
        result = KE::lexicographical_compare(member, dataBegin, dataEnd,
                                             compBegin, compEnd, m_binaryComp);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 3: {
        result =
            KE::lexicographical_compare(member, rowData, rowComp, m_binaryComp);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }
    }

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck = team_members_have_matching_result(
        member, result, m_resultsView(rowIndex));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(rowIndex) = intraTeamCheck;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(const TestCaseType testCase, std::size_t numTeams,
            std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level lexicographical_compare
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

  // create a view to compare it with dataView. If testCase == ViewsAreEqual,
  // compView is a copy of dataView. If testCase == FirstIsLess, we want the
  // dataView to be lexicographically less (and compView - greater). If testCase
  // == FirstIsGreater, we want the dataView to be lexicographically greater
  // (and compView - less).
  auto compEqualView   = create_deep_copyable_compatible_clone(dataView);
  auto compEqualView_h = create_mirror_view(Kokkos::HostSpace(), compEqualView);
  Kokkos::deep_copy(compEqualView_h, dataViewBeforeOp_h);
  const auto middle = numCols / 2;
  switch (testCase) {
    case TestCaseType::ViewsAreEqual: {
      // Do nothing - deep_copy was already done
      break;
    }

    case TestCaseType::FirstIsLess: {
      for (std::size_t i = 0; i < compEqualView_h.extent(0); ++i) {
        compEqualView_h(i, middle) += 1;
      }

      break;
    }

    case TestCaseType::FirstIsGreater: {
      for (std::size_t i = 0; i < compEqualView_h.extent(0); ++i) {
        compEqualView_h(i, middle) -= 1;
      }

      break;
    }
  }

  Kokkos::deep_copy(compEqualView, compEqualView_h);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // create the view to store results of equal()
  Kokkos::View<bool*> resultsView("resultsView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  LessFunctor<ValueType> binaryComp{};

  // use CTAD for functor
  TestFunctorA fnc(dataView, compEqualView, resultsView, intraTeamSentinelView,
                   apiId, binaryComp);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel and check
  // -----------------------------------------------
  auto resultsView_h           = create_host_space_copy(resultsView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);

  for (std::size_t i = 0; i < dataView.extent(0); ++i) {
    auto rowData = Kokkos::subview(dataViewBeforeOp_h, i, Kokkos::ALL());
    const auto dataBegin = KE::cbegin(rowData);
    const auto dataEnd   = KE::cend(rowData);

    auto rowComp         = Kokkos::subview(compEqualView_h, i, Kokkos::ALL());
    const auto compBegin = KE::cbegin(rowComp);
    const auto compEnd   = KE::cend(rowComp);

    ASSERT_TRUE(intraTeamSentinelView_h(i));
    switch (apiId) {
      case 0:
      case 1: {
        const bool result = std::lexicographical_compare(dataBegin, dataEnd,
                                                         compBegin, compEnd);

        switch (testCase) {
          case TestCaseType::ViewsAreEqual:
          case TestCaseType::FirstIsGreater: {
            EXPECT_FALSE(resultsView_h(i));
            ASSERT_EQ(result, resultsView_h(i));
            break;
          }

          case TestCaseType::FirstIsLess: {
            EXPECT_TRUE(resultsView_h(i));
            ASSERT_EQ(result, resultsView_h(i));
            break;
          }
        }

        break;
      }

      case 2:
      case 3: {
        const bool result = std::lexicographical_compare(
            dataBegin, dataEnd, compBegin, compEnd, binaryComp);

        switch (testCase) {
          case TestCaseType::ViewsAreEqual:
          case TestCaseType::FirstIsGreater: {
            EXPECT_FALSE(resultsView_h(i));
            ASSERT_EQ(result, resultsView_h(i));
            break;
          }

          case TestCaseType::FirstIsLess: {
            EXPECT_TRUE(resultsView_h(i));
            ASSERT_EQ(result, resultsView_h(i));
            break;
          }
        }

        break;
      }
    }
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios(const TestCaseType testCase) {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1, 2, 3}) {
        test_A<LayoutTag, ValueType>(testCase, numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_lexicographical_compare_team_test, views_are_equal) {
  constexpr TestCaseType testCaseType = TestCaseType::ViewsAreEqual;
  run_all_scenarios<DynamicTag, double>(testCaseType);
  run_all_scenarios<StridedTwoRowsTag, int>(testCaseType);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(testCaseType);
}

TEST(std_algorithms_lexicographical_compare_team_test, first_view_is_less) {
  constexpr TestCaseType testCaseType = TestCaseType::FirstIsLess;
  run_all_scenarios<DynamicTag, double>(testCaseType);
  run_all_scenarios<StridedTwoRowsTag, int>(testCaseType);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(testCaseType);
}

TEST(std_algorithms_lexicographical_compare_team_test, first_view_is_greater) {
  constexpr TestCaseType testCaseType = TestCaseType::FirstIsGreater;
  run_all_scenarios<DynamicTag, double>(testCaseType);
  run_all_scenarios<StridedTwoRowsTag, int>(testCaseType);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(testCaseType);
}

}  // namespace TeamLexicographicalCompare
}  // namespace stdalgos
}  // namespace Test
