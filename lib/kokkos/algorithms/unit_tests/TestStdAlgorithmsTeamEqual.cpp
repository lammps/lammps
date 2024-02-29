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
namespace TeamEqual {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct EqualFunctor {
  KOKKOS_INLINE_FUNCTION bool operator()(const ValueType& lhs,
                                         const ValueType& rhs) const {
    return lhs == rhs;
  }
};

template <class DataViewType, class CompViewType, class ResultsViewType,
          class IntraTeamSentinelView, class BinaryPredType>
struct TestFunctorA {
  DataViewType m_dataView;
  CompViewType m_compView;
  ResultsViewType m_resultsView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;
  BinaryPredType m_binaryPred;

  TestFunctorA(const DataViewType dataView, const CompViewType compView,
               const ResultsViewType resultsView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick,
               BinaryPredType binaryPred)
      : m_dataView(dataView),
        m_compView(compView),
        m_resultsView(resultsView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick),
        m_binaryPred(binaryPred) {}

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
        result = KE::equal(member, dataBegin, dataEnd, compBegin);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 1: {
        result = KE::equal(member, dataBegin, dataEnd, compBegin, m_binaryPred);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 2: {
        result = KE::equal(member, rowData, rowComp);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 3: {
        result = KE::equal(member, rowData, rowComp, m_binaryPred);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });

        break;
      }

      case 4: {
        result = KE::equal(member, dataBegin, dataEnd, compBegin, compEnd);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 5: {
        result = KE::equal(member, dataBegin, dataEnd, compBegin, compEnd,
                           m_binaryPred);
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
void test_A(const bool viewsAreEqual, std::size_t numTeams, std::size_t numCols,
            int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level equal
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

  // create the view to store results of equal()
  Kokkos::View<bool*> resultsView("resultsView", numTeams);
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
    auto rowData = Kokkos::subview(dataViewBeforeOp_h, i, Kokkos::ALL());
    const auto dataBegin = KE::cbegin(rowData);
    const auto dataEnd   = KE::cend(rowData);

    auto rowComp         = Kokkos::subview(compView_h, i, Kokkos::ALL());
    const auto compBegin = KE::cbegin(rowComp);
    const auto compEnd   = KE::cend(rowComp);

    ASSERT_TRUE(intraTeamSentinelView_h(i));
    switch (apiId) {
      case 0:
      case 2: {
        const bool result = std::equal(dataBegin, dataEnd, compBegin);

        if (viewsAreEqual) {
          EXPECT_TRUE(resultsView_h(i));
        } else {
          ASSERT_EQ(result, resultsView_h(i));
        }

        break;
      }

      case 1:
      case 3: {
        const bool result =
            std::equal(dataBegin, dataEnd, compBegin, binaryPred);

        if (viewsAreEqual) {
          EXPECT_TRUE(resultsView_h(i));
        } else {
          ASSERT_EQ(result, resultsView_h(i));
        }

        break;
      }

      case 4: {
        const bool result = std::equal(dataBegin, dataEnd, compBegin, compEnd);

        if (viewsAreEqual) {
          EXPECT_TRUE(resultsView_h(i));
        } else {
          ASSERT_EQ(result, resultsView_h(i));
        }

        break;
      }

      case 5: {
        const bool result =
            std::equal(dataBegin, dataEnd, compBegin, compEnd, binaryPred);

        if (viewsAreEqual) {
          EXPECT_TRUE(resultsView_h(i));
        } else {
          ASSERT_EQ(result, resultsView_h(i));
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
      for (int apiId : {0, 1, 2, 3, 4, 5}) {
        test_A<LayoutTag, ValueType>(viewsAreEqual, numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_equal_team_test, views_are_equal) {
  constexpr bool viewsAreEqual = true;
  run_all_scenarios<DynamicTag, double>(viewsAreEqual);
  run_all_scenarios<StridedTwoRowsTag, int>(viewsAreEqual);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(viewsAreEqual);
}

TEST(std_algorithms_equal_team_test, views_are_not_equal) {
  constexpr bool viewsAreEqual = false;
  run_all_scenarios<DynamicTag, double>(viewsAreEqual);
  run_all_scenarios<StridedTwoRowsTag, int>(viewsAreEqual);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(viewsAreEqual);
}

}  // namespace TeamEqual
}  // namespace stdalgos
}  // namespace Test
