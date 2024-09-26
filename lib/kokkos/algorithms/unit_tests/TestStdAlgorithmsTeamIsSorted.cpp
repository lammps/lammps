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
namespace TeamIsSorted {

namespace KE = Kokkos::Experimental;

template <class ViewType, class ReturnViewType, class IntraTeamSentinelView>
struct TestFunctorA {
  ViewType m_view;
  ReturnViewType m_returnsView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const ViewType view, const ReturnViewType returnsView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_view(view),
        m_returnsView(returnsView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowView        = Kokkos::subview(m_view, myRowIndex, Kokkos::ALL());
    bool result           = false;

    if (m_apiPick == 0) {
      result =
          KE::is_sorted(member, KE::cbegin(myRowView), KE::cend(myRowView));
      Kokkos::single(Kokkos::PerTeam(member),
                     [=, *this]() { m_returnsView(myRowIndex) = result; });
    } else if (m_apiPick == 1) {
      result = KE::is_sorted(member, myRowView);
      Kokkos::single(Kokkos::PerTeam(member),
                     [=, *this]() { m_returnsView(myRowIndex) = result; });
    }
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    else if (m_apiPick == 2) {
      using value_type = typename ViewType::value_type;
      result = KE::is_sorted(member, KE::cbegin(myRowView), KE::cend(myRowView),
                             CustomLessThanComparator<value_type>{});
      Kokkos::single(Kokkos::PerTeam(member),
                     [=, *this]() { m_returnsView(myRowIndex) = result; });
    } else if (m_apiPick == 3) {
      using value_type = typename ViewType::value_type;
      result           = KE::is_sorted(member, myRowView,
                             CustomLessThanComparator<value_type>{});
      Kokkos::single(Kokkos::PerTeam(member),
                     [=, *this]() { m_returnsView(myRowIndex) = result; });
    }
#endif

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck = team_members_have_matching_result(
        member, result, m_returnsView(myRowIndex));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(myRowIndex) = intraTeamCheck;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId,
            bool makeDataSortedOnPurpose) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level is_sorted
   */

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // construct in memory space associated with default exespace
  auto dataView =
      create_view<ValueType>(LayoutTag{}, numTeams, numCols, "dataView");

  // dataView might not deep copyable (e.g. strided layout) so to
  // randomize it, we make a new view that is for sure deep copyable,
  // modify it on the host, deep copy to device and then launch
  // a kernel to copy to dataView
  auto dataView_dc =
      create_deep_copyable_compatible_view_with_same_extent(dataView);
  auto dataView_dc_h = create_mirror_view(Kokkos::HostSpace(), dataView_dc);

  if (makeDataSortedOnPurpose) {
    for (std::size_t i = 0; i < dataView_dc_h.extent(0); ++i) {
      for (std::size_t j = 0; j < dataView_dc_h.extent(1); ++j) {
        dataView_dc_h(i, j) = ValueType(j);
      }
    }
  } else {
    // randomly fill the view
    using rand_pool =
        Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>;
    rand_pool pool(45234977);
    Kokkos::fill_random(dataView_dc_h, pool, ValueType{5}, ValueType{1545});
  }

  // copy to dataView_dc and then to dataView
  Kokkos::deep_copy(dataView_dc, dataView_dc_h);
  // use CTAD
  CopyFunctorRank2 F1(dataView_dc, dataView);
  Kokkos::parallel_for("copy", dataView.extent(0) * dataView.extent(1), F1);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // to verify that things work, each team stores the result
  // and then we check that these match what we expect
  Kokkos::View<bool*> returnView("returnView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(dataView, returnView, intraTeamSentinelView, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  auto returnView_h            = create_host_space_copy(returnView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  for (std::size_t i = 0; i < dataView_dc_h.extent(0); ++i) {
    auto myRow = Kokkos::subview(dataView_dc_h, i, Kokkos::ALL());

    bool stdResult;
    if (apiId <= 1) {
      stdResult = std::is_sorted(KE::cbegin(myRow), KE::cend(myRow));
    } else {
      stdResult = std::is_sorted(KE::cbegin(myRow), KE::cend(myRow),
                                 CustomLessThanComparator<ValueType>{});
    }

    // our result must match std
    EXPECT_TRUE(stdResult == returnView_h(i));

    // check also since we know in advance when data is really sorted.
    // note that we have to be careful because when we have only
    // 0, 1 columns, then the data is sorted by definition
    // and when we have 2 columns it is very likely it is sorted
    // so only do the following check for large enough cols count
    if (numCols <= 1) {
      EXPECT_TRUE(stdResult == true);
    } else if (numCols > 10) {
      EXPECT_TRUE(stdResult == makeDataSortedOnPurpose);
    }
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }

  // dataView should remain unchanged
  auto dataViewAfterOp_h = create_host_space_copy(dataView);
  expect_equal_host_views(dataView_dc_h, dataViewAfterOp_h);
}

template <class LayoutTag, class ValueType>
void run_all_scenarios(bool makeDataSortedOnPurpose) {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 5153}) {
#ifndef KOKKOS_ENABLE_OPENMPTARGET
      for (int apiId : {0, 1, 2, 3}) {
#else
      for (int apiId : {0, 1}) {
#endif
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId,
                                     makeDataSortedOnPurpose);
      }
    }
  }
}

TEST(std_algorithms_is_sorted_team_test,
     test_data_almost_certainly_not_sorted) {
  run_all_scenarios<DynamicTag, double>(false);
  run_all_scenarios<StridedTwoRowsTag, double>(false);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(false);
}

TEST(std_algorithms_is_sorted_team_test, test_data_certainly_sorted) {
  run_all_scenarios<DynamicTag, double>(true);
  run_all_scenarios<StridedTwoRowsTag, double>(true);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(true);
}

}  // namespace TeamIsSorted
}  // namespace stdalgos
}  // namespace Test
