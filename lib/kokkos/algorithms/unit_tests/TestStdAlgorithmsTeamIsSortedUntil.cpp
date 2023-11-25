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
namespace TeamIsSortedUntil {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<int> {
  using dist_type = std::uniform_int_distribution<int>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist(int a, int b, std::size_t seedIn) : m_dist(a, b) {
    m_gen.seed(seedIn);
  }

  int operator()() { return m_dist(m_gen); }
};

template <class ViewType, class DistancesViewType, class IntraTeamSentinelView>
struct TestFunctorA {
  ViewType m_view;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const ViewType view, const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_view(view),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto leagueRank = member.league_rank();
    const auto myRowIndex = leagueRank;
    auto myRowView        = Kokkos::subview(m_view, myRowIndex, Kokkos::ALL());
    ptrdiff_t resultDist  = 0;

    if (m_apiPick == 0) {
      auto it    = KE::is_sorted_until(member, KE::cbegin(myRowView),
                                    KE::cend(myRowView));
      resultDist = KE::distance(KE::cbegin(myRowView), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    } else if (m_apiPick == 1) {
      auto it    = KE::is_sorted_until(member, myRowView);
      resultDist = KE::distance(KE::begin(myRowView), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    }
#if not defined KOKKOS_ENABLE_OPENMPTARGET
    else if (m_apiPick == 2) {
      using value_type = typename ViewType::value_type;
      auto it          = KE::is_sorted_until(member, KE::cbegin(myRowView),
                                    KE::cend(myRowView),
                                    CustomLessThanComparator<value_type>{});
      resultDist       = KE::distance(KE::cbegin(myRowView), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    }

    else if (m_apiPick == 3) {
      using value_type = typename ViewType::value_type;
      auto it          = KE::is_sorted_until(member, myRowView,
                                    CustomLessThanComparator<value_type>{});
      resultDist       = KE::distance(KE::begin(myRowView), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    }
#endif

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
void test_A(std::size_t numTeams, std::size_t numCols, int apiId,
            const std::string& sIn) {
  /* description:
     use a rank-2 view and run a team-level is_sorted_until
     for various trivial and non trivial scenarios
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

  if (sIn == "trivialEmpty" || sIn == "trivialOneElement") {
    // do not do anything
  }

  else if (sIn == "nontrivialUntilLast") {
    for (std::size_t i = 0; i < dataView_dc_h.extent(0); ++i) {
      for (std::size_t j = 0; j < dataView_dc_h.extent(1); ++j) {
        dataView_dc_h(i, j) = ValueType(j);
      }
    }
  }

  else if (sIn == "nontrivialRandom") {
    // randomly fill the view
    Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace> pool(
        452377);
    // make the range tight so that we have low likelihood of having
    // sorted data by chance for small num of cols
    Kokkos::fill_random(dataView_dc_h, pool, ValueType(113), ValueType(120));

    /* pick randomly a location 0 < a < numCols-1/2
       annd fill data so that:
       - from 0 to a: is sorted
       - from a to a+3: not sorted
       - from a+3 to numCols-1: is sorted
       this allows us to exercise that the algorithm returns
       the larest sorted interval starting from 0
    */
    assert(numCols > 10);
    const std::size_t midPoint = numCols / 2;

    UnifDist<int> randPoolA(0, midPoint, 3432779);
    for (std::size_t i = 0; i < dataView_dc_h.extent(0); ++i) {
      const std::size_t a = randPoolA();
      for (std::size_t j = 0; j < a; ++j) {
        dataView_dc_h(i, j) = ValueType(j);
      }
      for (std::size_t j = a + 3; j < numCols; ++j) {
        dataView_dc_h(i, j) = ValueType(j);
      }
    }
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

  // each team stores the distance of the returned iterator from the
  // beginning of the interval that team operates on and then we check
  // that these distances match the std result
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(dataView, distancesView, intraTeamSentinelView, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  for (std::size_t i = 0; i < dataView_dc_h.extent(0); ++i) {
    auto myRow = Kokkos::subview(dataView_dc_h, i, Kokkos::ALL());

    std::size_t stdDistance = 0;
    if (apiId <= 1) {
      auto it     = std::is_sorted_until(KE::cbegin(myRow), KE::cend(myRow));
      stdDistance = KE::distance(KE::cbegin(myRow), it);
    } else {
      auto it     = std::is_sorted_until(KE::cbegin(myRow), KE::cend(myRow),
                                     CustomLessThanComparator<ValueType>{});
      stdDistance = KE::distance(KE::cbegin(myRow), it);
    }
    ASSERT_EQ(stdDistance, distancesView_h(i));
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }

  // dataView should remain unchanged
  auto dataViewAfterOp_h = create_host_space_copy(dataView);
  expect_equal_host_views(dataView_dc_h, dataViewAfterOp_h);
}

template <class LayoutTag, class ValueType>
void run_all_scenarios(const std::string& name, const std::vector<int>& cols) {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : cols) {
#if not defined KOKKOS_ENABLE_OPENMPTARGET
      for (int apiId : {0, 1, 2, 3}) {
#else
      for (int apiId : {0, 1}) {
#endif
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId, name);
      }
    }
  }
}

TEST(std_algorithms_is_sorted_until_team_test, test_trivialA) {
  const std::string name      = "trivialEmpty";
  const std::vector<int> cols = {0};

  run_all_scenarios<DynamicTag, double>(name, cols);
  run_all_scenarios<StridedTwoRowsTag, double>(name, cols);
  run_all_scenarios<StridedThreeRowsTag, int>(name, cols);
}

TEST(std_algorithms_is_sorted_until_team_test, test_trivialB) {
  const std::string name      = "trivialOneElement";
  const std::vector<int> cols = {1};
  run_all_scenarios<DynamicTag, double>(name, cols);
  run_all_scenarios<StridedTwoRowsTag, double>(name, cols);
  run_all_scenarios<StridedThreeRowsTag, int>(name, cols);
}

TEST(std_algorithms_is_sorted_until_team_test, test_nontrivialA) {
  const std::string name      = "nontrivialUntilLast";
  const std::vector<int> cols = {13, 101, 1444, 5153};
  run_all_scenarios<DynamicTag, double>(name, cols);
  run_all_scenarios<StridedTwoRowsTag, double>(name, cols);
  run_all_scenarios<StridedThreeRowsTag, int>(name, cols);
}

TEST(std_algorithms_is_sorted_until_team_test, test_nontrivialB) {
  const std::string name      = "nontrivialRandom";
  const std::vector<int> cols = {13, 101, 1444, 5153};
  run_all_scenarios<DynamicTag, double>(name, cols);
  run_all_scenarios<StridedTwoRowsTag, double>(name, cols);
  run_all_scenarios<StridedThreeRowsTag, int>(name, cols);
}

}  // namespace TeamIsSortedUntil
}  // namespace stdalgos
}  // namespace Test
