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
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace TeamPartitionCopy {

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

template <class ValueType>
struct GreaterThanValueFunctor {
  ValueType m_val;

  KOKKOS_INLINE_FUNCTION
  GreaterThanValueFunctor(ValueType val) : m_val(val) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(ValueType val) const { return (val > m_val); }
};

template <class SourceViewType, class DestViewType, class DistancesViewType,
          class IntraTeamSentinelView, class ValueType>
struct TestFunctorA {
  SourceViewType m_sourceView;

  DestViewType m_destTrueView;
  DestViewType m_destFalseView;

  DistancesViewType m_distancesTrueView;
  DistancesViewType m_distancesFalseView;
  IntraTeamSentinelView m_intraTeamSentinelView;

  ValueType m_threshold;
  int m_apiPick;

  TestFunctorA(const SourceViewType sourceView, const DestViewType destTrueView,
               const DestViewType destFalseView,
               const DistancesViewType distancesTrueView,
               const DistancesViewType distancesFalseView,
               const IntraTeamSentinelView intraTeamSentinelView,
               ValueType threshold, int apiPick)
      : m_sourceView(sourceView),
        m_destTrueView(destTrueView),
        m_destFalseView(destFalseView),
        m_distancesTrueView(distancesTrueView),
        m_distancesFalseView(distancesFalseView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_threshold(threshold),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowViewFrom =
        Kokkos::subview(m_sourceView, myRowIndex, Kokkos::ALL());

    auto myRowViewDestTrue =
        Kokkos::subview(m_destTrueView, myRowIndex, Kokkos::ALL());
    auto myRowViewDestFalse =
        Kokkos::subview(m_destFalseView, myRowIndex, Kokkos::ALL());

    ptrdiff_t resultDist1 = 0;
    ptrdiff_t resultDist2 = 0;

    GreaterThanValueFunctor predicate(m_threshold);
    if (m_apiPick == 0) {
      const auto result = KE::partition_copy(
          member, KE::cbegin(myRowViewFrom), KE::cend(myRowViewFrom),
          KE::begin(myRowViewDestTrue), KE::begin(myRowViewDestFalse),
          predicate);
      resultDist1 = KE::distance(KE::begin(myRowViewDestTrue), result.first);
      resultDist2 = KE::distance(KE::begin(myRowViewDestFalse), result.second);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesTrueView(myRowIndex)  = resultDist1;
        m_distancesFalseView(myRowIndex) = resultDist2;
      });
    }

    else if (m_apiPick == 1) {
      const auto result =
          KE::partition_copy(member, myRowViewFrom, myRowViewDestTrue,
                             myRowViewDestFalse, predicate);
      resultDist1 = KE::distance(KE::begin(myRowViewDestTrue), result.first);
      resultDist2 = KE::distance(KE::begin(myRowViewDestFalse), result.second);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesTrueView(myRowIndex)  = resultDist1;
        m_distancesFalseView(myRowIndex) = resultDist2;
      });
    }

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck1 = team_members_have_matching_result(
        member, resultDist1, m_distancesTrueView(myRowIndex));
    const bool intraTeamCheck2 = team_members_have_matching_result(
        member, resultDist2, m_distancesFalseView(myRowIndex));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(myRowIndex) = intraTeamCheck1 && intraTeamCheck2;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId,
            const std::string& sIn) {
  /* description:
     use a rank-2 view randomly filled with values in a range (a,b)
     and run a team-level partition_copy with predicate = IsGreaterThanValue
     where threshold is set to a number larger than b above
   */
  const auto threshold           = static_cast<ValueType>(1103);
  const auto valueForSureGreater = static_cast<ValueType>(2103);
  const auto valueForSureSmaller = static_cast<ValueType>(111);

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // construct in memory space associated with default exespace
  auto sourceView =
      create_view<ValueType>(LayoutTag{}, numTeams, numCols, "sourceView");

  // sourceView might not deep copyable (e.g. strided layout) so to
  // randomize it, we make a new view that is for sure deep copyable,
  // modify it on the host, deep copy to device and then launch
  // a kernel to copy to sourceView
  auto sourceView_dc =
      create_deep_copyable_compatible_view_with_same_extent(sourceView);
  auto sourceView_dc_h = create_mirror_view(Kokkos::HostSpace(), sourceView_dc);

  if (sIn == "trivialEmpty") {
    // do nothing
  }

  else if (sIn == "allTrue") {
    // randomly fill with values greater than threshold
    // so that all elements in each row satisfy the predicate
    // so this counts as being partitioned
    Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace> pool(
        452377);
    Kokkos::fill_random(sourceView_dc_h, pool, ValueType(2001),
                        ValueType(2501));
  }

  else if (sIn == "allFalse") {
    // randomly fill the view with values smaller than threshold
    // and even in this case each row counts as partitioned
    Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace> pool(
        452377);
    Kokkos::fill_random(sourceView_dc_h, pool, ValueType(0), ValueType(101));
  }

  else if (sIn == "random") {
    // randomly select a location and make all values before that
    // larger than threshol and all values after to be smaller than threshold
    // so that this picked location does partition the range
    UnifDist<int> indexProducer(0, numCols - 1, 3432779);
    for (std::size_t i = 0; i < sourceView_dc_h.extent(0); ++i) {
      const std::size_t a = indexProducer();
      for (std::size_t j = 0; j < a; ++j) {
        sourceView_dc_h(i, j) = valueForSureGreater;
      }
      for (std::size_t j = a; j < numCols; ++j) {
        sourceView_dc_h(i, j) = valueForSureSmaller;
      }
    }
  }

  // copy to sourceView_dc and then to sourceView
  Kokkos::deep_copy(sourceView_dc, sourceView_dc_h);
  // use CTAD
  CopyFunctorRank2 F1(sourceView_dc, sourceView);
  Kokkos::parallel_for("copy", sourceView.extent(0) * sourceView.extent(1), F1);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // create the destination views
  Kokkos::View<ValueType**> destTrueView("destViewTrue", numTeams, numCols);
  Kokkos::View<ValueType**> destFalseView("destViewFalse", numTeams, numCols);

  // to verify that things work, each team stores the result
  // and then we check that these match what we expect
  Kokkos::View<std::size_t*> distancesTrueView("distancesTrue", numTeams);
  Kokkos::View<std::size_t*> distancesFalseView("distancesFalse", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(sourceView, destTrueView, destFalseView, distancesTrueView,
                   distancesFalseView, intraTeamSentinelView, threshold, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  auto distancesTrueView_h     = create_host_space_copy(distancesTrueView);
  auto distancesFalseView_h    = create_host_space_copy(distancesFalseView);
  auto sourceViewAfterOp_h     = create_host_space_copy(sourceView);
  auto destTrueViewAfterOp_h   = create_host_space_copy(destTrueView);
  auto destFalseViewAfterOp_h  = create_host_space_copy(destFalseView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);

  Kokkos::View<ValueType**, Kokkos::HostSpace> stdDestTrueView(
      "stdDestTrueView", numTeams, numCols);
  Kokkos::View<ValueType**, Kokkos::HostSpace> stdDestFalseView(
      "stdDestFalseView", numTeams, numCols);
  GreaterThanValueFunctor predicate(threshold);

  for (std::size_t i = 0; i < sourceView_dc_h.extent(0); ++i) {
    auto myRowSource    = Kokkos::subview(sourceView_dc_h, i, Kokkos::ALL());
    auto myRowDestTrue  = Kokkos::subview(stdDestTrueView, i, Kokkos::ALL());
    auto myRowDestFalse = Kokkos::subview(stdDestFalseView, i, Kokkos::ALL());

    const auto stdResult = std::partition_copy(
        KE::cbegin(myRowSource), KE::cend(myRowSource),
        KE::begin(myRowDestTrue), KE::begin(myRowDestFalse), predicate);
    // our result must match std
    const std::size_t stdDistanceTrue =
        KE::distance(KE::begin(myRowDestTrue), stdResult.first);
    const std::size_t stdDistanceFalse =
        KE::distance(KE::begin(myRowDestFalse), stdResult.second);
    ASSERT_EQ(stdDistanceTrue, distancesTrueView_h(i));
    ASSERT_EQ(stdDistanceFalse, distancesFalseView_h(i));
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }

  expect_equal_host_views(sourceView_dc_h, sourceViewAfterOp_h);
  expect_equal_host_views(destTrueViewAfterOp_h, stdDestTrueView);
  expect_equal_host_views(destFalseViewAfterOp_h, stdDestFalseView);
}

template <class LayoutTag, class ValueType>
void run_all_scenarios(const std::string& name, const std::vector<int>& cols) {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : cols) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId, name);
      }
    }
  }
}

TEST(std_algorithms_partition_copy_team_test, empty) {
  const std::string name      = "trivialEmpty";
  const std::vector<int> cols = {0};
  run_all_scenarios<DynamicTag, double>(name, cols);
  run_all_scenarios<StridedTwoRowsTag, double>(name, cols);
  run_all_scenarios<StridedThreeRowsTag, int>(name, cols);
}

TEST(std_algorithms_partition_copy_team_test, all_true) {
  const std::string name      = "allTrue";
  const std::vector<int> cols = {13, 101, 1444, 5153};
  run_all_scenarios<DynamicTag, double>(name, cols);
  run_all_scenarios<StridedTwoRowsTag, double>(name, cols);
  run_all_scenarios<StridedThreeRowsTag, int>(name, cols);
}

TEST(std_algorithms_partition_copy_team_test, all_false) {
  const std::string name      = "allFalse";
  const std::vector<int> cols = {13, 101, 1444, 5153};
  run_all_scenarios<DynamicTag, double>(name, cols);
  run_all_scenarios<StridedTwoRowsTag, double>(name, cols);
  run_all_scenarios<StridedThreeRowsTag, int>(name, cols);
}

TEST(std_algorithms_partition_copy_team_test, random) {
  const std::string name      = "random";
  const std::vector<int> cols = {13, 101, 1444, 5153};
  run_all_scenarios<DynamicTag, double>(name, cols);
  run_all_scenarios<StridedTwoRowsTag, double>(name, cols);
  run_all_scenarios<StridedThreeRowsTag, int>(name, cols);
}

}  // namespace TeamPartitionCopy
}  // namespace stdalgos
}  // namespace Test
