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
namespace TeamRemoveCopy {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<int> {
  using dist_type = std::uniform_int_distribution<int>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist(int b, std::size_t seedIn) : m_dist(0, b) { m_gen.seed(seedIn); }

  int operator()() { return m_dist(m_gen); }
};

template <class SourceViewType, class DestViewType, class DistancesViewType,
          class IntraTeamSentinelView, class ValueType>
struct TestFunctorA {
  SourceViewType m_sourceView;
  DestViewType m_destView;
  ValueType m_targetValue;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const SourceViewType sourceView, const DestViewType destView,
               ValueType targetVal, const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_sourceView(sourceView),
        m_destView(destView),
        m_targetValue(targetVal),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowViewFrom =
        Kokkos::subview(m_sourceView, myRowIndex, Kokkos::ALL());
    auto myRowViewDest = Kokkos::subview(m_destView, myRowIndex, Kokkos::ALL());
    ptrdiff_t resultDist = 0;

    if (m_apiPick == 0) {
      auto it    = KE::remove_copy(member, KE::cbegin(myRowViewFrom),
                                KE::cend(myRowViewFrom),
                                KE::begin(myRowViewDest), m_targetValue);
      resultDist = KE::distance(KE::begin(myRowViewDest), it);
      Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
        m_distancesView(myRowIndex) = resultDist;
      });
    } else if (m_apiPick == 1) {
      auto it =
          KE::remove_copy(member, myRowViewFrom, myRowViewDest, m_targetValue);
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
     set a random subset of each row of a rank-2 view
     to a target value and run a team-level KE::remove_copy
     to a destination view with one team per row to remove all those elements.
   */

  const auto targetVal = static_cast<ValueType>(531);

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // construct in memory space associated with default exespace
  auto sourceView =
      create_view<ValueType>(LayoutTag{}, numTeams, numCols, "dataView");

  // sourceView might not deep copyable (e.g. strided layout) so to fill it
  // we make a new view that is for sure deep copyable, modify it on the host
  // deep copy to device and then launch copy kernel to sourceView
  auto sourceView_dc =
      create_deep_copyable_compatible_view_with_same_extent(sourceView);
  auto sourceView_dc_h = create_mirror_view(Kokkos::HostSpace(), sourceView_dc);

  Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace> pool(
      45234399);
  Kokkos::fill_random(sourceView_dc_h, pool, ValueType(0), ValueType(1177));

  // for each row, randomly select columns, fill with targetVal
  std::vector<std::size_t> perRowRealCount(numTeams);
  const std::size_t maxColInd = numCols > 0 ? numCols - 1 : 0;
  UnifDist<int> colCountProducer(maxColInd, 3123377);
  UnifDist<int> colIndicesProducer(maxColInd, 455225);
  for (std::size_t i = 0; i < sourceView_dc_h.extent(0); ++i) {
    const std::size_t currCount = colCountProducer();
    std::vector<std::size_t> colIndForThisRow(currCount);
    for (std::size_t j = 0; j < currCount; ++j) {
      const auto colInd          = colIndicesProducer();
      sourceView_dc_h(i, colInd) = targetVal;
      colIndForThisRow[j]        = colInd;
    }

    // note that we need to count how many elements are equal
    // to targetVal because the sourceView was origianlly filled
    // with random values so it could be that we have more matches
    // than what we manually set above
    std::size_t realCount = 0;
    for (std::size_t j = 0; j < sourceView_dc_h.extent(1); ++j) {
      if (sourceView_dc_h(i, j) == targetVal) {
        realCount++;
      }
    }
    perRowRealCount[i] = realCount;
  }

  // copy to sourceView_dc and then to sourceView
  Kokkos::deep_copy(sourceView_dc, sourceView_dc_h);
  CopyFunctorRank2 F1(sourceView_dc, sourceView);
  Kokkos::parallel_for("copy", sourceView.extent(0) * sourceView.extent(1), F1);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // create the destination view
  Kokkos::View<ValueType**> destView("destView", numTeams, numCols);

  // each team stores the distance of the returned iterator from the
  // beginning of the interval that team operates on and then we check
  // that these distances match the std result
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // use CTAD for functor
  TestFunctorA fnc(sourceView, destView, targetVal, distancesView,
                   intraTeamSentinelView, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check against std
  // -----------------------------------------------
  auto destViewAfterOp_h       = create_host_space_copy(destView);
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  Kokkos::View<ValueType**, Kokkos::HostSpace> stdDestView("stdDestView",
                                                           numTeams, numCols);

  for (std::size_t i = 0; i < destViewAfterOp_h.extent(0); ++i) {
    auto rowFrom = Kokkos::subview(sourceView_dc_h, i, Kokkos::ALL());
    auto rowDest = Kokkos::subview(stdDestView, i, Kokkos::ALL());

    auto stdIt = std::remove_copy(KE::cbegin(rowFrom), KE::cend(rowFrom),
                                  KE::begin(rowDest), targetVal);
    const std::size_t stdDistance = KE::distance(KE::begin(rowDest), stdIt);

    EXPECT_TRUE(distancesView_h(i) == stdDistance);
    // EXPECT_TRUE(distancesView_h(i) == numCols - perRowRealCount[i]);
    for (std::size_t j = 0; j < distancesView_h(i); ++j) {
      EXPECT_TRUE(destViewAfterOp_h(i, j) == stdDestView(i, j));
    }
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios() {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 8113}) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_remove_copy_team_test, test) {
// FIXME_OPENMPTARGET
#if defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(KOKKOS_ARCH_INTEL_GPU)
  GTEST_SKIP() << "the test is known to fail with OpenMPTarget on Intel GPUs";
#endif
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamRemoveCopy
}  // namespace stdalgos
}  // namespace Test
