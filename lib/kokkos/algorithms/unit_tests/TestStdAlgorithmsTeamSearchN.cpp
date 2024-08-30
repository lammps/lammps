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
namespace TeamSearchN {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct EqualFunctor {
  KOKKOS_INLINE_FUNCTION
  bool operator()(const ValueType& lhs, const ValueType& rhs) const {
    return lhs == rhs;
  }
};

template <class DataViewType, class SearchedValuesViewType,
          class DistancesViewType, class IntraTeamSentinelView,
          class BinaryPredType>
struct TestFunctorA {
  DataViewType m_dataView;
  std::size_t m_seqSize;
  SearchedValuesViewType m_searchedValuesView;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  BinaryPredType m_binaryPred;
  int m_apiPick;

  TestFunctorA(const DataViewType dataView, std::size_t seqSize,
               const SearchedValuesViewType searchedValuesView,
               const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView,
               BinaryPredType binaryPred, int apiPick)
      : m_dataView(dataView),
        m_seqSize(seqSize),
        m_searchedValuesView(searchedValuesView),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_binaryPred(binaryPred),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowViewFrom = Kokkos::subview(m_dataView, myRowIndex, Kokkos::ALL());
    auto rowFromBegin  = KE::begin(myRowViewFrom);
    auto rowFromEnd    = KE::end(myRowViewFrom);
    const auto searchedValue = m_searchedValuesView(myRowIndex);
    ptrdiff_t resultDist     = 0;

    switch (m_apiPick) {
      case 0: {
        const auto it = KE::search_n(member, rowFromBegin, rowFromEnd,
                                     m_seqSize, searchedValue);
        resultDist    = KE::distance(rowFromBegin, it);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_distancesView(myRowIndex) = resultDist;
        });

        break;
      }

      case 1: {
        const auto it =
            KE::search_n(member, myRowViewFrom, m_seqSize, searchedValue);
        resultDist = KE::distance(rowFromBegin, it);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_distancesView(myRowIndex) = resultDist;
        });

        break;
      }

      case 2: {
        const auto it = KE::search_n(member, rowFromBegin, rowFromEnd,
                                     m_seqSize, searchedValue, m_binaryPred);
        resultDist    = KE::distance(rowFromBegin, it);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_distancesView(myRowIndex) = resultDist;
        });

        break;
      }

      case 3: {
        const auto it = KE::search_n(member, myRowViewFrom, m_seqSize,
                                     searchedValue, m_binaryPred);
        resultDist    = KE::distance(rowFromBegin, it);
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
void test_A(const bool sequencesExist, std::size_t numTeams,
            std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level search_n
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

  // If sequencesExist == true we need to inject some sequence of count test
  // value into dataView. If sequencesExist == false we set searchedVal to a
  // value that is not present in dataView

  const std::size_t halfCols = (numCols > 1) ? ((numCols + 1) / 2) : (1);
  const std::size_t seqSize  = (numCols > 1) ? (std::log2(numCols)) : (1);

  Kokkos::View<ValueType*> searchedValuesView("searchedValuesView", numTeams);
  auto searchedValuesView_h = create_host_space_copy(searchedValuesView);

  // dataView might not deep copyable (e.g. strided layout) so to prepare it
  // correclty, we make a new view that is for sure deep copyable, modify it
  // on the host, deep copy to device and then launch a kernel to copy to
  // dataView
  auto dataView_dc =
      create_deep_copyable_compatible_view_with_same_extent(dataView);
  auto dataView_dc_h = create_mirror_view(Kokkos::HostSpace(), dataView_dc);

  if (sequencesExist) {
    const std::size_t dataBegin = halfCols - seqSize;
    for (std::size_t i = 0; i < searchedValuesView.extent(0); ++i) {
      const ValueType searchedVal = dataView_dc_h(i, dataBegin);
      searchedValuesView_h(i)     = searchedVal;

      for (std::size_t j = dataBegin + 1; j < seqSize; ++j) {
        dataView_dc_h(i, j) = searchedVal;
      }
    }

    // copy to dataView_dc and then to dataView
    Kokkos::deep_copy(dataView_dc, dataView_dc_h);

    CopyFunctorRank2 cpFun(dataView_dc, dataView);
    Kokkos::parallel_for("copy", dataView.extent(0) * dataView.extent(1),
                         cpFun);
  } else {
    using rand_pool =
        Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>;
    rand_pool pool(lowerBound * upperBound);
    Kokkos::fill_random(searchedValuesView_h, pool, upperBound, upperBound * 2);
  }

  Kokkos::deep_copy(searchedValuesView, searchedValuesView_h);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // search_n returns an iterator so to verify that it is correct each team
  // stores the distance of the returned iterator from the beginning of the
  // interval that team operates on and then we check that these distances match
  // the std result
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  EqualFunctor<ValueType> binaryPred;

  // use CTAD for functor
  TestFunctorA fnc(dataView, seqSize, searchedValuesView, distancesView,
                   intraTeamSentinelView, binaryPred, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel and check
  // -----------------------------------------------
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);

  for (std::size_t i = 0; i < dataView.extent(0); ++i) {
    auto rowFrom = Kokkos::subview(dataView_dc_h, i, Kokkos::ALL());

    const auto rowFromBegin = KE::cbegin(rowFrom);
    const auto rowFromEnd   = KE::cend(rowFrom);

    const ValueType searchedVal = searchedValuesView_h(i);

    const std::size_t beginEndDist = KE::distance(rowFromBegin, rowFromEnd);

    switch (apiId) {
      case 0:
      case 1: {
        const auto it =
            std::search_n(rowFromBegin, rowFromEnd, seqSize, searchedVal);
        const std::size_t stdDistance = KE::distance(rowFromBegin, it);

        if (sequencesExist) {
          EXPECT_LT(distancesView_h(i), beginEndDist);
        } else {
          ASSERT_EQ(distancesView_h(i), beginEndDist);
        }

        ASSERT_EQ(stdDistance, distancesView_h(i));
        ASSERT_TRUE(intraTeamSentinelView_h(i));
        break;
      }

      case 2:
      case 3: {
        const auto it = std::search_n(rowFromBegin, rowFromEnd, seqSize,
                                      searchedVal, binaryPred);
        const std::size_t stdDistance = KE::distance(rowFromBegin, it);

        if (sequencesExist) {
          EXPECT_LT(distancesView_h(i), beginEndDist);
        } else {
          ASSERT_EQ(distancesView_h(i), beginEndDist);
        }

        ASSERT_EQ(stdDistance, distancesView_h(i));
        ASSERT_TRUE(intraTeamSentinelView_h(i));

        break;
      }
    }
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios(const bool sequencesExist) {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(sequencesExist, numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_search_n_team_test, sequences_of_equal_elements_exist) {
  constexpr bool sequencesExist = true;

  run_all_scenarios<DynamicTag, double>(sequencesExist);
  run_all_scenarios<StridedTwoRowsTag, int>(sequencesExist);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(sequencesExist);
}

TEST(std_algorithms_search_n_team_test,
     sequences_of_equal_elements_probably_does_not_exist) {
  constexpr bool sequencesExist = false;

  run_all_scenarios<DynamicTag, double>(sequencesExist);
  run_all_scenarios<StridedTwoRowsTag, int>(sequencesExist);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(sequencesExist);
}

}  // namespace TeamSearchN
}  // namespace stdalgos
}  // namespace Test
