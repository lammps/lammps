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
namespace TeamFindIfNot {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct GreaterEqualFunctor {
  ValueType m_val;

  KOKKOS_INLINE_FUNCTION
  GreaterEqualFunctor(ValueType val) : m_val(val) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(ValueType val) const { return (val >= m_val); }
};

template <class DataViewType, class GreaterThanValuesViewType,
          class DistancesViewType, class IntraTeamSentinelView>
struct TestFunctorA {
  DataViewType m_dataView;
  GreaterThanValuesViewType m_greaterThanValuesView;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;

  TestFunctorA(const DataViewType dataView,
               const GreaterThanValuesViewType greaterThanValuesView,
               DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick)
      : m_dataView(dataView),
        m_greaterThanValuesView(greaterThanValuesView),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowViewFrom = Kokkos::subview(m_dataView, myRowIndex, Kokkos::ALL());
    const auto val     = m_greaterThanValuesView(myRowIndex);
    // FIXME_INTEL
#if defined(KOKKOS_COMPILER_INTEL) && (1900 == KOKKOS_COMPILER_INTEL)
    GreaterEqualFunctor<
        typename GreaterThanValuesViewType::non_const_value_type>
        unaryPred{val};
#else
    GreaterEqualFunctor unaryPred{val};
#endif
    ptrdiff_t resultDist = 0;

    switch (m_apiPick) {
      case 0: {
        auto it    = KE::find_if_not(member, KE::cbegin(myRowViewFrom),
                                  KE::cend(myRowViewFrom), unaryPred);
        resultDist = KE::distance(KE::cbegin(myRowViewFrom), it);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_distancesView(myRowIndex) = resultDist;
        });

        break;
      }

      case 1: {
        auto it    = KE::find_if_not(member, myRowViewFrom, unaryPred);
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
void test_A(const bool predicatesReturnTrue, std::size_t numTeams,
            std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level find_if_not
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

  // find_if_not returns an iterator so to verify that it is correct
  // each team stores the distance of the returned iterator from the
  // beginning of the interval that team operates on and then we check
  // that these distances match the std result
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  // If predicatesReturnTrue == true, we want to ensure that for each dataView's
  // row find_if_not always returns end iterator. To do that,
  // GreaterEqualFunctor predicate created for each row must return true for
  // every value in that row, so it needs to compare each value with value
  // smaller than lowerBound.
  //
  // If predicatesReturnTrue == false we want to ensure the opposite -
  // GreaterEqualFunctor needs to return false for every value of each
  // dataView's row, so the predicate is constructed with value randomly picked
  // from range [upperBound, upperBound*2).
  Kokkos::View<ValueType*> greaterEqualValuesView("greaterEqualValuesView",
                                                  numTeams);
  auto greaterEqualValuesView_h =
      create_mirror_view(Kokkos::HostSpace(), greaterEqualValuesView);

  using rand_pool =
      Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>;
  rand_pool pool(lowerBound * upperBound);

  if (predicatesReturnTrue) {
    Kokkos::fill_random(greaterEqualValuesView_h, pool, 0, lowerBound);
  } else {
    Kokkos::fill_random(greaterEqualValuesView_h, pool, upperBound,
                        upperBound * 2);
  }

  Kokkos::deep_copy(greaterEqualValuesView, greaterEqualValuesView_h);

  // use CTAD for functor
  TestFunctorA fnc(dataView, greaterEqualValuesView, distancesView,
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
    const auto val          = greaterEqualValuesView_h(i);
    // FIXME_INTEL
#if defined(KOKKOS_COMPILER_INTEL) && (1900 == KOKKOS_COMPILER_INTEL)
    const GreaterEqualFunctor<ValueType> unaryPred{val};
#else
    const GreaterEqualFunctor unaryPred{val};
#endif

    auto it = std::find_if_not(rowFromBegin, rowFromEnd, unaryPred);

    const std::size_t stdDistance      = KE::distance(rowFromBegin, it);
    const std::size_t beginEndDistance = KE::distance(rowFromBegin, rowFromEnd);

    if (predicatesReturnTrue) {
      ASSERT_EQ(stdDistance, beginEndDistance);
    } else {
      EXPECT_LT(stdDistance, beginEndDistance);
    }

    ASSERT_EQ(stdDistance, distancesView_h(i));
    ASSERT_TRUE(intraTeamSentinelView_h(i));
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios(const bool predicatesReturnTrue) {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(predicatesReturnTrue, numTeams, numCols,
                                     apiId);
      }
    }
  }
}

TEST(std_algorithms_find_if_not_team_test, predicates_return_true) {
  constexpr bool predicatesReturnTrue = true;

  run_all_scenarios<DynamicTag, double>(predicatesReturnTrue);
  run_all_scenarios<StridedTwoRowsTag, int>(predicatesReturnTrue);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(predicatesReturnTrue);
}

TEST(std_algorithms_find_if_not_team_test, predicates_return_false) {
  constexpr bool predicatesReturnTrue = false;

  run_all_scenarios<DynamicTag, double>(predicatesReturnTrue);
  run_all_scenarios<StridedTwoRowsTag, int>(predicatesReturnTrue);
  run_all_scenarios<StridedThreeRowsTag, unsigned>(predicatesReturnTrue);
}

}  // namespace TeamFindIfNot
}  // namespace stdalgos
}  // namespace Test
