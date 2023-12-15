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

#if not defined KOKKOS_ENABLE_OPENMPTARGET

namespace Test {
namespace stdalgos {
namespace TeamTransformInclusiveScan {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct PlusFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& lhs, const ValueType& rhs) const {
    return lhs + rhs;
  }
};

template <class ValueType>
struct MultipliesByTwoFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& value) const { return value * 2; }
};

template <class SourceViewType, class DestViewType, class DistancesViewType,
          class IntraTeamSentinelView, class InitValuesViewType,
          class BinaryOpType, class UnaryOpType>
struct TestFunctorA {
  SourceViewType m_sourceView;
  DestViewType m_destView;
  DistancesViewType m_distancesView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  InitValuesViewType m_initValuesView;
  BinaryOpType m_binaryOp;
  UnaryOpType m_unaryOp;
  int m_apiPick;

  TestFunctorA(const SourceViewType sourceView, const DestViewType destView,
               const DistancesViewType distancesView,
               const IntraTeamSentinelView intraTeamSentinelView,
               const InitValuesViewType initValuesView, BinaryOpType binaryOp,
               UnaryOpType unaryOp, int apiPick)
      : m_sourceView(sourceView),
        m_destView(destView),
        m_distancesView(distancesView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_initValuesView(initValuesView),
        m_binaryOp(binaryOp),
        m_unaryOp(unaryOp),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto rowIndex = member.league_rank();

    auto srcRow      = Kokkos::subview(m_sourceView, rowIndex, Kokkos::ALL());
    const auto first = KE::cbegin(srcRow);
    const auto last  = KE::cend(srcRow);
    auto destRow     = Kokkos::subview(m_destView, rowIndex, Kokkos::ALL());
    auto firstDest   = KE::begin(destRow);

    const auto initVal   = m_initValuesView(rowIndex);
    ptrdiff_t resultDist = 0;

    switch (m_apiPick) {
      case 0: {
        auto it = KE::transform_inclusive_scan(member, first, last, firstDest,
                                               m_binaryOp, m_unaryOp);
        resultDist = KE::distance(firstDest, it);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this] { m_distancesView(rowIndex) = resultDist; });

        break;
      }

      case 1: {
        auto it    = KE::transform_inclusive_scan(member, srcRow, destRow,
                                               m_binaryOp, m_unaryOp);
        resultDist = KE::distance(firstDest, it);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this] { m_distancesView(rowIndex) = resultDist; });

        break;
      }

      case 2: {
        auto it = KE::transform_inclusive_scan(member, first, last, firstDest,
                                               m_binaryOp, m_unaryOp, initVal);
        resultDist = KE::distance(firstDest, it);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this] { m_distancesView(rowIndex) = resultDist; });

        break;
      }

      case 3: {
        auto it    = KE::transform_inclusive_scan(member, srcRow, destRow,
                                               m_binaryOp, m_unaryOp, initVal);
        resultDist = KE::distance(firstDest, it);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this] { m_distancesView(rowIndex) = resultDist; });

        break;
      }
    }

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck = team_members_have_matching_result(
        member, resultDist, m_distancesView(rowIndex));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(rowIndex) = intraTeamCheck;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level transform_inclusive_scan
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

  auto [sourceView, sourceViewBeforeOp_h] = create_random_view_and_host_clone(
      LayoutTag{}, numTeams, numCols, bounds, "sourceView");

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // create the destination view
  Kokkos::View<ValueType**> destView("destView", numTeams, numCols);

  // tranform_inclusive_scan returns an iterator so to verify that it is correct
  // each team stores the distance of the returned iterator from the beginning
  // of the interval that team operates on and then we check that these
  // distances match the std result
  Kokkos::View<std::size_t*> distancesView("distancesView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  PlusFunctor<ValueType> binaryOp;
  MultipliesByTwoFunctor<ValueType> unaryOp;

  // Create view of reduce init values to be used by test cases
  Kokkos::View<ValueType*, Kokkos::DefaultHostExecutionSpace> initValuesView_h(
      "initValuesView_h", numTeams);
  using rand_pool =
      Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>;
  rand_pool pool(lowerBound * upperBound);
  Kokkos::fill_random(initValuesView_h, pool, lowerBound, upperBound);

  // use CTAD for functor
  auto initValuesView =
      Kokkos::create_mirror_view_and_copy(space_t(), initValuesView_h);
  TestFunctorA fnc(sourceView, destView, distancesView, intraTeamSentinelView,
                   initValuesView, binaryOp, unaryOp, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel and check
  // -----------------------------------------------
  auto distancesView_h         = create_host_space_copy(distancesView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);
  Kokkos::View<ValueType**, Kokkos::HostSpace> stdDestView("stdDestView",
                                                           numTeams, numCols);

  for (std::size_t i = 0; i < sourceView.extent(0); ++i) {
    auto srcRow    = Kokkos::subview(sourceViewBeforeOp_h, i, Kokkos::ALL());
    auto first     = KE::cbegin(srcRow);
    auto last      = KE::cend(srcRow);
    auto destRow   = Kokkos::subview(stdDestView, i, Kokkos::ALL());
    auto firstDest = KE::begin(destRow);
    auto initValue = initValuesView_h(i);

    ASSERT_TRUE(intraTeamSentinelView_h(i));

// libstdc++ as provided by GCC 8 does not have transform_inclusive_scan and
// for GCC 9.1, 9.2 fails to compile for missing overload not accepting policy
#if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE <= 9)
#define transform_inclusive_scan testing_transform_inclusive_scan
#else
#define transform_inclusive_scan std::transform_inclusive_scan
#endif

    switch (apiId) {
      case 0:
      case 1: {
        const auto it =
            transform_inclusive_scan(first, last, firstDest, binaryOp, unaryOp);
        const std::size_t stdDistance = KE::distance(firstDest, it);
        ASSERT_EQ(stdDistance, distancesView_h(i));

        break;
      }

      case 2:
      case 3: {
        const auto it = transform_inclusive_scan(first, last, firstDest,
                                                 binaryOp, unaryOp, initValue);
        const std::size_t stdDistance = KE::distance(firstDest, it);
        ASSERT_EQ(stdDistance, distancesView_h(i));

        break;
      }
    }
  }
#undef transform_inclusive_scan

  auto dataViewAfterOp_h = create_host_space_copy(destView);
  expect_equal_host_views(stdDestView, dataViewAfterOp_h);
}

template <class LayoutTag, class ValueType>
void run_all_scenarios() {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1, 2, 3}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_transform_inclusive_scan_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamTransformInclusiveScan
}  // namespace stdalgos
}  // namespace Test

#endif
