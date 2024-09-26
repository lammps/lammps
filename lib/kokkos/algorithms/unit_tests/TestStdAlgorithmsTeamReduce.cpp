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

#ifndef KOKKOS_ENABLE_OPENMPTARGET

namespace Test {
namespace stdalgos {
namespace TeamReduce {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct PlusFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& lhs, const ValueType& rhs) const {
    return lhs + rhs;
  }
};

template <class DataViewType, class ReductionInitValuesViewType,
          class ReduceResultsViewType, class IntraTeamSentinelView,
          class BinaryPredType>
struct TestFunctorA {
  DataViewType m_dataView;
  ReductionInitValuesViewType m_reductionInitValuesView;
  ReduceResultsViewType m_reduceResultsView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  int m_apiPick;
  BinaryPredType m_binaryPred;

  TestFunctorA(const DataViewType dataView,
               const ReductionInitValuesViewType reductionInitValuesView,
               const ReduceResultsViewType reduceResultsView,
               const IntraTeamSentinelView intraTeamSentinelView, int apiPick,
               BinaryPredType binaryPred)
      : m_dataView(dataView),
        m_reductionInitValuesView(reductionInitValuesView),
        m_reduceResultsView(reduceResultsView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_apiPick(apiPick),
        m_binaryPred(binaryPred) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();

    auto myRowViewFrom = Kokkos::subview(m_dataView, myRowIndex, Kokkos::ALL());

    const auto rowFromBegin     = KE::cbegin(myRowViewFrom);
    const auto rowFromEnd       = KE::cend(myRowViewFrom);
    const auto initReductionVal = m_reductionInitValuesView(myRowIndex);
    typename ReduceResultsViewType::non_const_value_type result = 0;

    switch (m_apiPick) {
      case 0: {
        result = KE::reduce(member, rowFromBegin, rowFromEnd);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_reduceResultsView(myRowIndex) = result;
        });
        break;
      }

      case 1: {
        result = KE::reduce(member, myRowViewFrom);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_reduceResultsView(myRowIndex) = result;
        });
        break;
      }

      case 2: {
        result = KE::reduce(member, rowFromBegin, rowFromEnd, initReductionVal);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_reduceResultsView(myRowIndex) = result;
        });
        break;
      }

      case 3: {
        result = KE::reduce(member, myRowViewFrom, initReductionVal);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_reduceResultsView(myRowIndex) = result;
        });
        break;
      }

      case 4: {
        result = KE::reduce(member, rowFromBegin, rowFromEnd, initReductionVal,
                            m_binaryPred);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_reduceResultsView(myRowIndex) = result;
        });
        break;
      }

      case 5: {
        result =
            KE::reduce(member, myRowViewFrom, initReductionVal, m_binaryPred);
        Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
          m_reduceResultsView(myRowIndex) = result;
        });
        break;
      }
    }

    // store result of checking if all members have their local
    // values matching the one stored in m_distancesView
    member.team_barrier();
    const bool intraTeamCheck = team_members_have_matching_result(
        member, result, m_reduceResultsView(myRowIndex));
    Kokkos::single(Kokkos::PerTeam(member), [=, *this]() {
      m_intraTeamSentinelView(myRowIndex) = intraTeamCheck;
    });
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level reduce
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

  // Create view of reduce init values to be used by test cases
  Kokkos::View<ValueType*, Kokkos::DefaultHostExecutionSpace>
      reductionInitValuesView_h("reductionInitValuesView_h", numTeams);
  using rand_pool =
      Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>;
  rand_pool pool(lowerBound * upperBound);
  Kokkos::fill_random(reductionInitValuesView_h, pool, lowerBound, upperBound);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // to verify that things work, each team stores the result of its reduce
  // call, and then we check that these match what we expect
  Kokkos::View<ValueType*> reduceResultsView("reduceResultsView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  PlusFunctor<ValueType> binaryPred;

  // use CTAD for functor
  auto reductionInitValuesView =
      Kokkos::create_mirror_view_and_copy(space_t(), reductionInitValuesView_h);
  TestFunctorA fnc(dataView, reductionInitValuesView, reduceResultsView,
                   intraTeamSentinelView, apiId, binaryPred);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel and check
  // -----------------------------------------------

  auto reduceResultsView_h     = create_host_space_copy(reduceResultsView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);

  for (std::size_t i = 0; i < dataView.extent(0); ++i) {
    auto rowFrom = Kokkos::subview(dataViewBeforeOp_h, i, Kokkos::ALL());

    const auto rowFromBegin = KE::cbegin(rowFrom);
    const auto rowFromEnd   = KE::cend(rowFrom);
    const auto initVal      = reductionInitValuesView_h(i);

    ASSERT_TRUE(intraTeamSentinelView_h(i));

    // libstdc++ as provided by GCC 8 does not have reduce, transform_reduce,
    // exclusive_scan, inclusive_scan, transform_exclusive_scan,
    // transform_inclusive_scan and for GCC 9.1, 9.2 fails to compile them for
    // missing overload not accepting policy

#if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE <= 9)
#define reduce testing_reduce
#else
#define reduce std::reduce
#endif

    switch (apiId) {
      case 0:
      case 1: {
        const ValueType result = reduce(rowFromBegin, rowFromEnd);
        if constexpr (std::is_floating_point_v<ValueType>) {
          EXPECT_FLOAT_EQ(result, reduceResultsView_h(i));
        } else {
          ASSERT_EQ(result, reduceResultsView_h(i));
        }

        break;
      }

      case 2:
      case 3: {
        const ValueType result = reduce(rowFromBegin, rowFromEnd, initVal);
        if constexpr (std::is_floating_point_v<ValueType>) {
          EXPECT_FLOAT_EQ(result, reduceResultsView_h(i));
        } else {
          ASSERT_EQ(result, reduceResultsView_h(i));
        }

        break;
      }

      case 4:
      case 5: {
        const ValueType result =
            reduce(rowFromBegin, rowFromEnd, initVal, binaryPred);
        if constexpr (std::is_floating_point_v<ValueType>) {
          EXPECT_FLOAT_EQ(result, reduceResultsView_h(i));
        } else {
          ASSERT_EQ(result, reduceResultsView_h(i));
        }

        break;
      }
    }

#undef reduce
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios() {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1, 2, 3, 4, 5}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_reduce_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamReduce
}  // namespace stdalgos
}  // namespace Test

#endif
