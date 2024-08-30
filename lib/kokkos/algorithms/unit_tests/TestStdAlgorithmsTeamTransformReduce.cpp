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
namespace TeamTransformReduce {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct PlusFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& lhs, const ValueType& rhs) const {
    return lhs + rhs;
  }
};

template <class ValueType>
struct MultipliesFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& lhs, const ValueType& rhs) const {
    return lhs * rhs;
  }
};

template <class ValueType>
struct PlusOneFunctor {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()(const ValueType& val) const { return val + 1; };
};

template <class FirstDataViewType, class SecondDataViewType,
          class InitValuesViewType, class ResultsViewType,
          class IntraTeamSentinelView, class BinaryJoinerType,
          class BinaryTransformType, class UnaryTransformType>
struct TestFunctorA {
  FirstDataViewType m_firstDataView;
  SecondDataViewType m_secondDataView;
  InitValuesViewType m_initValuesView;
  ResultsViewType m_resultsView;
  IntraTeamSentinelView m_intraTeamSentinelView;
  BinaryJoinerType m_binaryJoiner;
  BinaryTransformType m_binaryTransform;
  UnaryTransformType m_unaryTransform;
  int m_apiPick;

  TestFunctorA(const FirstDataViewType firstDataView,
               const SecondDataViewType secondDataview,
               const InitValuesViewType initValuesView,
               const ResultsViewType resultsView,
               const IntraTeamSentinelView intraTeamSentinelView,
               BinaryJoinerType binaryJoiner,
               BinaryTransformType binaryTransform,
               UnaryTransformType unaryTransform, int apiPick)
      : m_firstDataView(firstDataView),
        m_secondDataView(secondDataview),
        m_initValuesView(initValuesView),
        m_resultsView(resultsView),
        m_intraTeamSentinelView(intraTeamSentinelView),
        m_binaryJoiner(binaryJoiner),
        m_binaryTransform(binaryTransform),
        m_unaryTransform(unaryTransform),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const int rowIndex = member.league_rank();

    auto firstDataRow =
        Kokkos::subview(m_firstDataView, rowIndex, Kokkos::ALL());
    auto firstDataRowBegin = KE::cbegin(firstDataRow);
    auto firstDataRowEnd   = KE::cend(firstDataRow);

    auto secondDataRow =
        Kokkos::subview(m_secondDataView, rowIndex, Kokkos::ALL());
    auto secondDataRowBegin = KE::cbegin(secondDataRow);

    const auto initVal = m_initValuesView(rowIndex);
    typename InitValuesViewType::non_const_value_type result = 0;

    switch (m_apiPick) {
      case 0: {
        result =
            KE::transform_reduce(member, firstDataRowBegin, firstDataRowEnd,
                                 secondDataRowBegin, initVal);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 1: {
        result =
            KE::transform_reduce(member, firstDataRow, secondDataRow, initVal);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 2: {
        result = KE::transform_reduce(
            member, firstDataRowBegin, firstDataRowEnd, secondDataRowBegin,
            initVal, m_binaryJoiner, m_binaryTransform);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 3: {
        result =
            KE::transform_reduce(member, firstDataRow, secondDataRow, initVal,
                                 m_binaryJoiner, m_binaryTransform);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 4: {
        result =
            KE::transform_reduce(member, firstDataRowBegin, firstDataRowEnd,
                                 initVal, m_binaryJoiner, m_unaryTransform);
        Kokkos::single(Kokkos::PerTeam(member),
                       [=, *this]() { m_resultsView(rowIndex) = result; });
        break;
      }

      case 5: {
        result = KE::transform_reduce(member, firstDataRow, initVal,
                                      m_binaryJoiner, m_unaryTransform);
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
void test_A(std::size_t numTeams, std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level transform_reduce
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

  auto [firstDataView, firstDataViewBeforeOp_h] =
      create_random_view_and_host_clone(LayoutTag{}, numTeams, numCols, bounds,
                                        "firstDataView");
  auto [secondDataView, secondDataViewBeforeOp_h] =
      create_random_view_and_host_clone(LayoutTag{}, numTeams, numCols, bounds,
                                        "secondDataView");

  // Create view of init values to be used by test cases
  Kokkos::View<ValueType*, Kokkos::DefaultHostExecutionSpace> initValuesView_h(
      "initValuesView_h", numTeams);
  using rand_pool =
      Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace>;
  rand_pool pool(lowerBound * upperBound);
  Kokkos::fill_random(initValuesView_h, pool, lowerBound, upperBound);

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // to verify that things work, each team stores the result of its
  // transform_reduce call, and then we check that these match what we expect
  Kokkos::View<ValueType*> resultsView("resultsView", numTeams);
  // sentinel to check if all members of the team compute the same result
  Kokkos::View<bool*> intraTeamSentinelView("intraTeamSameResult", numTeams);

  PlusFunctor<ValueType> binaryJoiner;
  MultipliesFunctor<ValueType> binaryTransform;
  PlusOneFunctor<ValueType> unaryTransform;

  // use CTAD for functor
  auto initValuesView =
      Kokkos::create_mirror_view_and_copy(space_t(), initValuesView_h);
  TestFunctorA fnc(firstDataView, secondDataView, initValuesView, resultsView,
                   intraTeamSentinelView, binaryJoiner, binaryTransform,
                   unaryTransform, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel and check
  // -----------------------------------------------

  auto resultsView_h           = create_host_space_copy(resultsView);
  auto intraTeamSentinelView_h = create_host_space_copy(intraTeamSentinelView);

  for (std::size_t i = 0; i < firstDataView.extent(0); ++i) {
    auto firstDataRow =
        Kokkos::subview(firstDataViewBeforeOp_h, i, Kokkos::ALL());

    const auto firstDataRowBegin = KE::cbegin(firstDataRow);
    const auto firstDataRowEnd   = KE::cend(firstDataRow);

    auto secondDataRow =
        Kokkos::subview(secondDataViewBeforeOp_h, i, Kokkos::ALL());

    const auto secondDataRowBegin = KE::cbegin(secondDataRow);

    const auto initVal = initValuesView_h(i);

    ASSERT_TRUE(intraTeamSentinelView_h(i));

// libstdc++ as provided by GCC 8 does not have reduce, transform_reduce,
// exclusive_scan, inclusive_scan, transform_exclusive_scan,
// transform_inclusive_scan and for GCC 9.1, 9.2 fails to compile them for
// missing overload not accepting policy
#if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE <= 9)
#define transform_reduce testing_transform_reduce
#else
#define transform_reduce std::transform_reduce
#endif

    switch (apiId) {
      case 0:
      case 1: {
        const auto result = transform_reduce(firstDataRowBegin, firstDataRowEnd,
                                             secondDataRowBegin, initVal);

        if constexpr (std::is_floating_point_v<ValueType>) {
          EXPECT_FLOAT_EQ(result, resultsView_h(i));
        } else {
          ASSERT_EQ(result, resultsView_h(i));
        }

        break;
      }

      case 2:
      case 3: {
        const ValueType result = transform_reduce(
            firstDataRowBegin, firstDataRowEnd, secondDataRowBegin, initVal,
            binaryJoiner, binaryTransform);

        if constexpr (std::is_floating_point_v<ValueType>) {
          EXPECT_FLOAT_EQ(result, resultsView_h(i));
        } else {
          ASSERT_EQ(result, resultsView_h(i));
        }

        break;
      }

      case 4:
      case 5: {
        const ValueType result =
            transform_reduce(firstDataRowBegin, firstDataRowEnd, initVal,
                             binaryJoiner, unaryTransform);

        if constexpr (std::is_floating_point_v<ValueType>) {
          EXPECT_FLOAT_EQ(result, resultsView_h(i));
        } else {
          ASSERT_EQ(result, resultsView_h(i));
        }

        break;
      }
    }

#undef transform_reduce
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

TEST(std_algorithms_transform_reduce_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamTransformReduce
}  // namespace stdalgos
}  // namespace Test

#endif
