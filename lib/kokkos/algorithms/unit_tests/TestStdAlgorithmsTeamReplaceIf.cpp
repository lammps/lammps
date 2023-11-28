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
namespace TeamReplaceIf {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct GreaterThanValueFunctor {
  ValueType m_val;

  KOKKOS_INLINE_FUNCTION
  GreaterThanValueFunctor(ValueType val) : m_val(val) {}

  KOKKOS_INLINE_FUNCTION
  bool operator()(ValueType val) const { return (val > m_val); }
};

template <class ViewType, class ValueType>
struct TestFunctorA {
  ViewType m_view;
  ValueType m_threshold;
  ValueType m_newVal;
  int m_apiPick;

  TestFunctorA(const ViewType view, ValueType threshold, ValueType newVal,
               int apiPick)
      : m_view(view),
        m_threshold(threshold),
        m_newVal(newVal),
        m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowView        = Kokkos::subview(m_view, myRowIndex, Kokkos::ALL());

    GreaterThanValueFunctor predicate(m_threshold);
    if (m_apiPick == 0) {
      KE::replace_if(member, KE::begin(myRowView), KE::end(myRowView),
                     predicate, m_newVal);
    } else if (m_apiPick == 1) {
      KE::replace_if(member, myRowView, predicate, m_newVal);
    }
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId) {
  /* description:
     use a rank-2 view randomly filled with values,
     and run a team-level replace_if where the values strictly greater
     than a threshold are replaced with a new value.
   */
  const auto threshold = static_cast<ValueType>(151);
  const auto newVal    = static_cast<ValueType>(1);

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // create a view in the memory space associated with default exespace
  // with as many rows as the number of teams and fill it with random
  // values from an arbitrary range
  auto [dataView, cloneOfDataViewBeforeOp_h] =
      create_random_view_and_host_clone(
          LayoutTag{}, numTeams, numCols,
          Kokkos::pair<ValueType, ValueType>{5, 523}, "dataView");

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());
  // use CTAD for functor
  TestFunctorA fnc(dataView, threshold, newVal, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // run cpp-std kernel
  // -----------------------------------------------
  Kokkos::View<ValueType**, Kokkos::HostSpace> stdDataView("stdDataView",
                                                           numTeams, numCols);
  // ensure that we use the same data to run the std algo on
  for (std::size_t i = 0; i < dataView.extent(0); ++i) {
    for (std::size_t j = 0; j < dataView.extent(1); ++j) {
      stdDataView(i, j) = cloneOfDataViewBeforeOp_h(i, j);
    }
  }
  GreaterThanValueFunctor predicate(threshold);
  for (std::size_t i = 0; i < dataView.extent(0); ++i) {
    auto thisRow = Kokkos::subview(stdDataView, i, Kokkos::ALL());
    std::replace_if(KE::begin(thisRow), KE::end(thisRow), predicate, newVal);
  }

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  auto dataViewAfterOp_h = create_host_space_copy(dataView);
  expect_equal_host_views(stdDataView, dataViewAfterOp_h);
}

template <class LayoutTag, class ValueType>
void run_all_scenarios() {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 8153}) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_replace_if_team_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamReplaceIf
}  // namespace stdalgos
}  // namespace Test
