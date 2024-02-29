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
namespace TeamGenerate {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct Generator {
  KOKKOS_INLINE_FUNCTION
  ValueType operator()() const { return static_cast<ValueType>(23); }
};

template <class ViewType>
struct TestFunctorA {
  ViewType m_view;
  int m_apiPick;

  TestFunctorA(const ViewType view, int apiPick)
      : m_view(view), m_apiPick(apiPick) {}

  template <class MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType& member) const {
    const auto myRowIndex = member.league_rank();
    auto myRowView        = Kokkos::subview(m_view, myRowIndex, Kokkos::ALL());

    using value_type = typename ViewType::value_type;
    if (m_apiPick == 0) {
      KE::generate(member, KE::begin(myRowView), KE::end(myRowView),
                   Generator<value_type>());
    } else if (m_apiPick == 1) {
      KE::generate(member, myRowView, Generator<value_type>());
    }
  }
};

template <class LayoutTag, class ValueType>
void test_A(std::size_t numTeams, std::size_t numCols, int apiId) {
  /* description:
     randomly fill a view and then do a team-level generate
     with one team per row to assign to each element a value
     produced via a generator functor
   */

  // -----------------------------------------------
  // prepare data
  // -----------------------------------------------
  // create a view in the memory space associated with default exespace
  // with as many rows as the number of teams and fill it with random
  // values from an arbitrary range. Pick range so that it does NOT
  // contain the value produced by the generator (see top of file)
  // otherwise test check below is ill-posed
  auto [dataView, cloneOfDataViewBeforeOp_h] =
      create_random_view_and_host_clone(
          LayoutTag{}, numTeams, numCols,
          Kokkos::pair<ValueType, ValueType>{105, 523}, "dataView");

  // -----------------------------------------------
  // launch kokkos kernel
  // -----------------------------------------------
  using space_t = Kokkos::DefaultExecutionSpace;
  Kokkos::TeamPolicy<space_t> policy(numTeams, Kokkos::AUTO());

  // use CTAD for functor
  TestFunctorA fnc(dataView, apiId);
  Kokkos::parallel_for(policy, fnc);

  // -----------------------------------------------
  // check
  // -----------------------------------------------
  auto dataViewAfterOp_h = create_host_space_copy(dataView);
  for (std::size_t i = 0; i < dataViewAfterOp_h.extent(0); ++i) {
    for (std::size_t j = 0; j < dataViewAfterOp_h.extent(1); ++j) {
      EXPECT_TRUE(dataViewAfterOp_h(i, j) == static_cast<ValueType>(23));
      EXPECT_TRUE(dataViewAfterOp_h(i, j) != cloneOfDataViewBeforeOp_h(i, j));
    }
  }
}

template <class LayoutTag, class ValueType>
void run_all_scenarios() {
  for (int numTeams : teamSizesToTest) {
    for (const auto& numCols : {0, 1, 2, 13, 101, 1444, 51153}) {
      for (int apiId : {0, 1}) {
        test_A<LayoutTag, ValueType>(numTeams, numCols, apiId);
      }
    }
  }
}

TEST(std_algorithms_generate_team_test, test_unary_op) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoRowsTag, int>();
  run_all_scenarios<StridedThreeRowsTag, unsigned>();
}

}  // namespace TeamGenerate
}  // namespace stdalgos
}  // namespace Test
