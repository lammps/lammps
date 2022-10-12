/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

// purpose of this test is to check that the reducers used
// to implement some std algorithms work independently of the order

namespace Test {

enum class StdReducersTestEnumOrder { LeftToRight, RightToLeft, Random };

std::string order_to_string(StdReducersTestEnumOrder value) {
  switch (value) {
    case StdReducersTestEnumOrder::LeftToRight: return "LeftToRight";
    case StdReducersTestEnumOrder::RightToLeft: return "RightToLeft";
    case StdReducersTestEnumOrder::Random: return "Random";
  }
  return {};
}

auto create_host_view_with_reduction_order_indices(
    std::size_t extent, StdReducersTestEnumOrder enum_value) {
  using view_t = Kokkos::View<int*, Kokkos::HostSpace>;
  view_t result("v", extent);

  if (enum_value == StdReducersTestEnumOrder::LeftToRight) {
    result(0) = 0;
    result(1) = 1;
    result(2) = 2;
    result(3) = 3;
    result(4) = 4;
    result(5) = 5;
    result(6) = 6;
    result(7) = 7;
    result(8) = 8;
    result(9) = 9;
  } else if (enum_value == StdReducersTestEnumOrder::RightToLeft) {
    result(0) = 9;
    result(1) = 8;
    result(2) = 7;
    result(3) = 6;
    result(4) = 5;
    result(5) = 4;
    result(6) = 3;
    result(7) = 2;
    result(8) = 1;
    result(9) = 0;
  } else if (enum_value == StdReducersTestEnumOrder::Random) {
    result(0) = 0;
    result(1) = 8;
    result(2) = 3;
    result(3) = 2;
    result(4) = 9;
    result(5) = 4;
    result(6) = 6;
    result(7) = 1;
    result(8) = 7;
    result(9) = 5;
  } else {
    throw std::runtime_error("test: Invalid enum");
  }

  return result;
}

template <int flag, class ExeSpace, class IndexType, class ViewType>
auto run_min_or_max_test(ViewType view, StdReducersTestEnumOrder enValue) {
  static_assert(std::is_same<ExeSpace, Kokkos::HostSpace>::value,
                "test is only enabled for HostSpace");

  std::cout << "checking reduction with order: " << order_to_string(enValue)
            << "\n";

  using view_value_type = typename ViewType::value_type;
  using reducer_type    = std::conditional_t<
      (flag == 0), Kokkos::MaxFirstLoc<view_value_type, IndexType, ExeSpace>,
      Kokkos::MinFirstLoc<view_value_type, IndexType, ExeSpace> >;
  using reduction_value_type = typename reducer_type::value_type;

  reduction_value_type red_result;
  reducer_type reducer(red_result);
  EXPECT_TRUE(reducer.references_scalar());
  reducer.init(red_result);

  auto red_order =
      create_host_view_with_reduction_order_indices(view.extent(0), enValue);
  for (std::size_t i = 0; i < view.extent(0); ++i) {
    const auto index = red_order(i);
    reducer.join(red_result, reduction_value_type{view(index), index});
  }

  using return_type = Kokkos::pair<view_value_type, IndexType>;
  return return_type{red_result.val, red_result.loc};
}

TEST(std_algorithms_reducers, max_first_loc) {
  using hostspace = Kokkos::HostSpace;

  using view_t                 = Kokkos::View<double*, hostspace>;
  constexpr std::size_t extent = 10;
  view_t view_h("v", extent);
  view_h(0) = 0.;
  view_h(1) = 0.;
  view_h(2) = 0.;
  view_h(3) = 2.;
  view_h(4) = 2.;
  view_h(5) = 1.;
  view_h(6) = 1.;
  view_h(7) = 1.;
  view_h(8) = 1.;
  view_h(9) = 0.;

  using index_type                 = int;
  using view_value_type            = typename view_t::value_type;
  const view_value_type gold_value = 2.;
  const index_type gold_location   = 3;

  const auto pair1 = run_min_or_max_test<0, hostspace, index_type>(
      view_h, StdReducersTestEnumOrder::LeftToRight);
  ASSERT_EQ(pair1.first, gold_value);
  ASSERT_EQ(pair1.second, gold_location);

  const auto pair2 = run_min_or_max_test<0, hostspace, index_type>(
      view_h, StdReducersTestEnumOrder::RightToLeft);
  ASSERT_EQ(pair2.first, gold_value);
  ASSERT_EQ(pair2.second, gold_location);

  const auto pair3 = run_min_or_max_test<0, hostspace, index_type>(
      view_h, StdReducersTestEnumOrder::Random);
  ASSERT_EQ(pair3.first, gold_value);
  ASSERT_EQ(pair3.second, gold_location);
}

TEST(std_algorithms_reducers, min_first_loc) {
  using hostspace = Kokkos::HostSpace;

  using view_t                 = Kokkos::View<double*, hostspace>;
  constexpr std::size_t extent = 10;
  view_t view_h("v", extent);
  view_h(0) = 0.;
  view_h(1) = 0.;
  view_h(2) = 0.;
  view_h(3) = 2.;
  view_h(4) = 2.;
  view_h(5) = -1.;
  view_h(6) = -1.;
  view_h(7) = 1.;
  view_h(8) = 1.;
  view_h(9) = 0.;

  using index_type                 = int;
  using view_value_type            = typename view_t::value_type;
  const view_value_type gold_value = -1.;
  const index_type gold_location   = 5;

  const auto pair1 = run_min_or_max_test<1, hostspace, index_type>(
      view_h, StdReducersTestEnumOrder::LeftToRight);
  ASSERT_EQ(pair1.first, gold_value);
  ASSERT_EQ(pair1.second, gold_location);

  const auto pair2 = run_min_or_max_test<1, hostspace, index_type>(
      view_h, StdReducersTestEnumOrder::RightToLeft);
  ASSERT_EQ(pair2.first, gold_value);
  ASSERT_EQ(pair2.second, gold_location);

  const auto pair3 = run_min_or_max_test<1, hostspace, index_type>(
      view_h, StdReducersTestEnumOrder::Random);
  ASSERT_EQ(pair3.first, gold_value);
  ASSERT_EQ(pair3.second, gold_location);
}

template <class ExeSpace, class IndexType, class ViewType, class ValuesPair,
          class IndexPair>
void run_min_max_test(ViewType view, StdReducersTestEnumOrder enValue,
                      const ValuesPair gold_values, const IndexPair gold_locs) {
  static_assert(std::is_same<ExeSpace, Kokkos::HostSpace>::value,
                "test is only enabled for HostSpace");

  std::cout << "checking reduction with order: " << order_to_string(enValue)
            << "\n";

  using view_value_type = typename ViewType::value_type;
  using reducer_type =
      Kokkos::MinMaxFirstLastLoc<view_value_type, IndexType, ExeSpace>;
  using reduction_value_type = typename reducer_type::value_type;

  reduction_value_type red_result;
  reducer_type reducer(red_result);
  EXPECT_TRUE(reducer.references_scalar());
  reducer.init(red_result);

  auto red_order =
      create_host_view_with_reduction_order_indices(view.extent(0), enValue);
  for (std::size_t i = 0; i < view.extent(0); ++i) {
    const auto index = red_order(i);
    reducer.join(red_result,
                 reduction_value_type{view(index), view(index), index, index});
  }

  ASSERT_EQ(red_result.min_val, gold_values.first);
  ASSERT_EQ(red_result.max_val, gold_values.second);
  ASSERT_EQ(red_result.min_loc, gold_locs.first);
  ASSERT_EQ(red_result.max_loc, gold_locs.second);
}

TEST(std_algorithms_reducers, min_max_first_last_loc) {
  using hostspace = Kokkos::HostSpace;

  using view_t                 = Kokkos::View<double*, hostspace>;
  constexpr std::size_t extent = 10;
  view_t view_h("v", extent);
  view_h(0) = 0.;
  view_h(1) = 0.;
  view_h(2) = 0.;
  view_h(3) = 2.;
  view_h(4) = 2.;
  view_h(5) = -1.;
  view_h(6) = 1.;
  view_h(7) = -1.;
  view_h(8) = 2.;
  view_h(9) = 0.;

  using index_type      = int;
  using view_value_type = typename view_t::value_type;
  Kokkos::pair<view_value_type, view_value_type> gold_values = {-1., 2.};
  Kokkos::pair<index_type, index_type> gold_indices          = {5, 8};

  run_min_max_test<hostspace, index_type>(
      view_h, StdReducersTestEnumOrder::LeftToRight, gold_values, gold_indices);

  run_min_max_test<hostspace, index_type>(
      view_h, StdReducersTestEnumOrder::RightToLeft, gold_values, gold_indices);

  run_min_max_test<hostspace, index_type>(
      view_h, StdReducersTestEnumOrder::Random, gold_values, gold_indices);
}

}  // namespace Test
