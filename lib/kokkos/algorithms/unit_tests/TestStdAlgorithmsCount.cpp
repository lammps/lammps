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

#include <TestStdAlgorithmsCommon.hpp>
#include <std_algorithms/Kokkos_BeginEnd.hpp>
#include <std_algorithms/Kokkos_NonModifyingSequenceOperations.hpp>
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace Count {

namespace KE = Kokkos::Experimental;

template <class ViewType>
void test_count(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;

  view_host_space_t expected("count_expected", view.extent(0));
  compare_views(expected, view);

  {
    const value_t count_value = 0;
    const auto std_result =
        std::count(KE::cbegin(expected), KE::cend(expected), count_value);
    EXPECT_EQ(view.extent(0), size_t(std_result));

    // pass const iterators
    EXPECT_EQ(std_result, KE::count(exespace(), KE::cbegin(view),
                                    KE::cend(view), count_value));
    // pass view
    EXPECT_EQ(std_result, KE::count(exespace(), view, count_value));
  }

  {
    const value_t count_value = 13;
    const auto std_result =
        std::count(KE::cbegin(expected), KE::cend(expected), count_value);

    // pass iterators
    EXPECT_EQ(std_result, KE::count("label", exespace(), KE::begin(view),
                                    KE::end(view), count_value));
    // pass view
    EXPECT_EQ(std_result, KE::count("label", exespace(), view, count_value));
  }
}

template <class ViewType>
void test_count_if(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;

  view_host_space_t expected("count_expected", view.extent(0));
  compare_views(expected, view);

  // no positive elements (all zeroes)
  const auto predicate = IsPositiveFunctor<value_type>();
  EXPECT_EQ(0,
            std::count_if(KE::begin(expected), KE::end(expected), predicate));

  // pass iterators
  EXPECT_EQ(
      0, KE::count_if(exespace(), KE::begin(view), KE::end(view), predicate));
  // pass view
  EXPECT_EQ(0, KE::count_if(exespace(), view, predicate));

  fill_views_inc(view, expected);

  const auto std_result =
      std::count_if(KE::begin(expected), KE::end(expected), predicate);
  // pass const iterators
  EXPECT_EQ(std_result, KE::count_if("label", exespace(), KE::cbegin(view),
                                     KE::cend(view), predicate));
  // pass view
  EXPECT_EQ(std_result, KE::count_if("label", exespace(), view, predicate));
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  for (const auto& scenario : default_scenarios) {
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "count");
      test_count(view);
    }
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "count");
      test_count_if(view);
    }
  }
}

TEST(std_algorithms_count_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoTag, int>();
  run_all_scenarios<StridedThreeTag, unsigned>();
}

}  // namespace Count
}  // namespace stdalgos
}  // namespace Test
