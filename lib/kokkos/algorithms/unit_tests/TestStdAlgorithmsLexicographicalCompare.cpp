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
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace LexicographicalCompare {

namespace KE = Kokkos::Experimental;

template <class ViewType1, class ViewType2>
void test_lexicographical_compare(const ViewType1 view_1, ViewType2 view_2) {
  auto host_copy_1 = create_host_space_copy(view_1);
  auto host_copy_2 = create_host_space_copy(view_2);

  auto first_1 = KE::begin(view_1);
  auto last_1  = KE::end(view_1);
  auto first_2 = KE::begin(view_2);
  auto last_2  = KE::end(view_2);

  auto h_first_1 = KE::begin(host_copy_1);
  auto h_last_1  = KE::end(host_copy_1);
  auto h_first_2 = KE::begin(host_copy_2);
  auto h_last_2  = KE::end(host_copy_2);

  {
    // default comparator
    auto std_result =
        std::lexicographical_compare(h_first_1, h_last_1, h_first_2, h_last_2);

    // pass iterators
    EXPECT_EQ(std_result, KE::lexicographical_compare(exespace(), first_1,
                                                      last_1, first_2, last_2));
    EXPECT_EQ(std_result,
              KE::lexicographical_compare("label", exespace(), first_1, last_1,
                                          first_2, last_2));

    // pass views
    EXPECT_EQ(std_result,
              KE::lexicographical_compare(exespace(), view_1, view_2));
    EXPECT_EQ(std_result,
              KE::lexicographical_compare("label", exespace(), view_1, view_2));
  }

  {
    // custom comparator
    using value_t_1 = typename ViewType1::value_type;
    using value_t_2 = typename ViewType2::value_type;
    const auto custom_comparator =
        CustomLessThanComparator<value_t_1, value_t_2>();
    auto std_result = std::lexicographical_compare(
        h_first_1, h_last_1, h_first_2, h_last_2, custom_comparator);

    // pass iterators
    EXPECT_EQ(std_result,
              KE::lexicographical_compare(exespace(), first_1, last_1, first_2,
                                          last_2, custom_comparator));
    EXPECT_EQ(std_result,
              KE::lexicographical_compare("label", exespace(), first_1, last_1,
                                          first_2, last_2, custom_comparator));

    // pass views
    EXPECT_EQ(std_result, KE::lexicographical_compare(
                              exespace(), view_1, view_2, custom_comparator));
    EXPECT_EQ(std_result,
              KE::lexicographical_compare("label", exespace(), view_1, view_2,
                                          custom_comparator));
  }

  {
    // empty vs non-empty
    auto std_result =
        std::lexicographical_compare(h_first_1, h_first_1, h_first_2, h_last_2);
    EXPECT_EQ(std_result, KE::lexicographical_compare(
                              exespace(), first_1, first_1, first_2, last_2));
  }

  {
    // pass shorter range
    if (view_1.extent(0) > 1) {
      auto std_result = std::lexicographical_compare(h_first_1, h_last_1 - 1,
                                                     h_first_2, h_last_2);
      EXPECT_EQ(std_result,
                KE::lexicographical_compare(exespace(), first_1, last_1 - 1,
                                            first_2, last_2));
    }
  }

  {
    // first element smaller
    if (view_1.extent(0) > 1) {
      KE::fill(exespace(), first_1, first_1 + 1, 1);
      KE::fill(exespace(), first_2, first_2 + 1, 2);

      EXPECT_TRUE(KE::lexicographical_compare(exespace(), first_1, last_1,
                                              first_2, last_2));
    }
  }

  {
    // first element bigger, last element smaller
    if (view_1.extent(0) > 2) {
      KE::fill(exespace(), first_1, first_1 + 1, 2);
      KE::fill(exespace(), first_2, first_2 + 1, 1);

      KE::fill(exespace(), last_1 - 1, last_1, 1);
      KE::fill(exespace(), last_2 - 1, last_2, 2);

      EXPECT_FALSE(KE::lexicographical_compare(exespace(), first_1, last_1,
                                               first_2, last_2));
    }
  }
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  for (const auto& scenario : default_scenarios) {
    auto view1 = create_view<ValueType>(Tag{}, scenario.second,
                                        "lexicographical_compare_1");
    auto view2 = create_view<ValueType>(Tag{}, scenario.second,
                                        "lexicographical_compare_2");

    test_lexicographical_compare(view1, view2);
  }
}

TEST(std_algorithms_lexicographical_compare_test, test) {
// FIXME: should this disable only custom comparator tests?
#if !defined KOKKOS_ENABLE_OPENMPTARGET
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoTag, int>();
  run_all_scenarios<StridedThreeTag, unsigned>();
#endif
}

}  // namespace LexicographicalCompare
}  // namespace stdalgos
}  // namespace Test
