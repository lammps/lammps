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

namespace KE = Kokkos::Experimental;

namespace Test {
namespace stdalgos {

struct std_algorithms_min_max_element_test : std_algorithms_test {
  const int m_number_of_filling_cases = 5;

  void fillFixtureViews(int case_number) {
    static_view_t tmpView("tmpView");
    auto tmp_view_h = Kokkos::create_mirror_view(Kokkos::HostSpace(), tmpView);
    if (case_number == 1) {
      tmp_view_h(0) = 0;
      tmp_view_h(1) = 0;
      tmp_view_h(2) = 0;
      tmp_view_h(3) = 2;
      tmp_view_h(4) = 2;
      tmp_view_h(5) = 1;
      tmp_view_h(6) = 1;
      tmp_view_h(7) = 1;
      tmp_view_h(8) = 1;
      tmp_view_h(9) = 0;
    }

    else if (case_number == 2) {
      tmp_view_h(0) = 1;
      tmp_view_h(1) = 2;
      tmp_view_h(2) = 3;
      tmp_view_h(3) = 4;
      tmp_view_h(4) = 5;
      tmp_view_h(5) = 6;
      tmp_view_h(6) = 7;
      tmp_view_h(7) = 8;
      tmp_view_h(8) = 9;
      tmp_view_h(9) = 10;
    }

    else if (case_number == 3) {
      tmp_view_h(0) = 8;
      tmp_view_h(1) = 8;
      tmp_view_h(2) = -1;
      tmp_view_h(3) = -1;
      tmp_view_h(4) = 5;
      tmp_view_h(5) = 5;
      tmp_view_h(6) = 5;
      tmp_view_h(7) = 8;
      tmp_view_h(8) = 2;
      tmp_view_h(9) = 1;
    }

    else if (case_number == 4) {
      tmp_view_h(0) = 2;
      tmp_view_h(1) = 2;
      tmp_view_h(2) = 2;
      tmp_view_h(3) = 2;
      tmp_view_h(4) = 2;
      tmp_view_h(5) = 2;
      tmp_view_h(6) = 2;
      tmp_view_h(7) = 2;
      tmp_view_h(8) = 2;
      tmp_view_h(9) = 2;
    }

    else if (case_number == 5) {
      tmp_view_h(0) = 1;
      tmp_view_h(1) = 2;
      tmp_view_h(2) = 3;
      tmp_view_h(3) = 4;
      tmp_view_h(4) = 5;
      tmp_view_h(5) = 12;
      tmp_view_h(6) = 5;
      tmp_view_h(7) = 4;
      tmp_view_h(8) = 3;
      tmp_view_h(9) = 2;
    }

    else {
    }

    Kokkos::deep_copy(tmpView, tmp_view_h);
    copyInputViewToFixtureViews(tmpView);
  }

  Kokkos::pair<int, value_type> goldSolutionMaxElement(int caseNumber) {
    // returns {indexOfMaxElem, maxValue}
    if (caseNumber == 1) {
      return {3, 2};
    } else if (caseNumber == 2) {
      return {9, 10};
    } else if (caseNumber == 3) {
      return {0, 8};
    } else if (caseNumber == 4) {
      return {0, 2};
    } else if (caseNumber == 5) {
      return {5, 12};
    } else {
      return {};
    }
  }

  Kokkos::pair<int, value_type> goldSolutionMinElement(int caseNumber) {
    // returns {indexOfMinElem, minValue}
    if (caseNumber == 1) {
      return {0, 0};
    } else if (caseNumber == 2) {
      return {0, 1};
    } else if (caseNumber == 3) {
      return {2, -1};
    } else if (caseNumber == 4) {
      return {0, 2};
    } else if (caseNumber == 5) {
      return {0, 1};
    } else {
      return {};
    }
  }

  Kokkos::pair<Kokkos::pair<int, value_type>, Kokkos::pair<int, value_type>>
  goldSolutionMinMaxElement(int caseNumber) {
    // returns {{indexOfMinElem, minValue}, {indexOfMaxElem, maxValue}}
    // remember that for min it finds the first smallest element
    // remember that for max it finds the last biggest element
    if (caseNumber == 1) {
      return {{0, 0}, {4, 2}};
    } else if (caseNumber == 2) {
      return {{0, 1}, {9, 10}};
    } else if (caseNumber == 3) {
      return {{2, -1}, {7, 8}};
    } else if (caseNumber == 4) {
      return {{0, 2}, {9, 2}};
    } else if (caseNumber == 5) {
      return {{0, 1}, {5, 12}};
    } else {
      return {};
    }
  }

  template <class ViewType>
  void test_max_element_non_trivial_data(ViewType view);
  template <class ViewType>
  void test_min_element_non_trivial_data(ViewType view);
  template <class ViewType>
  void test_minmax_element_non_trivial_data(ViewType view);

  template <class ViewType>
  void test_max_element_non_trivial_data_custom_comp(ViewType view);
  template <class ViewType>
  void test_min_element_non_trivial_data_custom_comp(ViewType view);
  template <class ViewType>
  void test_minmax_element_non_trivial_data_custom_comp(ViewType view);
};

template <class IndexType, class ValueType, class ItType, class TestedViewType>
void std_algo_min_max_test_verify(Kokkos::pair<IndexType, ValueType> goldPair,
                                  const ItType result,
                                  TestedViewType testedView) {
  // check that iterator is pointing to right element
  EXPECT_EQ(result - KE::begin(testedView), goldPair.first);

  // create a view for the result to copy into it the iterator's value
  using result_view_t = Kokkos::View<int>;
  result_view_t resultView("result");
  CopyFromIteratorFunctor<ItType, result_view_t> cf(result, resultView);
  Kokkos::parallel_for("_std_algo_copy", 1, cf);
  auto result_v_h =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), resultView);

  // use the host mirror of the result view to check that the values match
  EXPECT_EQ(result_v_h(), goldPair.second);
}

template <class GoldSolutionType, class ItType, class TestedViewType>
void std_algo_min_max_test_verify(const GoldSolutionType& goldSolution,
                                  const ItType itMin, const ItType itMax,
                                  TestedViewType testedView) {
  std_algo_min_max_test_verify(goldSolution.first, itMin, testedView);
  std_algo_min_max_test_verify(goldSolution.second, itMax, testedView);
}

template <class ViewType>
void test_max_element_trivial_data(ViewType view) {
  /* if we pass empty range, should return last */
  auto result = KE::max_element(exespace(), KE::cbegin(view), KE::cbegin(view));
  EXPECT_EQ(result, KE::cbegin(view));

  /* if we pass empty range, should return last */
  auto it0     = KE::cbegin(view) + 3;
  auto it1     = it0;
  auto result2 = KE::max_element(exespace(), it0, it1);
  EXPECT_EQ(result2, it1);
}

template <class ViewType>
void test_min_element_trivial_data(ViewType view) {
  /* if we pass empty range, should return last */
  auto result = KE::min_element(exespace(), KE::cbegin(view), KE::cbegin(view));
  EXPECT_EQ(result, KE::cbegin(view));

  /* if we pass empty range, should return last */
  auto it0     = KE::cbegin(view) + 3;
  auto it1     = it0;
  auto result2 = KE::min_element(exespace(), it0, it1);
  EXPECT_EQ(result2, it1);
}

template <class ViewType>
void test_minmax_element_empty_range(ViewType view) {
  auto result =
      KE::minmax_element(exespace(), KE::cbegin(view), KE::cbegin(view));
  EXPECT_EQ(result.first, KE::cbegin(view));
  EXPECT_EQ(result.second, KE::cbegin(view));
  auto it0     = KE::cbegin(view) + 3;
  auto it1     = it0;
  auto result2 = KE::minmax_element(exespace(), it0, it1);
  EXPECT_EQ(result2.first, it1);
  EXPECT_EQ(result2.second, it1);
}

template <class ViewType>
void std_algorithms_min_max_element_test::test_max_element_non_trivial_data(
    ViewType view) {
  for (int id = 1; id <= m_number_of_filling_cases; ++id) {
    fillFixtureViews(id);
    const auto gold_solution = goldSolutionMaxElement(id);
    // API accepting view
    {
      const auto result = KE::max_element(exespace(), view);
      std_algo_min_max_test_verify(gold_solution, result, view);
      const auto result2 = KE::max_element("MYCUSTOMLABEL1", exespace(), view);
      std_algo_min_max_test_verify(gold_solution, result2, view);
    }
    // API accepting iterators
    {
      const auto result =
          KE::max_element(exespace(), KE::begin(view), KE::end(view));
      std_algo_min_max_test_verify(gold_solution, result, view);
      const auto result2 = KE::max_element("MYCUSTOMLABEL2", exespace(),
                                           KE::begin(view), KE::end(view));
      std_algo_min_max_test_verify(gold_solution, result2, view);
    }
  }
}

template <class ViewType>
void std_algorithms_min_max_element_test::test_min_element_non_trivial_data(
    ViewType view) {
  for (int id = 1; id <= m_number_of_filling_cases; ++id) {
    fillFixtureViews(id);
    const auto goldPair = goldSolutionMinElement(id);
    // API accepting view
    {
      const auto result = KE::min_element(exespace(), view);
      std_algo_min_max_test_verify(goldPair, result, view);
      const auto result2 = KE::min_element("MYCUSTOMLABEL1", exespace(), view);
      std_algo_min_max_test_verify(goldPair, result2, view);
    }
    // API accepting iterators
    {
      const auto result =
          KE::min_element(exespace(), KE::begin(view), KE::end(view));
      std_algo_min_max_test_verify(goldPair, result, view);
      const auto result2 = KE::min_element("MYCUSTOMLABEL2", exespace(),
                                           KE::begin(view), KE::end(view));
      std_algo_min_max_test_verify(goldPair, result2, view);
    }
  }
}

template <class ViewType>
void std_algorithms_min_max_element_test::test_minmax_element_non_trivial_data(
    ViewType view) {
  for (int id = 1; id <= m_number_of_filling_cases; ++id) {
    fillFixtureViews(id);
    const auto gold = goldSolutionMinMaxElement(id);
    {
      auto result = KE::minmax_element(exespace(), view);
      std_algo_min_max_test_verify(gold, result.first, result.second, view);
      const auto result2 =
          KE::minmax_element("MYCUSTOMLABEL1", exespace(), view);
      std_algo_min_max_test_verify(gold, result2.first, result2.second, view);
    }
    {
      const auto result =
          KE::minmax_element(exespace(), KE::begin(view), KE::end(view));
      std_algo_min_max_test_verify(gold, result.first, result.second, view);
      const auto result2 = KE::minmax_element("MYCUSTOMLABEL2", exespace(),
                                              KE::begin(view), KE::end(view));
      std_algo_min_max_test_verify(gold, result2.first, result2.second, view);
    }
  }
}

#if !defined KOKKOS_ENABLE_OPENMPTARGET
template <class ViewType>
void std_algorithms_min_max_element_test::
    test_max_element_non_trivial_data_custom_comp(ViewType view) {
  for (int id = 1; id <= m_number_of_filling_cases; ++id) {
    fillFixtureViews(id);
    const auto goldPair = goldSolutionMaxElement(id);
    CustomLessThanComparator<value_type, value_type> comp;
    // API accepting view
    {
      const auto result = KE::max_element(exespace(), view, comp);
      std_algo_min_max_test_verify(goldPair, result, view);
      const auto result2 =
          KE::max_element("MYCUSTOMLABEL3", exespace(), view, comp);
      std_algo_min_max_test_verify(goldPair, result2, view);
    }
    // API accepting iterators
    {
      const auto result =
          KE::max_element(exespace(), KE::begin(view), KE::end(view), comp);
      std_algo_min_max_test_verify(goldPair, result, view);
      const auto result2 = KE::max_element(
          "MYCUSTOMLABEL4", exespace(), KE::begin(view), KE::end(view), comp);
      std_algo_min_max_test_verify(goldPair, result2, view);
    }
  }
}

template <class ViewType>
void std_algorithms_min_max_element_test::
    test_min_element_non_trivial_data_custom_comp(ViewType view) {
  for (int id = 1; id <= m_number_of_filling_cases; ++id) {
    fillFixtureViews(id);
    const auto goldPair = goldSolutionMinElement(id);
    CustomLessThanComparator<value_type, value_type> comp;
    // API accepting view
    {
      const auto result = KE::min_element(exespace(), view, comp);
      std_algo_min_max_test_verify(goldPair, result, view);
      const auto result2 =
          KE::min_element("MYCUSTOMLABEL3", exespace(), view, comp);
      std_algo_min_max_test_verify(goldPair, result2, view);
    }
    // API accepting iterators
    {
      const auto result =
          KE::min_element(exespace(), KE::begin(view), KE::end(view), comp);
      std_algo_min_max_test_verify(goldPair, result, view);
      const auto result2 = KE::min_element(
          "MYCUSTOMLABEL4", exespace(), KE::begin(view), KE::end(view), comp);
      std_algo_min_max_test_verify(goldPair, result2, view);
    }
  }
}

template <class ViewType>
void std_algorithms_min_max_element_test::
    test_minmax_element_non_trivial_data_custom_comp(ViewType view) {
  for (int id = 1; id <= m_number_of_filling_cases; ++id) {
    fillFixtureViews(id);
    const auto goldPair = goldSolutionMinMaxElement(id);
    CustomLessThanComparator<value_type, value_type> comp;
    {
      const auto result = KE::minmax_element(exespace(), view, comp);
      std_algo_min_max_test_verify(goldPair, result.first, result.second, view);
      const auto result2 =
          KE::minmax_element("MYCUSTOMLABEL3", exespace(), view, comp);
      std_algo_min_max_test_verify(goldPair, result2.first, result2.second,
                                   view);
    }
    {
      const auto result =
          KE::minmax_element(exespace(), KE::begin(view), KE::end(view), comp);
      std_algo_min_max_test_verify(goldPair, result.first, result.second, view);
      const auto result2 = KE::minmax_element(
          "MYCUSTOMLABEL4", exespace(), KE::begin(view), KE::end(view), comp);
      std_algo_min_max_test_verify(goldPair, result2.first, result2.second,
                                   view);
    }
  }
}
#endif

// trivial case
TEST_F(std_algorithms_min_max_element_test, min_element_empty_range) {
  test_min_element_trivial_data(m_static_view);
  test_min_element_trivial_data(m_dynamic_view);
  test_min_element_trivial_data(m_strided_view);
}

TEST_F(std_algorithms_min_max_element_test, max_element_empty_range) {
  test_max_element_trivial_data(m_static_view);
  test_max_element_trivial_data(m_dynamic_view);
  test_max_element_trivial_data(m_strided_view);
}

// non-trivial data
TEST_F(std_algorithms_min_max_element_test, min_element_non_trivial_data) {
  test_min_element_non_trivial_data(m_static_view);
  test_min_element_non_trivial_data(m_dynamic_view);
  test_min_element_non_trivial_data(m_strided_view);
}

TEST_F(std_algorithms_min_max_element_test, max_element_non_trivial_data) {
  test_max_element_non_trivial_data(m_static_view);
  test_max_element_non_trivial_data(m_dynamic_view);
  test_max_element_non_trivial_data(m_strided_view);
}

#if !defined KOKKOS_ENABLE_OPENMPTARGET
// non-trivial data, custom comp
TEST_F(std_algorithms_min_max_element_test,
       min_element_non_trivial_data_custom_comp) {
  test_min_element_non_trivial_data_custom_comp(m_static_view);
  test_min_element_non_trivial_data_custom_comp(m_dynamic_view);
  test_min_element_non_trivial_data_custom_comp(m_strided_view);
}

TEST_F(std_algorithms_min_max_element_test,
       max_element_non_trivial_data_custom_comp) {
  test_max_element_non_trivial_data_custom_comp(m_static_view);
  test_max_element_non_trivial_data_custom_comp(m_dynamic_view);
  test_max_element_non_trivial_data_custom_comp(m_strided_view);
}
#endif

#if defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(KOKKOS_COMPILER_CLANG) && \
    (KOKKOS_COMPILER_CLANG >= 1300)
TEST_F(std_algorithms_min_max_element_test, minmax_element_empty_range) {
  test_minmax_element_empty_range(m_static_view);
  test_minmax_element_empty_range(m_dynamic_view);
  test_minmax_element_empty_range(m_strided_view);
}

TEST_F(std_algorithms_min_max_element_test, minmax_element_non_trivial_data) {
  test_minmax_element_non_trivial_data(m_static_view);
  test_minmax_element_non_trivial_data(m_dynamic_view);
  test_minmax_element_non_trivial_data(m_strided_view);
}
#endif

#if !defined KOKKOS_ENABLE_OPENMPTARGET
// OpenMPTarget does not yet support custom comparator
TEST_F(std_algorithms_min_max_element_test,
       minmax_element_non_trivial_data_custom_comp) {
  test_minmax_element_non_trivial_data_custom_comp(m_static_view);
  test_minmax_element_non_trivial_data_custom_comp(m_dynamic_view);
  test_minmax_element_non_trivial_data_custom_comp(m_strided_view);
}
#endif

}  // namespace stdalgos
}  // namespace Test
