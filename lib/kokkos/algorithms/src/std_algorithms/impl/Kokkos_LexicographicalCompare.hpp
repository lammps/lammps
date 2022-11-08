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

#ifndef KOKKOS_STD_ALGORITHMS_LEXICOGRAPHICAL_COMPARE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_LEXICOGRAPHICAL_COMPARE_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class IteratorType1, class IteratorType2,
          class ComparatorType>
struct StdCompareFunctor {
  IteratorType1 m_it1;
  IteratorType2 m_it2;
  ComparatorType m_predicate;

  KOKKOS_FUNCTION
  void operator()(IndexType /* i is unused */, int& lsum) const {
    if (m_predicate(*m_it1, *m_it2)) {
      lsum = 1;
    }
  }

  KOKKOS_FUNCTION
  StdCompareFunctor(IteratorType1 _it1, IteratorType2 _it2,
                    ComparatorType _predicate)
      : m_it1(std::move(_it1)),
        m_it2(std::move(_it2)),
        m_predicate(std::move(_predicate)) {}
};

template <class IndexType, class IteratorType1, class IteratorType2,
          class ReducerType, class ComparatorType>
struct StdLexicographicalCompareFunctor {
  using red_value_type = typename ReducerType::value_type;
  IteratorType1 m_first1;
  IteratorType2 m_first2;
  ReducerType m_reducer;
  ComparatorType m_comparator;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    const auto& my_value1 = m_first1[i];
    const auto& my_value2 = m_first2[i];

    bool different = m_comparator(my_value1, my_value2) ||
                     m_comparator(my_value2, my_value1);
    auto rv =
        different
            ? red_value_type{i}
            : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdLexicographicalCompareFunctor(IteratorType1 _first1, IteratorType2 _first2,
                                   ReducerType _reducer, ComparatorType _comp)
      : m_first1(std::move(_first1)),
        m_first2(std::move(_first2)),
        m_reducer(std::move(_reducer)),
        m_comparator(std::move(_comp)) {}
};

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ComparatorType>
bool lexicographical_compare_impl(const std::string& label,
                                  const ExecutionSpace& ex,
                                  IteratorType1 first1, IteratorType1 last1,
                                  IteratorType2 first2, IteratorType2 last2,
                                  ComparatorType comp) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);
  Impl::expect_valid_range(first2, last2);

  // aliases
  using index_type           = typename IteratorType1::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;

  // run
  const auto d1    = Kokkos::Experimental::distance(first1, last1);
  const auto d2    = Kokkos::Experimental::distance(first2, last2);
  const auto range = Kokkos::min(d1, d2);
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  using func1_t =
      StdLexicographicalCompareFunctor<index_type, IteratorType1, IteratorType2,
                                       reducer_type, ComparatorType>;

  ::Kokkos::parallel_reduce(label, RangePolicy<ExecutionSpace>(ex, 0, range),
                            func1_t(first1, first2, reducer, comp), reducer);

  // fence not needed because reducing into scalar
  // no mismatch
  if (red_result.min_loc_true ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    auto new_last1 = first1 + range;
    auto new_last2 = first2 + range;
    bool is_prefix = (new_last1 == last1) && (new_last2 != last2);
    return is_prefix;
  }

  // check mismatched
  int less      = 0;
  auto it1      = first1 + red_result.min_loc_true;
  auto it2      = first2 + red_result.min_loc_true;
  using func2_t = StdCompareFunctor<index_type, IteratorType1, IteratorType2,
                                    ComparatorType>;
  ::Kokkos::parallel_reduce(label, RangePolicy<ExecutionSpace>(ex, 0, 1),
                            func2_t(it1, it2, comp), less);

  // fence not needed because reducing into scalar
  return static_cast<bool>(less);
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
bool lexicographical_compare_impl(const std::string& label,
                                  const ExecutionSpace& ex,
                                  IteratorType1 first1, IteratorType1 last1,
                                  IteratorType2 first2, IteratorType2 last2) {
  using value_type_1 = typename IteratorType1::value_type;
  using value_type_2 = typename IteratorType2::value_type;
  using predicate_t =
      Impl::StdAlgoLessThanBinaryPredicate<value_type_1, value_type_2>;
  return lexicographical_compare_impl(label, ex, first1, last1, first2, last2,
                                      predicate_t());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
