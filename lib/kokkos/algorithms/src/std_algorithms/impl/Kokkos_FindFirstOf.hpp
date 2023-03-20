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

#ifndef KOKKOS_STD_ALGORITHMS_FIND_FIRST_OF_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_FIND_FIRST_OF_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class IteratorType1, class IteratorType2,
          class ReducerType, class PredicateType>
struct StdFindFirstOfFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType1 m_first;
  IteratorType2 m_s_first;
  IteratorType2 m_s_last;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    namespace KE        = ::Kokkos::Experimental;
    const auto& myvalue = m_first[i];
    bool found          = false;

    const auto search_count = KE::distance(m_s_first, m_s_last);
    for (IndexType k = 0; k < search_count; ++k) {
      if (m_p(myvalue, m_s_first[k])) {
        found = true;
        break;
      }
    }

    const auto rv =
        found ? red_value_type{i}
              : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdFindFirstOfFunctor(IteratorType1 first, IteratorType2 s_first,
                        IteratorType2 s_last, ReducerType reducer,
                        PredicateType p)
      : m_first(std::move(first)),
        m_s_first(std::move(s_first)),
        m_s_last(std::move(s_last)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
IteratorType1 find_first_of_impl(const std::string& label,
                                 const ExecutionSpace& ex, IteratorType1 first,
                                 IteratorType1 last, IteratorType2 s_first,
                                 IteratorType2 s_last,
                                 const BinaryPredicateType& pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, s_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, s_first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(s_first, s_last);

  if ((s_first == s_last) || (first == last)) {
    return last;
  }

  using index_type           = typename IteratorType1::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t = StdFindFirstOfFunctor<index_type, IteratorType1, IteratorType2,
                                       reducer_type, BinaryPredicateType>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      func_t(first, s_first, s_last, reducer, pred), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  if (red_result.min_loc_true ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    // if here, nothing found
    return last;
  } else {
    // a location has been found
    return first + red_result.min_loc_true;
  }
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType1 find_first_of_impl(const std::string& label,
                                 const ExecutionSpace& ex, IteratorType1 first,
                                 IteratorType1 last, IteratorType2 s_first,
                                 IteratorType2 s_last) {
  using value_type1    = typename IteratorType1::value_type;
  using value_type2    = typename IteratorType2::value_type;
  using predicate_type = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return find_first_of_impl(label, ex, first, last, s_first, s_last,
                            predicate_type());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
