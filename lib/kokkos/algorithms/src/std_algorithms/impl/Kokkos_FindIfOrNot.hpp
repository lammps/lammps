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

#ifndef KOKKOS_STD_ALGORITHMS_FIND_IF_AND_FIND_IF_NOT_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_FIND_IF_AND_FIND_IF_NOT_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <bool is_find_if, class IndexType, class IteratorType,
          class ReducerType, class PredicateType>
struct StdFindIfOrNotFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType m_first;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    const auto& my_value = m_first[i];

    // if doing find_if, look for when predicate is true
    // if doing find_if_not, look for when predicate is false
    const bool found_condition = is_find_if ? m_p(my_value) : !m_p(my_value);

    auto rv =
        found_condition
            ? red_value_type{i}
            : red_value_type{::Kokkos::reduction_identity<IndexType>::min()};

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdFindIfOrNotFunctor(IteratorType first, ReducerType reducer,
                        PredicateType p)
      : m_first(std::move(first)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <bool is_find_if, class ExecutionSpace, class IteratorType,
          class PredicateType>
IteratorType find_if_or_not_impl(const std::string& label,
                                 const ExecutionSpace& ex, IteratorType first,
                                 IteratorType last, PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(
      ex, first);  // only need one It per type
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return last;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t = StdFindIfOrNotFunctor<is_find_if, index_type, IteratorType,
                                       reducer_type, PredicateType>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer, pred), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  if (red_result.min_loc_true ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    // here, it means a valid loc has not been found,
    return last;
  } else {
    // a location has been found
    return first + red_result.min_loc_true;
  }
}

template <class ExecutionSpace, class InputIterator, class T>
InputIterator find_impl(const std::string& label, ExecutionSpace ex,
                        InputIterator first, InputIterator last,
                        const T& value) {
  return find_if_or_not_impl<true>(
      label, ex, first, last,
      ::Kokkos::Experimental::Impl::StdAlgoEqualsValUnaryPredicate<T>(value));
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
