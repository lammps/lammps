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

#ifndef KOKKOS_STD_ALGORITHMS_SEARCH_N_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_SEARCH_N_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include "Kokkos_AllOfAnyOfNoneOf.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class IteratorType, class SizeType, class ValueType,
          class ReducerType, class PredicateType>
struct StdSearchNFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType m_first;
  IteratorType m_last;
  SizeType m_count;
  ValueType m_value;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    namespace KE = ::Kokkos::Experimental;
    auto myit    = m_first + i;
    bool found   = true;

    for (SizeType k = 0; k < m_count; ++k) {
      // note that we add this EXPECT to check if we are in a valid range
      // but I think we can remove this beceause the guarantee we don't go
      // out of bounds is taken care of at the calling site
      // where we launch the par-reduce.
      KOKKOS_EXPECTS((myit + k) < m_last);

      if (!m_p(myit[k], m_value)) {
        found = false;
        break;
      }
    }

    // FIXME_NVHPC using a ternary operator causes problems
    red_value_type rv = {::Kokkos::reduction_identity<IndexType>::min()};
    if (found) {
      rv.min_loc_true = i;
    }

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdSearchNFunctor(IteratorType first, IteratorType last, SizeType count,
                    ValueType value, ReducerType reducer, PredicateType p)
      : m_first(std::move(first)),
        m_last(std::move(last)),
        m_count(std::move(count)),
        m_value(std::move(value)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType, class BinaryPredicateType>
IteratorType search_n_impl(const std::string& label, const ExecutionSpace& ex,
                           IteratorType first, IteratorType last,
                           SizeType count, const ValueType& value,
                           const BinaryPredicateType& pred) {
  // checks
  static_assert_random_access_and_accessible(ex, first);
  expect_valid_range(first, last);
  KOKKOS_EXPECTS((std::ptrdiff_t)count >= 0);

  // count should not be larger than the range [first, last)
  namespace KE            = ::Kokkos::Experimental;
  const auto num_elements = KE::distance(first, last);
  // cast things to avoid compiler warning
  KOKKOS_EXPECTS((std::size_t)num_elements >= (std::size_t)count);

  if (first == last) {
    return first;
  }

  // special case where num elements in [first, last) == count
  if ((std::size_t)num_elements == (std::size_t)count) {
    using equal_to_value = StdAlgoEqualsValUnaryPredicate<ValueType>;
    const auto satisfies =
        all_of_impl(label, ex, first, last, equal_to_value(value));
    return (satisfies) ? first : last;
  } else {
    // aliases
    using index_type           = typename IteratorType::difference_type;
    using reducer_type         = FirstLoc<index_type>;
    using reduction_value_type = typename reducer_type::value_type;
    using func_t =
        StdSearchNFunctor<index_type, IteratorType, SizeType, ValueType,
                          reducer_type, BinaryPredicateType>;

    // run
    reduction_value_type red_result;
    reducer_type reducer(red_result);

    // decide the size of the range policy of the par_red:
    // the last feasible index to start looking is the index
    // whose distance from the "last" is equal to count.
    // the +1 is because we need to include that location too.
    const auto range_size = num_elements - count + 1;

    // run par reduce
    ::Kokkos::parallel_reduce(
        label, RangePolicy<ExecutionSpace>(ex, 0, range_size),
        func_t(first, last, count, value, reducer, pred), reducer);

    // fence not needed because reducing into scalar

    // decide and return
    if (red_result.min_loc_true ==
        ::Kokkos::reduction_identity<index_type>::min()) {
      // location has not been found
      return last;
    } else {
      // location has been found
      return first + red_result.min_loc_true;
    }
  }
}

template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType>
IteratorType search_n_impl(const std::string& label, const ExecutionSpace& ex,
                           IteratorType first, IteratorType last,
                           SizeType count, const ValueType& value) {
  using iter_value_type = typename IteratorType::value_type;
  using predicate_type =
      StdAlgoEqualBinaryPredicate<iter_value_type, ValueType>;

  /* above we use <iter_value_type, ValueType> for the predicate_type
     to be consistent with the standard, which says:

     "
     The signature of the predicate function should be equivalent to:

        bool pred(const Type1 &a, const Type2 &b);

     The type Type1 must be such that an object of type ForwardIt can be
     dereferenced and then implicitly converted to Type1. The type Type2 must be
     such that an object of type T can be implicitly converted to Type2.
     "

     In our case, IteratorType = ForwardIt, and ValueType = T.
   */

  return search_n_impl(label, ex, first, last, count, value, predicate_type());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
