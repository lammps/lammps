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

#ifndef KOKKOS_STD_ALGORITHMS_MISMATCH_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_MISMATCH_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class IteratorType1, class IteratorType2,
          class ReducerType, class BinaryPredicateType>
struct StdMismatchRedFunctor {
  using red_value_type = typename ReducerType::value_type;

  IteratorType1 m_first1;
  IteratorType2 m_first2;
  ReducerType m_reducer;
  BinaryPredicateType m_predicate;

  KOKKOS_FUNCTION
  void operator()(const IndexType i, red_value_type& red_value) const {
    const auto& my_value1 = m_first1[i];
    const auto& my_value2 = m_first2[i];

    // FIXME_NVHPC using a ternary operator causes problems
    red_value_type rv = {i};
    if (m_predicate(my_value1, my_value2)) {
      rv = {::Kokkos::reduction_identity<IndexType>::min()};
    }

    m_reducer.join(red_value, rv);
  }

  KOKKOS_FUNCTION
  StdMismatchRedFunctor(IteratorType1 first1, IteratorType2 first2,
                        ReducerType reducer, BinaryPredicateType predicate)
      : m_first1(std::move(first1)),
        m_first2(std::move(first2)),
        m_reducer(std::move(reducer)),
        m_predicate(std::move(predicate)) {}
};

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class BinaryPredicateType>
::Kokkos::pair<IteratorType1, IteratorType2> mismatch_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
    IteratorType1 last1, IteratorType2 first2, IteratorType2 last2,
    BinaryPredicateType predicate) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first1, first2);
  Impl::static_assert_iterators_have_matching_difference_type(first1, first2);
  Impl::expect_valid_range(first1, last1);
  Impl::expect_valid_range(first2, last2);

  // aliases
  using return_type          = ::Kokkos::pair<IteratorType1, IteratorType2>;
  using index_type           = typename IteratorType1::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using functor_type =
      StdMismatchRedFunctor<index_type, IteratorType1, IteratorType2,
                            reducer_type, BinaryPredicateType>;

  // trivial case: note that this is important,
  // for OpenMPTarget, omitting special handling of
  // the trivial case was giving all sorts of strange stuff.
  const auto num_e1 = last1 - first1;
  const auto num_e2 = last2 - first2;
  if (num_e1 == 0 || num_e2 == 0) {
    return return_type(first1, first2);
  }

  // run
  const auto num_elemen_par_reduce = (num_e1 <= num_e2) ? num_e1 : num_e2;
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elemen_par_reduce),
      functor_type(first1, first2, reducer, std::move(predicate)), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  constexpr auto red_min = ::Kokkos::reduction_identity<index_type>::min();
  if (red_result.min_loc_true == red_min) {
    // in here means mismatch has not been found
    if (num_e1 == num_e2) {
      return return_type(last1, last2);
    } else if (num_e1 < num_e2) {
      return return_type(last1, first2 + num_e1);
    } else {
      return return_type(first1 + num_e2, last2);
    }
  } else {
    // in here means mismatch has been found
    return return_type(first1 + red_result.min_loc_true,
                       first2 + red_result.min_loc_true);
  }
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
::Kokkos::pair<IteratorType1, IteratorType2> mismatch_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType1 first1,
    IteratorType1 last1, IteratorType2 first2, IteratorType2 last2) {
  using value_type1 = typename IteratorType1::value_type;
  using value_type2 = typename IteratorType2::value_type;
  using pred_t      = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;
  return mismatch_impl(label, ex, first1, last1, first2, last2, pred_t());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
