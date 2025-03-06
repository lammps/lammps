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

#ifndef KOKKOS_STD_ALGORITHMS_ADJACENT_FIND_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_ADJACENT_FIND_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class ReducerType, class PredicateType>
struct StdAdjacentFindFunctor {
  using index_type     = typename IteratorType::difference_type;
  using red_value_type = typename ReducerType::value_type;

  IteratorType m_first;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& red_value) const {
    const auto& my_value   = m_first[i];
    const auto& next_value = m_first[i + 1];
    const bool are_equal   = m_p(my_value, next_value);

    // FIXME_NVHPC using a ternary operator causes problems
    red_value_type value = {::Kokkos::reduction_identity<index_type>::min()};
    if (are_equal) {
      value.min_loc_true = i;
    }

    m_reducer.join(red_value, value);
  }

  KOKKOS_FUNCTION
  StdAdjacentFindFunctor(IteratorType first, ReducerType reducer,
                         PredicateType p)
      : m_first(std::move(first)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

//
// exespace impl
//
template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType adjacent_find_exespace_impl(const std::string& label,
                                         const ExecutionSpace& ex,
                                         IteratorType first, IteratorType last,
                                         PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);

  if (num_elements <= 1) {
    return last;
  }

  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;

  reduction_value_type red_result;
  reducer_type reducer(red_result);

  // note that we use below num_elements-1 because
  // each index i in the reduction checks i and (i+1).
  ::Kokkos::parallel_reduce(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements - 1),
      // use CTAD
      StdAdjacentFindFunctor(first, reducer, pred), reducer);

  // fence not needed because reducing into scalar
  if (red_result.min_loc_true ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    return last;
  } else {
    return first + red_result.min_loc_true;
  }
}

template <class ExecutionSpace, class IteratorType>
IteratorType adjacent_find_exespace_impl(const std::string& label,
                                         const ExecutionSpace& ex,
                                         IteratorType first,
                                         IteratorType last) {
  using value_type     = typename IteratorType::value_type;
  using default_pred_t = StdAlgoEqualBinaryPredicate<value_type>;
  return adjacent_find_exespace_impl(label, ex, first, last, default_pred_t());
}

//
// team impl
//
template <class TeamHandleType, class IteratorType, class PredicateType>
KOKKOS_FUNCTION IteratorType
adjacent_find_team_impl(const TeamHandleType& teamHandle, IteratorType first,
                        IteratorType last, PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);

  if (num_elements <= 1) {
    return last;
  }

  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = FirstLoc<index_type>;
  using reduction_value_type = typename reducer_type::value_type;

  reduction_value_type red_result;
  reducer_type reducer(red_result);

  // note that we use below num_elements-1 because
  // each index i in the reduction checks i and (i+1).
  ::Kokkos::parallel_reduce(TeamThreadRange(teamHandle, 0, num_elements - 1),
                            // use CTAD
                            StdAdjacentFindFunctor(first, reducer, pred),
                            reducer);

  teamHandle.team_barrier();

  if (red_result.min_loc_true ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    return last;
  } else {
    return first + red_result.min_loc_true;
  }
}

template <class TeamHandleType, class IteratorType>
KOKKOS_FUNCTION IteratorType adjacent_find_team_impl(
    const TeamHandleType& teamHandle, IteratorType first, IteratorType last) {
  using value_type     = typename IteratorType::value_type;
  using default_pred_t = StdAlgoEqualBinaryPredicate<value_type>;
  return adjacent_find_team_impl(teamHandle, first, last, default_pred_t());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
