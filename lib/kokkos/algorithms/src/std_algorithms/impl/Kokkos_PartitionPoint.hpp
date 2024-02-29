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

#ifndef KOKKOS_STD_ALGORITHMS_PARTITION_POINT_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_PARTITION_POINT_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class ReducerType, class PredicateType>
struct StdPartitionPointFunctor {
  using red_value_type = typename ReducerType::value_type;
  using index_type     = typename IteratorType::difference_type;

  IteratorType m_first;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& redValue) const {
    const auto predicate_value = m_p(m_first[i]);

    // FIXME_NVHPC using a ternary operator causes problems
    red_value_type rv = {i};
    if (predicate_value) {
      rv = {::Kokkos::reduction_identity<index_type>::min()};
    }

    m_reducer.join(redValue, rv);
  }

  KOKKOS_FUNCTION
  StdPartitionPointFunctor(IteratorType first, ReducerType reducer,
                           PredicateType p)
      : m_first(std::move(first)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType partition_point_exespace_impl(const std::string& label,
                                           const ExecutionSpace& ex,
                                           IteratorType first,
                                           IteratorType last,
                                           PredicateType pred) {
  // locates the end of the first partition, that is, the first
  // element that does not satisfy p or last if all elements satisfy p.
  // Implementation below finds the first location where p is false.

  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return first;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = StdPartitionPoint<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t =
      StdPartitionPointFunctor<IteratorType, reducer_type, PredicateType>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer, pred), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  if (red_result.min_loc_false ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    // if all elements are true, return last
    return last;
  } else {
    return first + red_result.min_loc_false;
  }
}

template <class TeamHandleType, class IteratorType, class PredicateType>
KOKKOS_FUNCTION IteratorType
partition_point_team_impl(const TeamHandleType& teamHandle, IteratorType first,
                          IteratorType last, PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  if (first == last) {
    return first;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = StdPartitionPoint<index_type>;
  using reduction_value_type = typename reducer_type::value_type;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(TeamThreadRange(teamHandle, 0, num_elements),
                            StdPartitionPointFunctor(first, reducer, pred),
                            reducer);

  // fence not needed because reducing into scalar

  // decide and return
  if (red_result.min_loc_false ==
      ::Kokkos::reduction_identity<index_type>::min()) {
    // if all elements are true, return last
    return last;
  } else {
    return first + red_result.min_loc_false;
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
