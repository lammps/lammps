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

#ifndef KOKKOS_STD_ALGORITHMS_IS_PARTITIONED_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_IS_PARTITIONED_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class ReducerType, class PredicateType>
struct StdIsPartitionedFunctor {
  using red_value_type = typename ReducerType::value_type;
  using index_type     = typename IteratorType::difference_type;

  IteratorType m_first;
  ReducerType m_reducer;
  PredicateType m_p;

  KOKKOS_FUNCTION
  void operator()(const index_type i, red_value_type& redValue) const {
    const auto predicate_value = m_p(m_first[i]);
    constexpr index_type m_red_id_min =
        ::Kokkos::reduction_identity<index_type>::min();
    constexpr index_type m_red_id_max =
        ::Kokkos::reduction_identity<index_type>::max();

    // FIXME_NVHPC using a ternary operator causes problems
    red_value_type rv = {m_red_id_max, i};
    if (predicate_value) {
      rv = {i, m_red_id_min};
    }

    m_reducer.join(redValue, rv);
  }

  KOKKOS_FUNCTION
  StdIsPartitionedFunctor(IteratorType first, ReducerType reducer,
                          PredicateType p)
      : m_first(std::move(first)),
        m_reducer(std::move(reducer)),
        m_p(std::move(p)) {}
};

template <class ExecutionSpace, class IteratorType, class PredicateType>
bool is_partitioned_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last,
                         PredicateType pred) {
  // true if all elements in the range [first, last) that satisfy
  // the predicate "pred" appear before all elements that don't.
  // Also returns true if [first, last) is empty.
  // also true if all elements satisfy the predicate.

  // we implement it by finding:
  // - the max location where predicate is true  (max_loc_true)
  // - the min location where predicate is false (min_loc_false)
  // so the range is partitioned if max_loc_true < (min_loc_false)

  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // trivial case
  if (first == last) {
    return true;
  }

  // aliases
  using index_type           = typename IteratorType::difference_type;
  using reducer_type         = StdIsPartitioned<index_type>;
  using reduction_value_type = typename reducer_type::value_type;
  using func_t =
      StdIsPartitionedFunctor<IteratorType, reducer_type, PredicateType>;

  // run
  reduction_value_type red_result;
  reducer_type reducer(red_result);
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            func_t(first, reducer, pred), reducer);

  // fence not needed because reducing into scalar

  // decide and return
  constexpr index_type red_id_min =
      ::Kokkos::reduction_identity<index_type>::min();
  constexpr index_type red_id_max =
      ::Kokkos::reduction_identity<index_type>::max();

  if (red_result.max_loc_true != red_id_max &&
      red_result.min_loc_false != red_id_min) {
    return red_result.max_loc_true < red_result.min_loc_false;
  } else if (first + red_result.max_loc_true == --last) {
    return true;
  } else {
    return false;
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
