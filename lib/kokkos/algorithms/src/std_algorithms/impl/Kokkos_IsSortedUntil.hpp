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

#ifndef KOKKOS_STD_ALGORITHMS_IS_SORTED_UNTIL_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_IS_SORTED_UNTIL_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <std_algorithms/Kokkos_Find.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class ComparatorType, class ReducerType>
struct StdIsSortedUntilFunctor {
  using index_type = typename IteratorType::difference_type;
  using value_type = typename ReducerType::value_type;

  IteratorType m_first;
  ComparatorType m_comparator;
  ReducerType m_reducer;

  KOKKOS_FUNCTION
  void operator()(const index_type i, value_type& reduction_result) const {
    const auto& val_i   = m_first[i];
    const auto& val_ip1 = m_first[i + 1];
    if (m_comparator(val_ip1, val_i)) {
      m_reducer.join(reduction_result, i);
    }
  }

  KOKKOS_FUNCTION
  StdIsSortedUntilFunctor(IteratorType first, ComparatorType comparator,
                          ReducerType reducer)
      : m_first(std::move(first)),
        m_comparator(std::move(comparator)),
        m_reducer(std::move(reducer)) {}
};

template <class ExecutionSpace, class IteratorType, class ComparatorType>
IteratorType is_sorted_until_impl(const std::string& label,
                                  const ExecutionSpace& ex, IteratorType first,
                                  IteratorType last, ComparatorType comp) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);

  // trivial case
  if (num_elements <= 1) {
    return last;
  }

  /*
    Do a par_reduce computing the *min* index that breaks the sorting.
    If such an index is found, then the range is sorted until that element.
    If no such index is found, then the range is sorted until the end.
  */
  using index_type = typename IteratorType::difference_type;
  index_type reduction_result;
  ::Kokkos::Min<index_type> reducer(reduction_result);
  ::Kokkos::parallel_reduce(
      label,
      // use num_elements-1 because each index handles i and i+1
      RangePolicy<ExecutionSpace>(ex, 0, num_elements - 1),
      // use CTAD
      StdIsSortedUntilFunctor(first, comp, reducer), reducer);

  /* If the reduction result is equal to the initial value,
     it means the range is sorted until the end */
  index_type reduction_result_init;
  reducer.init(reduction_result_init);
  if (reduction_result == reduction_result_init) {
    return last;
  } else {
    /* If such an index is found, then the range is sorted until there and
       we need to return an iterator past the element found so do +1 */
    return first + (reduction_result + 1);
  }
}

template <class ExecutionSpace, class IteratorType>
IteratorType is_sorted_until_impl(const std::string& label,
                                  const ExecutionSpace& ex, IteratorType first,
                                  IteratorType last) {
  using value_type = typename IteratorType::value_type;
  using pred_t     = Impl::StdAlgoLessThanBinaryPredicate<value_type>;
  return is_sorted_until_impl(label, ex, first, last, pred_t());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
