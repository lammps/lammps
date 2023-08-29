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

#ifndef KOKKOS_STD_ALGORITHMS_UNIQUE_COPY_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_UNIQUE_COPY_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include "Kokkos_CopyCopyN.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IndexType, class InputIt, class OutputIt,
          class BinaryPredicateType>
struct StdUniqueCopyFunctor {
  InputIt m_first_from;
  InputIt m_last_from;
  OutputIt m_first_dest;
  BinaryPredicateType m_pred;

  KOKKOS_FUNCTION
  StdUniqueCopyFunctor(InputIt first_from, InputIt last_from,
                       OutputIt first_dest, BinaryPredicateType pred)
      : m_first_from(std::move(first_from)),
        m_last_from(std::move(last_from)),
        m_first_dest(std::move(first_dest)),
        m_pred(std::move(pred)) {}

  KOKKOS_FUNCTION
  void operator()(const IndexType i, IndexType& update,
                  const bool final_pass) const {
    const auto& val_i   = m_first_from[i];
    const auto& val_ip1 = m_first_from[i + 1];

    if (final_pass) {
      if (!m_pred(val_i, val_ip1)) {
        m_first_dest[update] = val_i;
      }
    }

    if (!m_pred(val_i, val_ip1)) {
      update += 1;
    }
  }
};

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class PredicateType>
OutputIterator unique_copy_impl(const std::string& label,
                                const ExecutionSpace& ex, InputIterator first,
                                InputIterator last, OutputIterator d_first,
                                PredicateType pred) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // branch for trivial vs non trivial case
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  if (num_elements == 0) {
    return d_first;
  } else if (num_elements == 1) {
    return Impl::copy_impl("Kokkos::copy_from_unique_copy", ex, first, last,
                           d_first);
  } else {
    // aliases
    using index_type = typename InputIterator::difference_type;
    using func_type  = StdUniqueCopyFunctor<index_type, InputIterator,
                                           OutputIterator, PredicateType>;

    // note here that we run scan for num_elements - 1
    // because of the way we implement this, the last element is always needed.
    // We avoid performing checks inside functor that we are within limits
    // and run a "safe" scan and then copy the last element.
    const auto scan_size = num_elements - 1;
    index_type count     = 0;
    ::Kokkos::parallel_scan(label,
                            RangePolicy<ExecutionSpace>(ex, 0, scan_size),
                            func_type(first, last, d_first, pred), count);

    return Impl::copy_impl("Kokkos::copy_from_unique_copy", ex,
                           first + scan_size, last, d_first + count);
  }
}

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator unique_copy_impl(const std::string& label,
                                const ExecutionSpace& ex, InputIterator first,
                                InputIterator last, OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // aliases
  using value_type1 = typename InputIterator::value_type;
  using value_type2 = typename OutputIterator::value_type;

  // default binary predicate uses ==
  using binary_pred_t = StdAlgoEqualBinaryPredicate<value_type1, value_type2>;

  // run
  return unique_copy_impl(label, ex, first, last, d_first, binary_pred_t());
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
