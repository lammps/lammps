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

#ifndef KOKKOS_STD_ALGORITHMS_ROTATE_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_ROTATE_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Move.hpp>
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class ExecutionSpace, class IteratorType>
IteratorType rotate_with_pivot_in_left_half(const std::string& label,
                                            const ExecutionSpace& ex,
                                            IteratorType first,
                                            IteratorType n_first,
                                            IteratorType last) {
  /*
    This impl is specific for when the n_first iterator points to
    an element that is before or equal to the middle of the range.

    If we have:

    | 0 | 1 | 2 | 1 | 4 | 5 | 2 | 2 | 10 | -3 | 1 | -6 | -5 | 8 | 9 | 11 | *
      ^           ^              mid					   ^
    first       n_first							  last

    In step 1, we create a temporary view with extent = distance(n_first, last)
    and *move* the elements from [n_first, last) to tmp view, such that
    tmp view becomes:

    | 1 | 4 | 5 | 2 | 2 | 10 | -3 | 1 | -6 | -5 | 8 | 9 | 11 |

    In step 2, we move the elements in [first, n_first)
    to the new position where they are supposed to end up.

    In step 3, we move the elements from the tmp view to
    the range starting at first.
   */

  namespace KE                     = ::Kokkos::Experimental;
  const auto num_elements_on_left  = KE::distance(first, n_first);
  const auto num_elements_on_right = KE::distance(n_first, last);

  // create helper tmp view
  using value_type    = typename IteratorType::value_type;
  using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
  tmp_view_type tmp_view("rotate_impl_for_pivot_in_left_half_impl",
                         num_elements_on_right);
  using tmp_readwrite_iterator_type = decltype(begin(tmp_view));

  // index_type is the same and needed in all steps
  using index_type = typename IteratorType::difference_type;

  // stage 1
  using step1_func_type =
      StdMoveFunctor<index_type, IteratorType, tmp_readwrite_iterator_type>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_on_right),
      step1_func_type(n_first, begin(tmp_view)));

  // stage 2
  using step2_func_type =
      StdMoveFunctor<index_type, IteratorType, IteratorType>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_on_left),
      step2_func_type(first, first + num_elements_on_right));

  // step 3
  using step3_func_type =
      StdMoveFunctor<index_type, tmp_readwrite_iterator_type, IteratorType>;
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
                         step3_func_type(begin(tmp_view), first));

  ex.fence("Kokkos::rotate: fence after operation");
  return first + (last - n_first);
}

template <class ExecutionSpace, class IteratorType>
IteratorType rotate_with_pivot_in_right_half(const std::string& label,
                                             const ExecutionSpace& ex,
                                             IteratorType first,
                                             IteratorType n_first,
                                             IteratorType last) {
  /*
    This impl is specific for when the n_first iterator points to
    an element that is after the middle of the range.

    If we have:

    | 0 | 1 | 2 | 1 | 4 | 5 | 2 | 2 | 10 | -3 | 1 | -6 | -5 | 8 | 9 | 11 | *
      ^                          mid            ^                          ^
    first                                    n_first			  last

    In step 1, we create a temporary view with extent = distance(first, n_first)
    and *move* the elements from [first, n_first) to tmp view,
    such that tmp view becomes:

    | 0 | 1 | 2 | 1 | 4 | 5 | 2 | 2 | 10 | -3 | 1 |

    In step 2, we move the elements in [n_first, last)
    to the beginning where they are supposed to end up.

    In step 3, we move the elements from the tmp view to
    the range starting at first.
   */

  namespace KE                     = ::Kokkos::Experimental;
  const auto num_elements_on_left  = KE::distance(first, n_first);
  const auto num_elements_on_right = KE::distance(n_first, last);

  // create helper tmp view
  using value_type    = typename IteratorType::value_type;
  using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
  tmp_view_type tmp_view("rotate_impl_for_pivot_in_left_half_impl",
                         num_elements_on_left);
  using tmp_readwrite_iterator_type = decltype(begin(tmp_view));

  // index_type is the same and needed in all steps
  using index_type = typename IteratorType::difference_type;

  // stage 1
  using step1_func_type =
      StdMoveFunctor<index_type, IteratorType, tmp_readwrite_iterator_type>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_on_left),
      step1_func_type(first, begin(tmp_view)));

  // stage 2
  using step2_func_type =
      StdMoveFunctor<index_type, IteratorType, IteratorType>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_on_right),
      step2_func_type(n_first, first));

  // step 3:
  using step3_func_type =
      StdMoveFunctor<index_type, tmp_readwrite_iterator_type, IteratorType>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
      step3_func_type(begin(tmp_view), first + num_elements_on_right));

  ex.fence("Kokkos::rotate: fence after operation");
  return first + (last - n_first);
}

template <class ExecutionSpace, class IteratorType>
IteratorType rotate_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType n_first,
                         IteratorType last) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(first, n_first);
  Impl::expect_valid_range(n_first, last);

  namespace KE                     = ::Kokkos::Experimental;
  const auto num_elements          = KE::distance(first, last);
  const auto n_distance_from_first = KE::distance(first, n_first);
  if (n_distance_from_first <= num_elements / 2) {
    return rotate_with_pivot_in_left_half(label, ex, first, n_first, last);
  } else {
    return rotate_with_pivot_in_right_half(label, ex, first, n_first, last);
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
