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

#ifndef KOKKOS_STD_ALGORITHMS_SHIFT_LEFT_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_SHIFT_LEFT_IMPL_HPP

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
IteratorType shift_left_impl(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, IteratorType last,
                             typename IteratorType::difference_type n) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);
  KOKKOS_EXPECTS(n >= 0);

  // handle trivial cases
  if (n == 0) {
    return last;
  }

  if (n >= Kokkos::Experimental::distance(first, last)) {
    return first;
  }

  /*
    Suppose that n = 5, and our [first,last) spans:

    | 0  | 1  |  2 | 1  | 2  | 1  | 2  | 2  | 10 | -3 | 1  | -6 | *
      ^                         				  ^
    first							 last

    shift_left modifies the range such that we have this data:
    | 1  | 2  | 2  | 10  | -3 | 1  | -6 | x | x  | x  | x  |  x | *
                                          ^
                                   return it pointing here


    and returns an iterator pointing to one past the new end.
    Note: elements marked x are in undefined state because have been moved.

    We implement this in two steps:
    step 1:
      we create a temporary view with extent = distance(first+n, last)
      and *move* assign the elements from [first+n, last) to tmp view, such that
      tmp view becomes:

      | 1  | 2  | 2  | 10  | -3 | 1  | -6 |

    step 2:
      move elements of tmp view back to range starting at first.
   */

  const auto num_elements_to_move =
      ::Kokkos::Experimental::distance(first + n, last);

  // create tmp view
  using value_type    = typename IteratorType::value_type;
  using tmp_view_type = Kokkos::View<value_type*, ExecutionSpace>;
  tmp_view_type tmp_view("shift_left_impl", num_elements_to_move);
  using tmp_readwrite_iterator_type = decltype(begin(tmp_view));

  using index_type = typename IteratorType::difference_type;

  // step 1
  using step1_func_type =
      StdMoveFunctor<index_type, IteratorType, tmp_readwrite_iterator_type>;
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements_to_move),
      step1_func_type(first + n, begin(tmp_view)));

  // step 2
  using step2_func_type =
      StdMoveFunctor<index_type, tmp_readwrite_iterator_type, IteratorType>;
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, tmp_view.extent(0)),
                         step2_func_type(begin(tmp_view), first));

  ex.fence("Kokkos::shift_left: fence after operation");

  return last - n;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
