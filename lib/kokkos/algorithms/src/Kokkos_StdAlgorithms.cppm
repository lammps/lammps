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

module;

#include <Kokkos_StdAlgorithms.hpp>

export module kokkos.std_algorithms;

export {
  namespace Kokkos::Experimental {
  using ::Kokkos::Experimental::adjacent_difference;
  using ::Kokkos::Experimental::adjacent_find;
  using ::Kokkos::Experimental::all_of;
  using ::Kokkos::Experimental::any_of;
  using ::Kokkos::Experimental::begin;
  using ::Kokkos::Experimental::cbegin;
  using ::Kokkos::Experimental::cend;
  using ::Kokkos::Experimental::copy;
  using ::Kokkos::Experimental::copy_backward;
  using ::Kokkos::Experimental::copy_if;
  using ::Kokkos::Experimental::copy_n;
  using ::Kokkos::Experimental::count;
  using ::Kokkos::Experimental::count_if;
  using ::Kokkos::Experimental::distance;
  using ::Kokkos::Experimental::end;
  using ::Kokkos::Experimental::equal;
  using ::Kokkos::Experimental::exclusive_scan;
  using ::Kokkos::Experimental::fill;
  using ::Kokkos::Experimental::fill_n;
  using ::Kokkos::Experimental::find;
  using ::Kokkos::Experimental::find_end;
  using ::Kokkos::Experimental::find_first_of;
  using ::Kokkos::Experimental::find_if;
  using ::Kokkos::Experimental::find_if_not;
  using ::Kokkos::Experimental::for_each;
  using ::Kokkos::Experimental::for_each_n;
  using ::Kokkos::Experimental::generate;
  using ::Kokkos::Experimental::generate_n;
  using ::Kokkos::Experimental::inclusive_scan;
  using ::Kokkos::Experimental::is_partitioned;
  using ::Kokkos::Experimental::is_sorted;
  using ::Kokkos::Experimental::is_sorted_until;
  using ::Kokkos::Experimental::iter_swap;
  using ::Kokkos::Experimental::lexicographical_compare;
  using ::Kokkos::Experimental::max_element;
  using ::Kokkos::Experimental::min_element;
  using ::Kokkos::Experimental::minmax_element;
  using ::Kokkos::Experimental::mismatch;
  using ::Kokkos::Experimental::move;
  using ::Kokkos::Experimental::move_backward;
  using ::Kokkos::Experimental::none_of;
  using ::Kokkos::Experimental::partition_copy;
  using ::Kokkos::Experimental::partition_point;
  using ::Kokkos::Experimental::reduce;
  using ::Kokkos::Experimental::remove;
  using ::Kokkos::Experimental::remove_copy;
  using ::Kokkos::Experimental::remove_copy_if;
  using ::Kokkos::Experimental::remove_if;
  using ::Kokkos::Experimental::replace;
  using ::Kokkos::Experimental::replace_copy;
  using ::Kokkos::Experimental::replace_copy_if;
  using ::Kokkos::Experimental::replace_if;
  using ::Kokkos::Experimental::reverse;
  using ::Kokkos::Experimental::reverse_copy;
  using ::Kokkos::Experimental::rotate;
  using ::Kokkos::Experimental::rotate_copy;
  using ::Kokkos::Experimental::search;
  using ::Kokkos::Experimental::search_n;
  using ::Kokkos::Experimental::shift_left;
  using ::Kokkos::Experimental::shift_right;
  using ::Kokkos::Experimental::swap_ranges;
  using ::Kokkos::Experimental::transform;
  using ::Kokkos::Experimental::transform_exclusive_scan;
  using ::Kokkos::Experimental::transform_inclusive_scan;
  using ::Kokkos::Experimental::transform_reduce;
  using ::Kokkos::Experimental::unique;
  using ::Kokkos::Experimental::unique_copy;
  }  // namespace Kokkos::Experimental
}
