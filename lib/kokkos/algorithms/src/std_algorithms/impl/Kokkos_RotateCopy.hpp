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

#ifndef KOKKOS_STD_ALGORITHMS_ROTATE_COPY_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_ROTATE_COPY_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class InputIterator, class OutputIterator>
struct StdRotateCopyFunctor {
  using index_type = typename InputIterator::difference_type;

  InputIterator m_first;
  InputIterator m_last;
  InputIterator m_first_n;
  OutputIterator m_dest_first;

  KOKKOS_FUNCTION
  void operator()(index_type i) const {
    const index_type shift = m_last - m_first_n;

    if (i < shift) {
      m_dest_first[i] = m_first_n[i];
    } else {
      m_dest_first[i] = m_first[i - shift];
    }
  }

  KOKKOS_FUNCTION
  StdRotateCopyFunctor(InputIterator first, InputIterator last,
                       InputIterator first_n, OutputIterator dest_first)
      : m_first(std::move(first)),
        m_last(std::move(last)),
        m_first_n(std::move(first_n)),
        m_dest_first(std::move(dest_first)) {}
};

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator rotate_copy_exespace_impl(
    const std::string& label, const ExecutionSpace& ex, InputIterator first,
    InputIterator n_first, InputIterator last, OutputIterator d_first) {
  /*
    algorithm is implemented as follows:

    first 	   n_first		last
    |		      |                  |
    o  o  o  o  o  o  o  o  o  o  o  o

    dest+0 -> first_n
    dest+1 -> first_n+1
    dest+2 -> first_n+2
    dest+3 -> first
    dest+4 -> first+1
    dest+5 -> first+2
    dest+6 -> first+3
    dest+7 -> first+4
    dest+8 -> first+5
    ...
    let shift = last - first_n;

    then we have:
    if (i < shift){
      *(dest_first + i) = *(first_n + i);
    }
    else{
      *(dest_first + i) = *(from + i - shift);
    }
  */

  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(first, n_first);
  Impl::expect_valid_range(n_first, last);

  if (first == last) {
    return d_first;
  }

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         StdRotateCopyFunctor(first, last, n_first, d_first));

  ex.fence("Kokkos::rotate_copy: fence after operation");

  // return
  return d_first + num_elements;
}

template <class TeamHandleType, class InputIterator, class OutputIterator>
KOKKOS_FUNCTION OutputIterator rotate_copy_team_impl(
    const TeamHandleType& teamHandle, InputIterator first,
    InputIterator n_first, InputIterator last, OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);
  Impl::expect_valid_range(first, n_first);
  Impl::expect_valid_range(n_first, last);

  if (first == last) {
    return d_first;
  }

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(TeamThreadRange(teamHandle, 0, num_elements),
                         StdRotateCopyFunctor(first, last, n_first, d_first));

  teamHandle.team_barrier();

  // return
  return d_first + num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
