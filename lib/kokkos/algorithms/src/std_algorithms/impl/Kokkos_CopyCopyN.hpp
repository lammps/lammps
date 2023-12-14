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

#ifndef KOKKOS_STD_ALGORITHMS_COPY_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_COPY_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class InputIterator, class OutputIterator>
struct StdCopyFunctor {
  // we can use difference type from InputIterator since
  // the calling functions below already static assert that
  // the iterators have matching difference type
  using index_type = typename InputIterator::difference_type;

  InputIterator m_first;
  OutputIterator m_dest_first;

  KOKKOS_FUNCTION
  void operator()(index_type i) const { m_dest_first[i] = m_first[i]; }

  KOKKOS_FUNCTION
  StdCopyFunctor(InputIterator _first, OutputIterator _dest_first)
      : m_first(std::move(_first)), m_dest_first(std::move(_dest_first)) {}
};

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator copy_exespace_impl(const std::string& label,
                                  const ExecutionSpace& ex, InputIterator first,
                                  InputIterator last, OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         // use CTAD
                         StdCopyFunctor(first, d_first));
  ex.fence("Kokkos::copy: fence after operation");

  // return
  return d_first + num_elements;
}

template <class ExecutionSpace, class InputIterator, class Size,
          class OutputIterator>
OutputIterator copy_n_exespace_impl(const std::string& label,
                                    const ExecutionSpace& ex,
                                    InputIterator first_from, Size count,
                                    OutputIterator first_dest) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);

  if (count > 0) {
    return copy_exespace_impl(label, ex, first_from, first_from + count,
                              first_dest);
  } else {
    return first_dest;
  }
}

//
// team-level impl
//
template <class TeamHandleType, class InputIterator, class OutputIterator>
KOKKOS_FUNCTION OutputIterator copy_team_impl(const TeamHandleType& teamHandle,
                                              InputIterator first,
                                              InputIterator last,
                                              OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(TeamThreadRange(teamHandle, 0, num_elements),
                         // use CTAD
                         StdCopyFunctor(first, d_first));
  teamHandle.team_barrier();

  // return
  return d_first + num_elements;
}

template <class TeamHandleType, class InputIterator, class Size,
          class OutputIterator>
KOKKOS_FUNCTION OutputIterator
copy_n_team_impl(const TeamHandleType& teamHandle, InputIterator first_from,
                 Size count, OutputIterator first_dest) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first_from,
                                                   first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);

  if (count > 0) {
    return copy_team_impl(teamHandle, first_from, first_from + count,
                          first_dest);
  } else {
    return first_dest;
  }
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
