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

#ifndef KOKKOS_STD_ALGORITHMS_REVERSE_COPY_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_REVERSE_COPY_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class InputIterator, class OutputIterator>
struct StdReverseCopyFunctor {
  using index_type = typename InputIterator::difference_type;
  static_assert(std::is_signed<index_type>::value,
                "Kokkos: StdReverseCopyFunctor requires signed index type");

  InputIterator m_last;
  OutputIterator m_dest_first;

  KOKKOS_FUNCTION
  void operator()(index_type i) const { m_dest_first[i] = m_last[-1 - i]; }

  KOKKOS_FUNCTION
  StdReverseCopyFunctor(InputIterator _last, OutputIterator _dest_first)
      : m_last(std::move(_last)), m_dest_first(std::move(_dest_first)) {}
};

template <class ExecutionSpace, class InputIterator, class OutputIterator>
OutputIterator reverse_copy_exespace_impl(const std::string& label,
                                          const ExecutionSpace& ex,
                                          InputIterator first,
                                          InputIterator last,
                                          OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         StdReverseCopyFunctor(last, d_first));
  ex.fence("Kokkos::reverse_copy: fence after operation");

  // return
  return d_first + num_elements;
}

template <class TeamHandleType, class InputIterator, class OutputIterator>
KOKKOS_FUNCTION OutputIterator
reverse_copy_team_impl(const TeamHandleType& teamHandle, InputIterator first,
                       InputIterator last, OutputIterator d_first) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first, d_first);
  Impl::static_assert_iterators_have_matching_difference_type(first, d_first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(TeamThreadRange(teamHandle, 0, num_elements),
                         StdReverseCopyFunctor(last, d_first));
  teamHandle.team_barrier();

  // return
  return d_first + num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
