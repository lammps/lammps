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

#ifndef KOKKOS_STD_ALGORITHMS_COUNT_IF_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_COUNT_IF_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class IteratorType, class Predicate>
struct StdCountIfFunctor {
  using index_type = typename IteratorType::difference_type;
  IteratorType m_first;
  Predicate m_predicate;

  KOKKOS_FUNCTION
  void operator()(index_type i, index_type& lsum) const {
    if (m_predicate(m_first[i])) {
      lsum++;
    }
  }

  KOKKOS_FUNCTION
  StdCountIfFunctor(IteratorType _first, Predicate _predicate)
      : m_first(std::move(_first)), m_predicate(std::move(_predicate)) {}
};

template <class ExecutionSpace, class IteratorType, class Predicate>
typename IteratorType::difference_type count_if_exespace_impl(
    const std::string& label, const ExecutionSpace& ex, IteratorType first,
    IteratorType last, Predicate predicate) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  typename IteratorType::difference_type count = 0;
  ::Kokkos::parallel_reduce(label,
                            RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                            // use CTAD
                            StdCountIfFunctor(first, predicate), count);
  ex.fence("Kokkos::count_if: fence after operation");

  return count;
}

template <class ExecutionSpace, class IteratorType, class T>
auto count_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last,
                         const T& value) {
  return count_if_exespace_impl(
      label, ex, first, last,
      ::Kokkos::Experimental::Impl::StdAlgoEqualsValUnaryPredicate<T>(value));
}

//
// team-level impl
//
template <class TeamHandleType, class IteratorType, class Predicate>
KOKKOS_FUNCTION typename IteratorType::difference_type count_if_team_impl(
    const TeamHandleType& teamHandle, IteratorType first, IteratorType last,
    Predicate predicate) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  typename IteratorType::difference_type count = 0;
  ::Kokkos::parallel_reduce(TeamThreadRange(teamHandle, 0, num_elements),
                            // use CTAD
                            StdCountIfFunctor(first, predicate), count);
  teamHandle.team_barrier();

  return count;
}

template <class TeamHandleType, class IteratorType, class T>
KOKKOS_FUNCTION auto count_team_impl(const TeamHandleType& teamHandle,
                                     IteratorType first, IteratorType last,
                                     const T& value) {
  return count_if_team_impl(
      teamHandle, first, last,
      ::Kokkos::Experimental::Impl::StdAlgoEqualsValUnaryPredicate<T>(value));
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
