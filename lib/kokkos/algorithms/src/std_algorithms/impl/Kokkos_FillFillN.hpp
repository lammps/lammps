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

#ifndef KOKKOS_STD_ALGORITHMS_FILL_AND_FILL_N_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_FILL_AND_FILL_N_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class InputIterator, class T>
struct StdFillFunctor {
  using index_type = typename InputIterator::difference_type;
  InputIterator m_first;
  T m_value;

  KOKKOS_FUNCTION
  void operator()(index_type i) const { m_first[i] = m_value; }

  KOKKOS_FUNCTION
  StdFillFunctor(InputIterator _first, T _value)
      : m_first(std::move(_first)), m_value(std::move(_value)) {}
};

//
// exespace impl
//
template <class ExecutionSpace, class IteratorType, class T>
void fill_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                        IteratorType first, IteratorType last, const T& value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  // run
  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(label,
                         RangePolicy<ExecutionSpace>(ex, 0, num_elements),
                         StdFillFunctor(first, value));
  ex.fence("Kokkos::fill: fence after operation");
}

template <class ExecutionSpace, class IteratorType, class SizeType, class T>
IteratorType fill_n_exespace_impl(const std::string& label,
                                  const ExecutionSpace& ex, IteratorType first,
                                  SizeType n, const T& value) {
  auto last = first + n;
  Impl::static_assert_random_access_and_accessible(ex, first);
  Impl::expect_valid_range(first, last);

  if (n <= 0) {
    return first;
  }

  fill_exespace_impl(label, ex, first, last, value);
  return last;
}

//
// team-level impl
//
template <class TeamHandleType, class IteratorType, class T>
KOKKOS_FUNCTION void fill_team_impl(const TeamHandleType& teamHandle,
                                    IteratorType first, IteratorType last,
                                    const T& value) {
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  const auto num_elements = Kokkos::Experimental::distance(first, last);
  ::Kokkos::parallel_for(TeamThreadRange(teamHandle, 0, num_elements),
                         StdFillFunctor(first, value));

  teamHandle.team_barrier();
}

template <class TeamHandleType, class IteratorType, class SizeType, class T>
KOKKOS_FUNCTION IteratorType fill_n_team_impl(const TeamHandleType& teamHandle,
                                              IteratorType first, SizeType n,
                                              const T& value) {
  auto last = first + n;
  Impl::static_assert_random_access_and_accessible(teamHandle, first);
  Impl::expect_valid_range(first, last);

  if (n <= 0) {
    return first;
  }

  fill_team_impl(teamHandle, first, last, value);
  return last;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
