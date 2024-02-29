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

#ifndef KOKKOS_STD_ALGORITHMS_REPLACE_COPY_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_REPLACE_COPY_IMPL_HPP

#include <Kokkos_Core.hpp>
#include "Kokkos_Constraints.hpp"
#include "Kokkos_HelperPredicates.hpp"
#include <std_algorithms/Kokkos_Distance.hpp>
#include <string>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <class InputIterator, class OutputIterator, class ValueType>
struct StdReplaceCopyFunctor {
  using index_type = typename InputIterator::difference_type;

  InputIterator m_first_from;
  OutputIterator m_first_dest;
  ValueType m_old_value;
  ValueType m_new_value;

  KOKKOS_FUNCTION
  void operator()(index_type i) const {
    const auto& myvalue_from = m_first_from[i];

    if (myvalue_from == m_old_value) {
      m_first_dest[i] = m_new_value;
    } else {
      m_first_dest[i] = myvalue_from;
    }
  }

  KOKKOS_FUNCTION
  StdReplaceCopyFunctor(InputIterator first_from, OutputIterator first_dest,
                        ValueType old_value, ValueType new_value)
      : m_first_from(std::move(first_from)),
        m_first_dest(std::move(first_dest)),
        m_old_value(std::move(old_value)),
        m_new_value(std::move(new_value)) {}
};

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType>
OutputIteratorType replace_copy_exespace_impl(const std::string& label,
                                              const ExecutionSpace& ex,
                                              InputIteratorType first_from,
                                              InputIteratorType last_from,
                                              OutputIteratorType first_dest,
                                              const ValueType& old_value,
                                              const ValueType& new_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(ex, first_from, first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_for(
      label, RangePolicy<ExecutionSpace>(ex, 0, num_elements),
      StdReplaceCopyFunctor(first_from, first_dest, old_value, new_value));
  ex.fence("Kokkos::replace_copy: fence after operation");

  // return
  return first_dest + num_elements;
}

template <class TeamHandleType, class InputIteratorType,
          class OutputIteratorType, class ValueType>
KOKKOS_FUNCTION OutputIteratorType replace_copy_team_impl(
    const TeamHandleType& teamHandle, InputIteratorType first_from,
    InputIteratorType last_from, OutputIteratorType first_dest,
    const ValueType& old_value, const ValueType& new_value) {
  // checks
  Impl::static_assert_random_access_and_accessible(teamHandle, first_from,
                                                   first_dest);
  Impl::static_assert_iterators_have_matching_difference_type(first_from,
                                                              first_dest);
  Impl::expect_valid_range(first_from, last_from);

  // run
  const auto num_elements =
      Kokkos::Experimental::distance(first_from, last_from);
  ::Kokkos::parallel_for(
      TeamThreadRange(teamHandle, 0, num_elements),
      StdReplaceCopyFunctor(first_from, first_dest, old_value, new_value));
  teamHandle.team_barrier();

  // return
  return first_dest + num_elements;
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
