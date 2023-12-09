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

#ifndef KOKKOS_STD_ALGORITHMS_SHIFT_RIGHT_HPP
#define KOKKOS_STD_ALGORITHMS_SHIFT_RIGHT_HPP

#include "impl/Kokkos_ShiftRight.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename IteratorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType shift_right(const ExecutionSpace& ex, IteratorType first,
                         IteratorType last,
                         typename IteratorType::difference_type n) {
  return Impl::shift_right_exespace_impl(
      "Kokkos::shift_right_iterator_api_default", ex, first, last, n);
}

template <
    typename ExecutionSpace, typename IteratorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType shift_right(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last,
                         typename IteratorType::difference_type n) {
  return Impl::shift_right_exespace_impl(label, ex, first, last, n);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto shift_right(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 typename decltype(begin(view))::difference_type n) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::shift_right_exespace_impl("Kokkos::shift_right_view_api_default",
                                         ex, begin(view), end(view), n);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto shift_right(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 typename decltype(begin(view))::difference_type n) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::shift_right_exespace_impl(label, ex, begin(view), end(view), n);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION IteratorType
shift_right(const TeamHandleType& teamHandle, IteratorType first,
            IteratorType last, typename IteratorType::difference_type n) {
  return Impl::shift_right_team_impl(teamHandle, first, last, n);
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto shift_right(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType, Properties...>& view,
    typename decltype(begin(view))::difference_type n) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::shift_right_team_impl(teamHandle, begin(view), end(view), n);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
