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

#ifndef KOKKOS_STD_ALGORITHMS_FILL_HPP
#define KOKKOS_STD_ALGORITHMS_FILL_HPP

#include "impl/Kokkos_FillFillN.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename IteratorType, typename T,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void fill(const ExecutionSpace& ex, IteratorType first, IteratorType last,
          const T& value) {
  Impl::fill_exespace_impl("Kokkos::fill_iterator_api_default", ex, first, last,
                           value);
}

template <
    typename ExecutionSpace, typename IteratorType, typename T,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void fill(const std::string& label, const ExecutionSpace& ex,
          IteratorType first, IteratorType last, const T& value) {
  Impl::fill_exespace_impl(label, ex, first, last, value);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename T,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void fill(const ExecutionSpace& ex,
          const ::Kokkos::View<DataType, Properties...>& view, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::fill_exespace_impl("Kokkos::fill_view_api_default", ex, begin(view),
                           end(view), value);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename T,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void fill(const std::string& label, const ExecutionSpace& ex,
          const ::Kokkos::View<DataType, Properties...>& view, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::fill_exespace_impl(label, ex, begin(view), end(view), value);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType, typename T,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION void fill(const TeamHandleType& th, IteratorType first,
                          IteratorType last, const T& value) {
  Impl::fill_team_impl(th, first, last, value);
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          typename T,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION void fill(const TeamHandleType& th,
                          const ::Kokkos::View<DataType, Properties...>& view,
                          const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::fill_team_impl(th, begin(view), end(view), value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
