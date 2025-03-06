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

#ifndef KOKKOS_STD_ALGORITHMS_ROTATE_HPP
#define KOKKOS_STD_ALGORITHMS_ROTATE_HPP

#include "impl/Kokkos_Rotate.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename IteratorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType rotate(const ExecutionSpace& ex, IteratorType first,
                    IteratorType n_first, IteratorType last) {
  return Impl::rotate_exespace_impl("Kokkos::rotate_iterator_api_default", ex,
                                    first, n_first, last);
}

template <
    typename ExecutionSpace, typename IteratorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType rotate(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType n_first,
                    IteratorType last) {
  return Impl::rotate_exespace_impl(label, ex, first, n_first, last);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto rotate(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            std::size_t n_location) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::rotate_exespace_impl("Kokkos::rotate_view_api_default", ex,
                                    begin(view), begin(view) + n_location,
                                    end(view));
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto rotate(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            std::size_t n_location) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::rotate_exespace_impl(label, ex, begin(view),
                                    begin(view) + n_location, end(view));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION IteratorType rotate(const TeamHandleType& teamHandle,
                                    IteratorType first, IteratorType n_first,
                                    IteratorType last) {
  return Impl::rotate_team_impl(teamHandle, first, n_first, last);
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto rotate(const TeamHandleType& teamHandle,
                            const ::Kokkos::View<DataType, Properties...>& view,
                            std::size_t n_location) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::rotate_team_impl(teamHandle, begin(view),
                                begin(view) + n_location, end(view));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
