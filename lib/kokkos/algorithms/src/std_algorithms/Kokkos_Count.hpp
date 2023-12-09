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

#ifndef KOKKOS_STD_ALGORITHMS_COUNT_HPP
#define KOKKOS_STD_ALGORITHMS_COUNT_HPP

#include "impl/Kokkos_CountCountIf.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename IteratorType, typename T,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
typename IteratorType::difference_type count(const ExecutionSpace& ex,
                                             IteratorType first,
                                             IteratorType last,
                                             const T& value) {
  return Impl::count_exespace_impl("Kokkos::count_iterator_api_default", ex,
                                   first, last, value);
}

template <
    typename ExecutionSpace, typename IteratorType, typename T,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
typename IteratorType::difference_type count(const std::string& label,
                                             const ExecutionSpace& ex,
                                             IteratorType first,
                                             IteratorType last,
                                             const T& value) {
  return Impl::count_exespace_impl(label, ex, first, last, value);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename T,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto count(const ExecutionSpace& ex,
           const ::Kokkos::View<DataType, Properties...>& v, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_exespace_impl("Kokkos::count_view_api_default", ex,
                                   KE::cbegin(v), KE::cend(v), value);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename T,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto count(const std::string& label, const ExecutionSpace& ex,
           const ::Kokkos::View<DataType, Properties...>& v, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_exespace_impl(label, ex, KE::cbegin(v), KE::cend(v),
                                   value);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//

template <typename TeamHandleType, typename IteratorType, typename T,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION typename IteratorType::difference_type count(
    const TeamHandleType& teamHandle, IteratorType first, IteratorType last,
    const T& value) {
  return Impl::count_team_impl(teamHandle, first, last, value);
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          typename T,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto count(const TeamHandleType& teamHandle,
                           const ::Kokkos::View<DataType, Properties...>& v,
                           const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::count_team_impl(teamHandle, KE::cbegin(v), KE::cend(v), value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
