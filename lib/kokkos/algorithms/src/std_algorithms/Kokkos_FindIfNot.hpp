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

#ifndef KOKKOS_STD_ALGORITHMS_FIND_IF_NOT_HPP
#define KOKKOS_STD_ALGORITHMS_FIND_IF_NOT_HPP

#include "impl/Kokkos_FindIfOrNot.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename IteratorType, typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType find_if_not(const ExecutionSpace& ex, IteratorType first,
                         IteratorType last, Predicate predicate) {
  return Impl::find_if_or_not_exespace_impl<false>(
      "Kokkos::find_if_not_iterator_api_default", ex, first, last,
      std::move(predicate));
}

template <
    typename ExecutionSpace, typename IteratorType, typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType find_if_not(const std::string& label, const ExecutionSpace& ex,
                         IteratorType first, IteratorType last,
                         Predicate predicate) {
  return Impl::find_if_or_not_exespace_impl<false>(label, ex, first, last,
                                                   std::move(predicate));
}

template <typename ExecutionSpace, typename DataType, typename... Properties,
          typename Predicate,
          std::enable_if_t<::Kokkos::is_execution_space<ExecutionSpace>::value,
                           int> = 0>
auto find_if_not(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v,
                 Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_if_or_not_exespace_impl<false>(
      "Kokkos::find_if_not_view_api_default", ex, KE::begin(v), KE::end(v),
      std::move(predicate));
}

template <typename ExecutionSpace, typename DataType, typename... Properties,
          typename Predicate,
          std::enable_if_t<::Kokkos::is_execution_space<ExecutionSpace>::value,
                           int> = 0>
auto find_if_not(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& v,
                 Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_if_or_not_exespace_impl<false>(
      label, ex, KE::begin(v), KE::end(v), std::move(predicate));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType, typename Predicate,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION IteratorType find_if_not(const TeamHandleType& teamHandle,
                                         IteratorType first, IteratorType last,
                                         Predicate predicate) {
  return Impl::find_if_or_not_team_impl<false>(teamHandle, first, last,
                                               std::move(predicate));
}

template <
    typename TeamHandleType, typename DataType, typename... Properties,
    typename Predicate,
    std::enable_if_t<::Kokkos::is_team_handle<TeamHandleType>::value, int> = 0>
KOKKOS_FUNCTION auto find_if_not(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType, Properties...>& v, Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_if_or_not_team_impl<false>(
      teamHandle, KE::begin(v), KE::end(v), std::move(predicate));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
