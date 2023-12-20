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

#ifndef KOKKOS_STD_ALGORITHMS_NONE_OF_HPP
#define KOKKOS_STD_ALGORITHMS_NONE_OF_HPP

#include "impl/Kokkos_AllOfAnyOfNoneOf.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename IteratorType, typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool none_of(const ExecutionSpace& ex, IteratorType first, IteratorType last,
             Predicate predicate) {
  return Impl::none_of_exespace_impl("Kokkos::none_of_iterator_api_default", ex,
                                     first, last, predicate);
}

template <
    typename ExecutionSpace, typename IteratorType, typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool none_of(const std::string& label, const ExecutionSpace& ex,
             IteratorType first, IteratorType last, Predicate predicate) {
  return Impl::none_of_exespace_impl(label, ex, first, last, predicate);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool none_of(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::none_of_exespace_impl("Kokkos::none_of_view_api_default", ex,
                                     KE::cbegin(v), KE::cend(v),
                                     std::move(predicate));
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename Predicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool none_of(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::none_of_exespace_impl(label, ex, KE::cbegin(v), KE::cend(v),
                                     std::move(predicate));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType, typename Predicate>
KOKKOS_FUNCTION
    std::enable_if_t<::Kokkos::is_team_handle<TeamHandleType>::value, bool>
    none_of(const TeamHandleType& teamHandle, IteratorType first,
            IteratorType last, Predicate predicate) {
  return Impl::none_of_team_impl(teamHandle, first, last, predicate);
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          typename Predicate>
KOKKOS_FUNCTION
    std::enable_if_t<::Kokkos::is_team_handle<TeamHandleType>::value, bool>
    none_of(const TeamHandleType& teamHandle,
            const ::Kokkos::View<DataType, Properties...>& v,
            Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::none_of_team_impl(teamHandle, KE::cbegin(v), KE::cend(v),
                                 std::move(predicate));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
