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

#ifndef KOKKOS_STD_ALGORITHMS_IS_SORTED_HPP
#define KOKKOS_STD_ALGORITHMS_IS_SORTED_HPP

#include "impl/Kokkos_IsSorted.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename IteratorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool is_sorted(const ExecutionSpace& ex, IteratorType first,
               IteratorType last) {
  return Impl::is_sorted_exespace_impl("Kokkos::is_sorted_iterator_api_default",
                                       ex, first, last);
}

template <
    typename ExecutionSpace, typename IteratorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool is_sorted(const std::string& label, const ExecutionSpace& ex,
               IteratorType first, IteratorType last) {
  return Impl::is_sorted_exespace_impl(label, ex, first, last);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool is_sorted(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_exespace_impl("Kokkos::is_sorted_view_api_default", ex,
                                       KE::cbegin(view), KE::cend(view));
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool is_sorted(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_exespace_impl(label, ex, KE::cbegin(view),
                                       KE::cend(view));
}

template <
    typename ExecutionSpace, typename IteratorType, typename ComparatorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool is_sorted(const ExecutionSpace& ex, IteratorType first, IteratorType last,
               ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);
  return Impl::is_sorted_exespace_impl("Kokkos::is_sorted_iterator_api_default",
                                       ex, first, last, std::move(comp));
}

template <
    typename ExecutionSpace, typename IteratorType, typename ComparatorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool is_sorted(const std::string& label, const ExecutionSpace& ex,
               IteratorType first, IteratorType last, ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);
  return Impl::is_sorted_exespace_impl(label, ex, first, last, std::move(comp));
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename ComparatorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool is_sorted(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view,
               ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_not_openmptarget(ex);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_exespace_impl("Kokkos::is_sorted_view_api_default", ex,
                                       KE::cbegin(view), KE::cend(view),
                                       std::move(comp));
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename ComparatorType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool is_sorted(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view,
               ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_not_openmptarget(ex);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_exespace_impl(label, ex, KE::cbegin(view),
                                       KE::cend(view), std::move(comp));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool is_sorted(const TeamHandleType& teamHandle,
                               IteratorType first, IteratorType last) {
  return Impl::is_sorted_team_impl(teamHandle, first, last);
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool is_sorted(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_team_impl(teamHandle, KE::cbegin(view),
                                   KE::cend(view));
}

template <typename TeamHandleType, typename IteratorType,
          typename ComparatorType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool is_sorted(const TeamHandleType& teamHandle,
                               IteratorType first, IteratorType last,
                               ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(teamHandle);
  return Impl::is_sorted_team_impl(teamHandle, first, last, std::move(comp));
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          typename ComparatorType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool is_sorted(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType, Properties...>& view, ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_not_openmptarget(teamHandle);

  namespace KE = ::Kokkos::Experimental;
  return Impl::is_sorted_team_impl(teamHandle, KE::cbegin(view), KE::cend(view),
                                   std::move(comp));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
