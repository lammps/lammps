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

#ifndef KOKKOS_STD_ALGORITHMS_FIND_FIRST_OF_HPP
#define KOKKOS_STD_ALGORITHMS_FIND_FIRST_OF_HPP

#include "impl/Kokkos_FindFirstOf.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//

// overload set 1: no binary predicate passed
template <
    typename ExecutionSpace, typename IteratorType1, typename IteratorType2,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType1 find_first_of(const ExecutionSpace& ex, IteratorType1 first,
                            IteratorType1 last, IteratorType2 s_first,
                            IteratorType2 s_last) {
  return Impl::find_first_of_exespace_impl(
      "Kokkos::find_first_of_iterator_api_default", ex, first, last, s_first,
      s_last);
}

template <
    typename ExecutionSpace, typename IteratorType1, typename IteratorType2,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType1 find_first_of(const std::string& label, const ExecutionSpace& ex,
                            IteratorType1 first, IteratorType1 last,
                            IteratorType2 s_first, IteratorType2 s_last) {
  return Impl::find_first_of_exespace_impl(label, ex, first, last, s_first,
                                           s_last);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto find_first_of(const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& view,
                   const ::Kokkos::View<DataType2, Properties2...>& s_view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_exespace_impl(
      "Kokkos::find_first_of_view_api_default", ex, KE::begin(view),
      KE::end(view), KE::begin(s_view), KE::end(s_view));
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto find_first_of(const std::string& label, const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& view,
                   const ::Kokkos::View<DataType2, Properties2...>& s_view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_exespace_impl(label, ex, KE::begin(view),
                                           KE::end(view), KE::begin(s_view),
                                           KE::end(s_view));
}

// overload set 2: binary predicate passed
template <
    typename ExecutionSpace, typename IteratorType1, typename IteratorType2,
    typename BinaryPredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType1 find_first_of(const ExecutionSpace& ex, IteratorType1 first,
                            IteratorType1 last, IteratorType2 s_first,
                            IteratorType2 s_last,
                            const BinaryPredicateType& pred) {
  return Impl::find_first_of_exespace_impl(
      "Kokkos::find_first_of_iterator_api_default", ex, first, last, s_first,
      s_last, pred);
}

template <
    typename ExecutionSpace, typename IteratorType1, typename IteratorType2,
    typename BinaryPredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType1 find_first_of(const std::string& label, const ExecutionSpace& ex,
                            IteratorType1 first, IteratorType1 last,
                            IteratorType2 s_first, IteratorType2 s_last,
                            const BinaryPredicateType& pred) {
  return Impl::find_first_of_exespace_impl(label, ex, first, last, s_first,
                                           s_last, pred);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename BinaryPredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto find_first_of(const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& view,
                   const ::Kokkos::View<DataType2, Properties2...>& s_view,
                   const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_exespace_impl(
      "Kokkos::find_first_of_view_api_default", ex, KE::begin(view),
      KE::end(view), KE::begin(s_view), KE::end(s_view), pred);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename BinaryPredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto find_first_of(const std::string& label, const ExecutionSpace& ex,
                   const ::Kokkos::View<DataType1, Properties1...>& view,
                   const ::Kokkos::View<DataType2, Properties2...>& s_view,
                   const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_exespace_impl(label, ex, KE::begin(view),
                                           KE::end(view), KE::begin(s_view),
                                           KE::end(s_view), pred);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//

// overload set 1: no binary predicate passed
template <typename TeamHandleType, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION IteratorType1 find_first_of(const TeamHandleType& teamHandle,
                                            IteratorType1 first,
                                            IteratorType1 last,
                                            IteratorType2 s_first,
                                            IteratorType2 s_last) {
  return Impl::find_first_of_team_impl(teamHandle, first, last, s_first,
                                       s_last);
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto find_first_of(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view,
    const ::Kokkos::View<DataType2, Properties2...>& s_view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_team_impl(teamHandle, KE::begin(view),
                                       KE::end(view), KE::begin(s_view),
                                       KE::end(s_view));
}

// overload set 2: binary predicate passed
template <typename TeamHandleType, typename IteratorType1,
          typename IteratorType2, typename BinaryPredicateType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>

KOKKOS_FUNCTION IteratorType1 find_first_of(const TeamHandleType& teamHandle,
                                            IteratorType1 first,
                                            IteratorType1 last,
                                            IteratorType2 s_first,
                                            IteratorType2 s_last,
                                            const BinaryPredicateType& pred) {
  return Impl::find_first_of_team_impl(teamHandle, first, last, s_first, s_last,
                                       pred);
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2,
          typename BinaryPredicateType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto find_first_of(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view,
    const ::Kokkos::View<DataType2, Properties2...>& s_view,
    const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(s_view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_first_of_team_impl(teamHandle, KE::begin(view),
                                       KE::end(view), KE::begin(s_view),
                                       KE::end(s_view), pred);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
