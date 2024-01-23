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

#ifndef KOKKOS_STD_ALGORITHMS_EQUAL_HPP
#define KOKKOS_STD_ALGORITHMS_EQUAL_HPP

#include "impl/Kokkos_Equal.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2> &&
                               Kokkos::is_execution_space_v<ExecutionSpace>,
                           int> = 0>
bool equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
           IteratorType2 first2) {
  return Impl::equal_exespace_impl("Kokkos::equal_iterator_api_default", ex,
                                   first1, last1, first2);
}

template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_execution_space_v<ExecutionSpace>,
                           int> = 0>
bool equal(const std::string& label, const ExecutionSpace& ex,
           IteratorType1 first1, IteratorType1 last1, IteratorType2 first2) {
  return Impl::equal_exespace_impl(label, ex, first1, last1, first2);
}

template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2, typename BinaryPredicateType,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_execution_space_v<ExecutionSpace>,
                           int> = 0>
bool equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
           IteratorType2 first2, BinaryPredicateType predicate) {
  return Impl::equal_exespace_impl("Kokkos::equal_iterator_api_default", ex,
                                   first1, last1, first2, std::move(predicate));
}

template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2, typename BinaryPredicateType,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_execution_space_v<ExecutionSpace>,
                           int> = 0>
bool equal(const std::string& label, const ExecutionSpace& ex,
           IteratorType1 first1, IteratorType1 last1, IteratorType2 first2,
           BinaryPredicateType predicate) {
  return Impl::equal_exespace_impl(label, ex, first1, last1, first2,
                                   std::move(predicate));
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool equal(const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_exespace_impl("Kokkos::equal_view_api_default", ex,
                                   KE::cbegin(view1), KE::cend(view1),
                                   KE::cbegin(view2));
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool equal(const std::string& label, const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_exespace_impl(label, ex, KE::cbegin(view1),
                                   KE::cend(view1), KE::cbegin(view2));
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename BinaryPredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool equal(const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2,
           BinaryPredicateType predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_exespace_impl("Kokkos::equal_view_api_default", ex,
                                   KE::cbegin(view1), KE::cend(view1),
                                   KE::cbegin(view2), std::move(predicate));
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename BinaryPredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool equal(const std::string& label, const ExecutionSpace& ex,
           const ::Kokkos::View<DataType1, Properties1...>& view1,
           ::Kokkos::View<DataType2, Properties2...>& view2,
           BinaryPredicateType predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_exespace_impl(label, ex, KE::cbegin(view1),
                                   KE::cend(view1), KE::cbegin(view2),
                                   std::move(predicate));
}

template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_execution_space_v<ExecutionSpace>,
                           int> = 0>
bool equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
           IteratorType2 first2, IteratorType2 last2) {
  return Impl::equal_exespace_impl("Kokkos::equal_iterator_api_default", ex,
                                   first1, last1, first2, last2);
}

template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_execution_space_v<ExecutionSpace>,
                           int> = 0>
bool equal(const std::string& label, const ExecutionSpace& ex,
           IteratorType1 first1, IteratorType1 last1, IteratorType2 first2,
           IteratorType2 last2) {
  return Impl::equal_exespace_impl(label, ex, first1, last1, first2, last2);
}

template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2, typename BinaryPredicateType,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_execution_space_v<ExecutionSpace>,
                           int> = 0>
bool equal(const ExecutionSpace& ex, IteratorType1 first1, IteratorType1 last1,
           IteratorType2 first2, IteratorType2 last2,
           BinaryPredicateType predicate) {
  return Impl::equal_exespace_impl("Kokkos::equal_iterator_api_default", ex,
                                   first1, last1, first2, last2,
                                   std::move(predicate));
}

template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2, typename BinaryPredicateType,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_execution_space_v<ExecutionSpace>,
                           int> = 0>
bool equal(const std::string& label, const ExecutionSpace& ex,
           IteratorType1 first1, IteratorType1 last1, IteratorType2 first2,
           IteratorType2 last2, BinaryPredicateType predicate) {
  return Impl::equal_exespace_impl(label, ex, first1, last1, first2, last2,
                                   std::move(predicate));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_team_handle_v<TeamHandleType>,
                           int> = 0>
KOKKOS_FUNCTION bool equal(const TeamHandleType& teamHandle,
                           IteratorType1 first1, IteratorType1 last1,
                           IteratorType2 first2) {
  return Impl::equal_team_impl(teamHandle, first1, last1, first2);
}

template <typename TeamHandleType, typename IteratorType1,
          typename IteratorType2, typename BinaryPredicateType,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_team_handle_v<TeamHandleType>,
                           int> = 0>
KOKKOS_FUNCTION bool equal(const TeamHandleType& teamHandle,
                           IteratorType1 first1, IteratorType1 last1,
                           IteratorType2 first2,
                           BinaryPredicateType predicate) {
  return Impl::equal_team_impl(teamHandle, first1, last1, first2,
                               std::move(predicate));
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool equal(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_team_impl(teamHandle, KE::cbegin(view1), KE::cend(view1),
                               KE::cbegin(view2));
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2,
          typename BinaryPredicateType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool equal(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    ::Kokkos::View<DataType2, Properties2...>& view2,
    BinaryPredicateType predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::equal_team_impl(teamHandle, KE::cbegin(view1), KE::cend(view1),
                               KE::cbegin(view2), std::move(predicate));
}

template <typename TeamHandleType, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_team_handle_v<TeamHandleType>,
                           int> = 0>
KOKKOS_FUNCTION bool equal(const TeamHandleType& teamHandle,
                           IteratorType1 first1, IteratorType1 last1,
                           IteratorType2 first2, IteratorType2 last2) {
  return Impl::equal_team_impl(teamHandle, first1, last1, first2, last2);
}

template <typename TeamHandleType, typename IteratorType1,
          typename IteratorType2, typename BinaryPredicateType,
          std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators_v<
                               IteratorType1, IteratorType2>&& ::Kokkos::
                               is_team_handle_v<TeamHandleType>,
                           int> = 0>
KOKKOS_FUNCTION bool equal(const TeamHandleType& teamHandle,
                           IteratorType1 first1, IteratorType1 last1,
                           IteratorType2 first2, IteratorType2 last2,
                           BinaryPredicateType predicate) {
  return Impl::equal_team_impl(teamHandle, first1, last1, first2, last2,
                               std::move(predicate));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
