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

#ifndef KOKKOS_STD_ALGORITHMS_LEXICOGRAPHICAL_COMPARE_HPP
#define KOKKOS_STD_ALGORITHMS_LEXICOGRAPHICAL_COMPARE_HPP

#include "impl/Kokkos_LexicographicalCompare.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    class ExecutionSpace, class IteratorType1, class IteratorType2,
    std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool lexicographical_compare(const ExecutionSpace& ex, IteratorType1 first1,
                             IteratorType1 last1, IteratorType2 first2,
                             IteratorType2 last2) {
  return Impl::lexicographical_compare_exespace_impl(
      "Kokkos::lexicographical_compare_iterator_api_default", ex, first1, last1,
      first2, last2);
}

template <
    class ExecutionSpace, class IteratorType1, class IteratorType2,
    std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool lexicographical_compare(const std::string& label, const ExecutionSpace& ex,
                             IteratorType1 first1, IteratorType1 last1,
                             IteratorType2 first2, IteratorType2 last2) {
  return Impl::lexicographical_compare_exespace_impl(label, ex, first1, last1,
                                                     first2, last2);
}

template <
    class ExecutionSpace, class DataType1, class... Properties1,
    class DataType2, class... Properties2,
    std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool lexicographical_compare(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    const ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_exespace_impl(
      "Kokkos::lexicographical_compare_view_api_default", ex, KE::cbegin(view1),
      KE::cend(view1), KE::cbegin(view2), KE::cend(view2));
}

template <
    class ExecutionSpace, class DataType1, class... Properties1,
    class DataType2, class... Properties2,
    std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool lexicographical_compare(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    const ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_exespace_impl(
      label, ex, KE::cbegin(view1), KE::cend(view1), KE::cbegin(view2),
      KE::cend(view2));
}

template <
    class ExecutionSpace, class IteratorType1, class IteratorType2,
    class ComparatorType,
    std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool lexicographical_compare(const ExecutionSpace& ex, IteratorType1 first1,
                             IteratorType1 last1, IteratorType2 first2,
                             IteratorType2 last2, ComparatorType comp) {
  return Impl::lexicographical_compare_exespace_impl(
      "Kokkos::lexicographical_compare_iterator_api_default", ex, first1, last1,
      first2, last2, comp);
}

template <
    class ExecutionSpace, class IteratorType1, class IteratorType2,
    class ComparatorType,
    std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool lexicographical_compare(const std::string& label, const ExecutionSpace& ex,
                             IteratorType1 first1, IteratorType1 last1,
                             IteratorType2 first2, IteratorType2 last2,
                             ComparatorType comp) {
  return Impl::lexicographical_compare_exespace_impl(label, ex, first1, last1,
                                                     first2, last2, comp);
}

template <
    class ExecutionSpace, class DataType1, class... Properties1,
    class DataType2, class... Properties2, class ComparatorType,
    std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool lexicographical_compare(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    const ::Kokkos::View<DataType2, Properties2...>& view2,
    ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_exespace_impl(
      "Kokkos::lexicographical_compare_view_api_default", ex, KE::cbegin(view1),
      KE::cend(view1), KE::cbegin(view2), KE::cend(view2), comp);
}

template <
    class ExecutionSpace, class DataType1, class... Properties1,
    class DataType2, class... Properties2, class ComparatorType,
    std::enable_if_t<Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
bool lexicographical_compare(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    const ::Kokkos::View<DataType2, Properties2...>& view2,
    ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_exespace_impl(
      label, ex, KE::cbegin(view1), KE::cend(view1), KE::cbegin(view2),
      KE::cend(view2), comp);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <class TeamHandleType, class IteratorType1, class IteratorType2,
          std::enable_if_t<Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool lexicographical_compare(const TeamHandleType& teamHandle,
                                             IteratorType1 first1,
                                             IteratorType1 last1,
                                             IteratorType2 first2,
                                             IteratorType2 last2) {
  return Impl::lexicographical_compare_team_impl(teamHandle, first1, last1,
                                                 first2, last2);
}

template <class TeamHandleType, class DataType1, class... Properties1,
          class DataType2, class... Properties2,
          std::enable_if_t<Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool lexicographical_compare(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    const ::Kokkos::View<DataType2, Properties2...>& view2) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_team_impl(
      teamHandle, KE::cbegin(view1), KE::cend(view1), KE::cbegin(view2),
      KE::cend(view2));
}

template <class TeamHandleType, class IteratorType1, class IteratorType2,
          class ComparatorType,
          std::enable_if_t<Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool lexicographical_compare(
    const TeamHandleType& teamHandle, IteratorType1 first1, IteratorType1 last1,
    IteratorType2 first2, IteratorType2 last2, ComparatorType comp) {
  return Impl::lexicographical_compare_team_impl(teamHandle, first1, last1,
                                                 first2, last2, comp);
}

template <class TeamHandleType, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ComparatorType,
          std::enable_if_t<Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION bool lexicographical_compare(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view1,
    const ::Kokkos::View<DataType2, Properties2...>& view2,
    ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view2);

  namespace KE = ::Kokkos::Experimental;
  return Impl::lexicographical_compare_team_impl(
      teamHandle, KE::cbegin(view1), KE::cend(view1), KE::cbegin(view2),
      KE::cend(view2), comp);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
