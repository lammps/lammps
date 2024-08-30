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

#ifndef KOKKOS_STD_ALGORITHMS_PARTITION_COPY_HPP
#define KOKKOS_STD_ALGORITHMS_PARTITION_COPY_HPP

#include "impl/Kokkos_PartitionCopy.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename InputIteratorType,
    typename OutputIteratorTrueType, typename OutputIteratorFalseType,
    typename PredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType> partition_copy(
    const ExecutionSpace& ex, InputIteratorType from_first,
    InputIteratorType from_last, OutputIteratorTrueType to_first_true,
    OutputIteratorFalseType to_first_false, PredicateType p) {
  return Impl::partition_copy_exespace_impl(
      "Kokkos::partition_copy_iterator_api_default", ex, from_first, from_last,
      to_first_true, to_first_false, std::move(p));
}

template <
    typename ExecutionSpace, typename InputIteratorType,
    typename OutputIteratorTrueType, typename OutputIteratorFalseType,
    typename PredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType> partition_copy(
    const std::string& label, const ExecutionSpace& ex,
    InputIteratorType from_first, InputIteratorType from_last,
    OutputIteratorTrueType to_first_true,
    OutputIteratorFalseType to_first_false, PredicateType p) {
  return Impl::partition_copy_exespace_impl(label, ex, from_first, from_last,
                                            to_first_true, to_first_false,
                                            std::move(p));
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename DataType3,
    typename... Properties3, typename PredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto partition_copy(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest_true,
    const ::Kokkos::View<DataType3, Properties3...>& view_dest_false,
    PredicateType p) {
  return Impl::partition_copy_exespace_impl(
      "Kokkos::partition_copy_view_api_default", ex, cbegin(view_from),
      cend(view_from), begin(view_dest_true), begin(view_dest_false),
      std::move(p));
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename DataType3,
    typename... Properties3, typename PredicateType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto partition_copy(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest_true,
    const ::Kokkos::View<DataType3, Properties3...>& view_dest_false,
    PredicateType p) {
  return Impl::partition_copy_exespace_impl(
      label, ex, cbegin(view_from), cend(view_from), begin(view_dest_true),
      begin(view_dest_false), std::move(p));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename InputIteratorType,
          typename OutputIteratorTrueType, typename OutputIteratorFalseType,
          typename PredicateType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION ::Kokkos::pair<OutputIteratorTrueType, OutputIteratorFalseType>
partition_copy(const TeamHandleType& teamHandle, InputIteratorType from_first,
               InputIteratorType from_last,
               OutputIteratorTrueType to_first_true,
               OutputIteratorFalseType to_first_false, PredicateType p) {
  return Impl::partition_copy_team_impl(teamHandle, from_first, from_last,
                                        to_first_true, to_first_false,
                                        std::move(p));
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename DataType3,
          typename... Properties3, typename PredicateType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto partition_copy(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest_true,
    const ::Kokkos::View<DataType3, Properties3...>& view_dest_false,
    PredicateType p) {
  return Impl::partition_copy_team_impl(teamHandle, cbegin(view_from),
                                        cend(view_from), begin(view_dest_true),
                                        begin(view_dest_false), std::move(p));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
