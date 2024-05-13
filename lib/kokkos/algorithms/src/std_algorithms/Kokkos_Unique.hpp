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

#ifndef KOKKOS_STD_ALGORITHMS_UNIQUE_HPP
#define KOKKOS_STD_ALGORITHMS_UNIQUE_HPP

#include "impl/Kokkos_Unique.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set1: default predicate, accepting execution space
//
template <typename ExecutionSpace, typename IteratorType,
          std::enable_if_t<Impl::is_iterator_v<IteratorType> &&
                               is_execution_space<ExecutionSpace>::value,
                           int> = 0>
IteratorType unique(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last) {
  return Impl::unique_exespace_impl("Kokkos::unique_iterator_api_default", ex,
                                    first, last);
}

template <typename ExecutionSpace, typename IteratorType,
          std::enable_if_t<Impl::is_iterator_v<IteratorType> &&
                               is_execution_space<ExecutionSpace>::value,
                           int> = 0>
IteratorType unique(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last) {
  return Impl::unique_exespace_impl(label, ex, first, last);
}

template <typename ExecutionSpace, typename DataType, typename... Properties,
          std::enable_if_t<is_execution_space<ExecutionSpace>::value, int> = 0>
auto unique(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::unique_exespace_impl("Kokkos::unique_view_api_default", ex,
                                    begin(view), end(view));
}

template <typename ExecutionSpace, typename DataType, typename... Properties,
          std::enable_if_t<is_execution_space<ExecutionSpace>::value, int> = 0>
auto unique(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::unique_exespace_impl(label, ex, begin(view), end(view));
}

//
// overload set2: custom predicate, accepting execution space
//
template <typename ExecutionSpace, typename IteratorType,
          typename BinaryPredicate,
          std::enable_if_t<is_execution_space<ExecutionSpace>::value, int> = 0>
IteratorType unique(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last, BinaryPredicate pred) {
  return Impl::unique_exespace_impl("Kokkos::unique_iterator_api_default", ex,
                                    first, last, pred);
}

template <typename ExecutionSpace, typename IteratorType,
          typename BinaryPredicate,
          std::enable_if_t<is_execution_space<ExecutionSpace>::value, int> = 0>
IteratorType unique(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last,
                    BinaryPredicate pred) {
  return Impl::unique_exespace_impl(label, ex, first, last, pred);
}

template <typename ExecutionSpace, typename DataType, typename... Properties,
          typename BinaryPredicate,
          std::enable_if_t<is_execution_space<ExecutionSpace>::value, int> = 0>
auto unique(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            BinaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::unique_exespace_impl("Kokkos::unique_view_api_default", ex,
                                    begin(view), end(view), std::move(pred));
}

template <typename ExecutionSpace, typename DataType, typename... Properties,
          typename BinaryPredicate,
          std::enable_if_t<is_execution_space<ExecutionSpace>::value, int> = 0>
auto unique(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            BinaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::unique_exespace_impl(label, ex, begin(view), end(view),
                                    std::move(pred));
}

//
// overload set3: default predicate, accepting team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType,
          std::enable_if_t<Impl::is_iterator_v<IteratorType> &&
                               is_team_handle<TeamHandleType>::value,
                           int> = 0>
KOKKOS_FUNCTION IteratorType unique(const TeamHandleType& teamHandle,
                                    IteratorType first, IteratorType last) {
  return Impl::unique_team_impl(teamHandle, first, last);
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          std::enable_if_t<is_team_handle<TeamHandleType>::value, int> = 0>
KOKKOS_FUNCTION auto unique(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType, Properties...>& view) {
  return Impl::unique_team_impl(teamHandle, begin(view), end(view));
}

//
// overload set4: custom predicate, accepting team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType,
          typename BinaryPredicate,
          std::enable_if_t<is_team_handle<TeamHandleType>::value, int> = 0>
KOKKOS_FUNCTION IteratorType unique(const TeamHandleType& teamHandle,
                                    IteratorType first, IteratorType last,
                                    BinaryPredicate pred) {
  return Impl::unique_team_impl(teamHandle, first, last, std::move(pred));
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          typename BinaryPredicate,
          std::enable_if_t<is_team_handle<TeamHandleType>::value, int> = 0>
KOKKOS_FUNCTION auto unique(const TeamHandleType& teamHandle,
                            const ::Kokkos::View<DataType, Properties...>& view,
                            BinaryPredicate pred) {
  return Impl::unique_team_impl(teamHandle, begin(view), end(view),
                                std::move(pred));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
