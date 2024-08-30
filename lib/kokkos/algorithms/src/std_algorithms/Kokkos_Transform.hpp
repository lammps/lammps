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

#ifndef KOKKOS_STD_ALGORITHMS_TRANSFORM_HPP
#define KOKKOS_STD_ALGORITHMS_TRANSFORM_HPP

#include "impl/Kokkos_Transform.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename InputIterator, typename OutputIterator,
    typename UnaryOperation,
    std::enable_if_t<Impl::are_iterators_v<InputIterator, OutputIterator> &&
                         is_execution_space_v<ExecutionSpace>,
                     int> = 0>
OutputIterator transform(const ExecutionSpace& ex, InputIterator first1,
                         InputIterator last1, OutputIterator d_first,
                         UnaryOperation unary_op) {
  return Impl::transform_exespace_impl("Kokkos::transform_iterator_api_default",
                                       ex, first1, last1, d_first,
                                       std::move(unary_op));
}

template <
    typename ExecutionSpace, typename InputIterator, typename OutputIterator,
    typename UnaryOperation,
    std::enable_if_t<Impl::are_iterators_v<InputIterator, OutputIterator> &&
                         is_execution_space_v<ExecutionSpace>,
                     int> = 0>
OutputIterator transform(const std::string& label, const ExecutionSpace& ex,
                         InputIterator first1, InputIterator last1,
                         OutputIterator d_first, UnaryOperation unary_op) {
  return Impl::transform_exespace_impl(label, ex, first1, last1, d_first,
                                       std::move(unary_op));
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename UnaryOperation,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
auto transform(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source,
               const ::Kokkos::View<DataType2, Properties2...>& dest,
               UnaryOperation unary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_exespace_impl("Kokkos::transform_view_api_default", ex,
                                       begin(source), end(source), begin(dest),
                                       std::move(unary_op));
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename UnaryOperation,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
auto transform(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source,
               const ::Kokkos::View<DataType2, Properties2...>& dest,
               UnaryOperation unary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_exespace_impl(label, ex, begin(source), end(source),
                                       begin(dest), std::move(unary_op));
}

template <
    typename ExecutionSpace, typename InputIterator1, typename InputIterator2,
    typename OutputIterator, typename BinaryOperation,
    std::enable_if_t<
        Impl::are_iterators_v<InputIterator1, InputIterator2, OutputIterator> &&
            is_execution_space_v<ExecutionSpace>,
        int> = 0>
OutputIterator transform(const ExecutionSpace& ex, InputIterator1 first1,
                         InputIterator1 last1, InputIterator2 first2,
                         OutputIterator d_first, BinaryOperation binary_op) {
  return Impl::transform_exespace_impl("Kokkos::transform_iterator_api_default",
                                       ex, first1, last1, first2, d_first,
                                       std::move(binary_op));
}

template <
    typename ExecutionSpace, typename InputIterator1, typename InputIterator2,
    typename OutputIterator, typename BinaryOperation,
    std::enable_if_t<
        Impl::are_iterators_v<InputIterator1, InputIterator2, OutputIterator> &&
            is_execution_space_v<ExecutionSpace>,
        int> = 0>
OutputIterator transform(const std::string& label, const ExecutionSpace& ex,
                         InputIterator1 first1, InputIterator1 last1,
                         InputIterator2 first2, OutputIterator d_first,
                         BinaryOperation binary_op) {
  return Impl::transform_exespace_impl(label, ex, first1, last1, first2,
                                       d_first, std::move(binary_op));
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename DataType3,
          typename... Properties3, typename BinaryOperation,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
auto transform(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source1,
               const ::Kokkos::View<DataType2, Properties2...>& source2,
               const ::Kokkos::View<DataType3, Properties3...>& dest,
               BinaryOperation binary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source2);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_exespace_impl(
      "Kokkos::transform_view_api_default", ex, begin(source1), end(source1),
      begin(source2), begin(dest), std::move(binary_op));
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename DataType3,
          typename... Properties3, typename BinaryOperation,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
auto transform(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source1,
               const ::Kokkos::View<DataType2, Properties2...>& source2,
               const ::Kokkos::View<DataType3, Properties3...>& dest,
               BinaryOperation binary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source2);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_exespace_impl(label, ex, begin(source1), end(source1),
                                       begin(source2), begin(dest),
                                       std::move(binary_op));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <
    typename TeamHandleType, typename InputIterator, typename OutputIterator,
    typename UnaryOperation,
    std::enable_if_t<Impl::are_iterators_v<InputIterator, OutputIterator> &&
                         is_team_handle_v<TeamHandleType>,
                     int> = 0>
KOKKOS_FUNCTION OutputIterator transform(const TeamHandleType& teamHandle,
                                         InputIterator first1,
                                         InputIterator last1,
                                         OutputIterator d_first,
                                         UnaryOperation unary_op) {
  return Impl::transform_team_impl(teamHandle, first1, last1, d_first,
                                   std::move(unary_op));
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename UnaryOperation,
          std::enable_if_t<is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto transform(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& source,
    const ::Kokkos::View<DataType2, Properties2...>& dest,
    UnaryOperation unary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_team_impl(teamHandle, begin(source), end(source),
                                   begin(dest), std::move(unary_op));
}

template <
    typename TeamHandleType, typename InputIterator1, typename InputIterator2,
    typename OutputIterator, typename BinaryOperation,
    std::enable_if_t<
        Impl::are_iterators_v<InputIterator1, InputIterator2, OutputIterator> &&
            is_team_handle_v<TeamHandleType>,
        int> = 0>
KOKKOS_FUNCTION OutputIterator transform(const TeamHandleType& teamHandle,
                                         InputIterator1 first1,
                                         InputIterator1 last1,
                                         InputIterator2 first2,
                                         OutputIterator d_first,
                                         BinaryOperation binary_op) {
  return Impl::transform_team_impl(teamHandle, first1, last1, first2, d_first,
                                   std::move(binary_op));
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename DataType3,
          typename... Properties3, typename BinaryOperation,
          std::enable_if_t<is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto transform(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& source1,
    const ::Kokkos::View<DataType2, Properties2...>& source2,
    const ::Kokkos::View<DataType3, Properties3...>& dest,
    BinaryOperation binary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source2);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_team_impl(teamHandle, begin(source1), end(source1),
                                   begin(source2), begin(dest),
                                   std::move(binary_op));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
