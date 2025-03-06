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

#ifndef KOKKOS_STD_ALGORITHMS_TRANSFORM_INCLUSIVE_SCAN_HPP
#define KOKKOS_STD_ALGORITHMS_TRANSFORM_INCLUSIVE_SCAN_HPP

#include "impl/Kokkos_TransformInclusiveScan.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//

// overload set 1 (no init value)
template <typename ExecutionSpace, typename InputIteratorType,
          typename OutputIteratorType, typename BinaryOpType,
          typename UnaryOpType,
          std::enable_if_t<
              Impl::are_iterators_v<InputIteratorType, OutputIteratorType>&& ::
                  Kokkos::is_execution_space_v<ExecutionSpace>,
              int> = 0>
OutputIteratorType transform_inclusive_scan(const ExecutionSpace& ex,
                                            InputIteratorType first,
                                            InputIteratorType last,
                                            OutputIteratorType first_dest,
                                            BinaryOpType binary_op,
                                            UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::transform_inclusive_scan_exespace_impl(
      "Kokkos::transform_inclusive_scan_custom_functors_iterator_api", ex,
      first, last, first_dest, binary_op, unary_op);
}

template <typename ExecutionSpace, typename InputIteratorType,
          typename OutputIteratorType, typename BinaryOpType,
          typename UnaryOpType,
          std::enable_if_t<
              Impl::are_iterators_v<InputIteratorType, OutputIteratorType>&& ::
                  Kokkos::is_execution_space_v<ExecutionSpace>,
              int> = 0>
OutputIteratorType transform_inclusive_scan(
    const std::string& label, const ExecutionSpace& ex, InputIteratorType first,
    InputIteratorType last, OutputIteratorType first_dest,
    BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::transform_inclusive_scan_exespace_impl(
      label, ex, first, last, first_dest, binary_op, unary_op);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename BinaryOpType,
    typename UnaryOpType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto transform_inclusive_scan(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_exespace_impl(
      "Kokkos::transform_inclusive_scan_custom_functors_view_api", ex,
      KE::cbegin(view_from), KE::cend(view_from), KE::begin(view_dest),
      binary_op, unary_op);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename BinaryOpType,
    typename UnaryOpType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto transform_inclusive_scan(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_exespace_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), binary_op, unary_op);
}

// overload set 2 (init value)
template <typename ExecutionSpace, typename InputIteratorType,
          typename OutputIteratorType, typename BinaryOpType,
          typename UnaryOpType, typename ValueType,

          std::enable_if_t<
              Impl::are_iterators_v<InputIteratorType, OutputIteratorType>&& ::
                  Kokkos::is_execution_space_v<ExecutionSpace>,
              int> = 0>
OutputIteratorType transform_inclusive_scan(
    const ExecutionSpace& ex, InputIteratorType first, InputIteratorType last,
    OutputIteratorType first_dest, BinaryOpType binary_op, UnaryOpType unary_op,
    ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(ex);
  static_assert(std::is_move_constructible_v<ValueType>,
                "ValueType must be move constructible.");

  return Impl::transform_inclusive_scan_exespace_impl(
      "Kokkos::transform_inclusive_scan_custom_functors_iterator_api", ex,
      first, last, first_dest, binary_op, unary_op, std::move(init_value));
}

template <typename ExecutionSpace, typename InputIteratorType,
          typename OutputIteratorType, typename BinaryOpType,
          typename UnaryOpType, typename ValueType,

          std::enable_if_t<
              Impl::are_iterators_v<InputIteratorType, OutputIteratorType>&& ::
                  Kokkos::is_execution_space_v<ExecutionSpace>,
              int> = 0>
OutputIteratorType transform_inclusive_scan(
    const std::string& label, const ExecutionSpace& ex, InputIteratorType first,
    InputIteratorType last, OutputIteratorType first_dest,
    BinaryOpType binary_op, UnaryOpType unary_op, ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(ex);
  static_assert(std::is_move_constructible_v<ValueType>,
                "ValueType must be move constructible.");

  return Impl::transform_inclusive_scan_exespace_impl(
      label, ex, first, last, first_dest, binary_op, unary_op,
      std::move(init_value));
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename BinaryOpType,
    typename UnaryOpType, typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto transform_inclusive_scan(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op, ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible_v<ValueType>,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_exespace_impl(
      "Kokkos::transform_inclusive_scan_custom_functors_view_api", ex,
      KE::cbegin(view_from), KE::cend(view_from), KE::begin(view_dest),
      binary_op, unary_op, std::move(init_value));
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename BinaryOpType,
    typename UnaryOpType, typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto transform_inclusive_scan(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op, ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible_v<ValueType>,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_exespace_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), binary_op, unary_op, std::move(init_value));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//

// overload set 1 (no init value)
template <typename TeamHandleType, typename InputIteratorType,
          typename OutputIteratorType, typename BinaryOpType,
          typename UnaryOpType,
          std::enable_if_t<
              Impl::are_iterators_v<InputIteratorType, OutputIteratorType> &&
                  Kokkos::is_team_handle_v<TeamHandleType>,
              int> = 0>
KOKKOS_FUNCTION OutputIteratorType transform_inclusive_scan(
    const TeamHandleType& teamHandle, InputIteratorType first,
    InputIteratorType last, OutputIteratorType first_dest,
    BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(teamHandle);

  return Impl::transform_inclusive_scan_team_impl(
      teamHandle, first, last, first_dest, binary_op, unary_op);
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename BinaryOpType,
          typename UnaryOpType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto transform_inclusive_scan(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(teamHandle);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_team_impl(
      teamHandle, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), binary_op, unary_op);
}

// overload set 2 (init value)
template <typename TeamHandleType, typename InputIteratorType,
          typename OutputIteratorType, typename BinaryOpType,
          typename UnaryOpType, typename ValueType,
          std::enable_if_t<
              Impl::are_iterators_v<InputIteratorType, OutputIteratorType> &&
                  Kokkos::is_team_handle_v<TeamHandleType>,
              int> = 0>
KOKKOS_FUNCTION OutputIteratorType transform_inclusive_scan(
    const TeamHandleType& teamHandle, InputIteratorType first,
    InputIteratorType last, OutputIteratorType first_dest,
    BinaryOpType binary_op, UnaryOpType unary_op, ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(teamHandle);
  static_assert(std::is_move_constructible_v<ValueType>,
                "ValueType must be move constructible.");

  return Impl::transform_inclusive_scan_team_impl(
      teamHandle, first, last, first_dest, binary_op, unary_op,
      std::move(init_value));
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename BinaryOpType,
          typename UnaryOpType, typename ValueType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto transform_inclusive_scan(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op, ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(teamHandle);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible_v<ValueType>,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_team_impl(
      teamHandle, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), binary_op, unary_op, std::move(init_value));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
