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

// overload set 1 (no init value)
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class BinaryOpType, class UnaryOpType>
std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators<
                     InputIteratorType, OutputIteratorType>::value,
                 OutputIteratorType>
transform_inclusive_scan(const ExecutionSpace& ex, InputIteratorType first,
                         InputIteratorType last, OutputIteratorType first_dest,
                         BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::transform_inclusive_scan_impl(
      "Kokkos::transform_inclusive_scan_custom_functors_iterator_api", ex,
      first, last, first_dest, binary_op, unary_op);
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class BinaryOpType, class UnaryOpType>
std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators<
                     InputIteratorType, OutputIteratorType>::value,
                 OutputIteratorType>
transform_inclusive_scan(const std::string& label, const ExecutionSpace& ex,
                         InputIteratorType first, InputIteratorType last,
                         OutputIteratorType first_dest, BinaryOpType binary_op,
                         UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::transform_inclusive_scan_impl(label, ex, first, last, first_dest,
                                             binary_op, unary_op);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryOpType,
          class UnaryOpType>
auto transform_inclusive_scan(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_impl(
      "Kokkos::transform_inclusive_scan_custom_functors_view_api", ex,
      KE::cbegin(view_from), KE::cend(view_from), KE::begin(view_dest),
      binary_op, unary_op);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryOpType,
          class UnaryOpType>
auto transform_inclusive_scan(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), binary_op, unary_op);
}

// overload set 2 (init value)
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class BinaryOpType, class UnaryOpType,
          class ValueType>
std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators<
                     InputIteratorType, OutputIteratorType>::value,
                 OutputIteratorType>
transform_inclusive_scan(const ExecutionSpace& ex, InputIteratorType first,
                         InputIteratorType last, OutputIteratorType first_dest,
                         BinaryOpType binary_op, UnaryOpType unary_op,
                         ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(ex);
  return Impl::transform_inclusive_scan_impl(
      "Kokkos::transform_inclusive_scan_custom_functors_iterator_api", ex,
      first, last, first_dest, binary_op, unary_op, init_value);
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class BinaryOpType, class UnaryOpType,
          class ValueType>
std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators<
                     InputIteratorType, OutputIteratorType>::value,
                 OutputIteratorType>
transform_inclusive_scan(const std::string& label, const ExecutionSpace& ex,
                         InputIteratorType first, InputIteratorType last,
                         OutputIteratorType first_dest, BinaryOpType binary_op,
                         UnaryOpType unary_op, ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(ex);
  return Impl::transform_inclusive_scan_impl(label, ex, first, last, first_dest,
                                             binary_op, unary_op, init_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryOpType,
          class UnaryOpType, class ValueType>
auto transform_inclusive_scan(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op, ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_impl(
      "Kokkos::transform_inclusive_scan_custom_functors_view_api", ex,
      KE::cbegin(view_from), KE::cend(view_from), KE::begin(view_dest),
      binary_op, unary_op, init_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class BinaryOpType,
          class UnaryOpType, class ValueType>
auto transform_inclusive_scan(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    BinaryOpType binary_op, UnaryOpType unary_op, ValueType init_value) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::transform_inclusive_scan_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), binary_op, unary_op, init_value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
