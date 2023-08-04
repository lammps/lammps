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

#ifndef KOKKOS_STD_ALGORITHMS_TRANSFORM_REDUCE_HPP
#define KOKKOS_STD_ALGORITHMS_TRANSFORM_REDUCE_HPP

#include "impl/Kokkos_TransformReduce.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

// ----------------------------
// overload set1:
// no custom functors passed, so equivalent to
// transform_reduce(first1, last1, first2, init, plus<>(), multiplies<>());
// ----------------------------
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType>
ValueType transform_reduce(const ExecutionSpace& ex, IteratorType1 first1,
                           IteratorType1 last1, IteratorType2 first2,
                           ValueType init_reduction_value) {
  return Impl::transform_reduce_default_functors_impl(
      "Kokkos::transform_reduce_default_functors_iterator_api", ex, first1,
      last1, first2, std::move(init_reduction_value));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType>
ValueType transform_reduce(const std::string& label, const ExecutionSpace& ex,
                           IteratorType1 first1, IteratorType1 last1,
                           IteratorType2 first2,
                           ValueType init_reduction_value) {
  return Impl::transform_reduce_default_functors_impl(
      label, ex, first1, last1, first2, std::move(init_reduction_value));
}

// overload1 accepting views
template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
ValueType transform_reduce(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& first_view,
    const ::Kokkos::View<DataType2, Properties2...>& second_view,
    ValueType init_reduction_value) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(first_view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(second_view);

  return Impl::transform_reduce_default_functors_impl(
      "Kokkos::transform_reduce_default_functors_iterator_api", ex,
      KE::cbegin(first_view), KE::cend(first_view), KE::cbegin(second_view),
      std::move(init_reduction_value));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
ValueType transform_reduce(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& first_view,
    const ::Kokkos::View<DataType2, Properties2...>& second_view,
    ValueType init_reduction_value) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(first_view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(second_view);

  return Impl::transform_reduce_default_functors_impl(
      label, ex, KE::cbegin(first_view), KE::cend(first_view),
      KE::cbegin(second_view), std::move(init_reduction_value));
}

//
// overload set2:
// accepts a custom transform and joiner functor
//

// Note the std refers to the arg BinaryReductionOp
// but in the Kokkos naming convention, it corresponds
// to a "joiner" that knows how to join two values
// NOTE: "joiner/transformer" need to be commutative.

// https://en.cppreference.com/w/cpp/algorithm/transform_reduce

// api accepting iterators
template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType, class BinaryJoinerType, class BinaryTransform>
ValueType transform_reduce(const ExecutionSpace& ex, IteratorType1 first1,
                           IteratorType1 last1, IteratorType2 first2,
                           ValueType init_reduction_value,
                           BinaryJoinerType joiner,
                           BinaryTransform transformer) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::transform_reduce_custom_functors_impl(
      "Kokkos::transform_reduce_custom_functors_iterator_api", ex, first1,
      last1, first2, std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2,
          class ValueType, class BinaryJoinerType, class BinaryTransform>
ValueType transform_reduce(const std::string& label, const ExecutionSpace& ex,
                           IteratorType1 first1, IteratorType1 last1,
                           IteratorType2 first2, ValueType init_reduction_value,
                           BinaryJoinerType joiner,
                           BinaryTransform transformer) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::transform_reduce_custom_functors_impl(
      label, ex, first1, last1, first2, std::move(init_reduction_value),
      std::move(joiner), std::move(transformer));
}

// accepting views
template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryJoinerType, class BinaryTransform>
ValueType transform_reduce(
    const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& first_view,
    const ::Kokkos::View<DataType2, Properties2...>& second_view,
    ValueType init_reduction_value, BinaryJoinerType joiner,
    BinaryTransform transformer) {
  namespace KE = ::Kokkos::Experimental;
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(first_view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(second_view);

  return Impl::transform_reduce_custom_functors_impl(
      "Kokkos::transform_reduce_custom_functors_view_api", ex,
      KE::cbegin(first_view), KE::cend(first_view), KE::cbegin(second_view),
      std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryJoinerType, class BinaryTransform>
ValueType transform_reduce(
    const std::string& label, const ExecutionSpace& ex,
    const ::Kokkos::View<DataType1, Properties1...>& first_view,
    const ::Kokkos::View<DataType2, Properties2...>& second_view,
    ValueType init_reduction_value, BinaryJoinerType joiner,
    BinaryTransform transformer) {
  namespace KE = ::Kokkos::Experimental;
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(first_view);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(second_view);

  return Impl::transform_reduce_custom_functors_impl(
      label, ex, KE::cbegin(first_view), KE::cend(first_view),
      KE::cbegin(second_view), std::move(init_reduction_value),
      std::move(joiner), std::move(transformer));
}

//
// overload set3:
//
// accepting iterators
template <class ExecutionSpace, class IteratorType, class ValueType,
          class BinaryJoinerType, class UnaryTransform>
// need this to avoid ambiguous call
std::enable_if_t<
    ::Kokkos::Experimental::Impl::are_iterators<IteratorType>::value, ValueType>
transform_reduce(const ExecutionSpace& ex, IteratorType first1,
                 IteratorType last1, ValueType init_reduction_value,
                 BinaryJoinerType joiner, UnaryTransform transformer) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::transform_reduce_custom_functors_impl(
      "Kokkos::transform_reduce_custom_functors_iterator_api", ex, first1,
      last1, std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

template <class ExecutionSpace, class IteratorType, class ValueType,
          class BinaryJoinerType, class UnaryTransform>
// need this to avoid ambiguous call
std::enable_if_t<
    ::Kokkos::Experimental::Impl::are_iterators<IteratorType>::value, ValueType>
transform_reduce(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first1, IteratorType last1,
                 ValueType init_reduction_value, BinaryJoinerType joiner,
                 UnaryTransform transformer) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::transform_reduce_custom_functors_impl(
      label, ex, first1, last1, std::move(init_reduction_value),
      std::move(joiner), std::move(transformer));
}

// accepting views
template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType, class BinaryJoinerType, class UnaryTransform>
ValueType transform_reduce(const ExecutionSpace& ex,
                           const ::Kokkos::View<DataType, Properties...>& view,
                           ValueType init_reduction_value,
                           BinaryJoinerType joiner,
                           UnaryTransform transformer) {
  namespace KE = ::Kokkos::Experimental;
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::transform_reduce_custom_functors_impl(
      "Kokkos::transform_reduce_custom_functors_view_api", ex, KE::cbegin(view),
      KE::cend(view), std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType, class BinaryJoinerType, class UnaryTransform>
ValueType transform_reduce(const std::string& label, const ExecutionSpace& ex,
                           const ::Kokkos::View<DataType, Properties...>& view,
                           ValueType init_reduction_value,
                           BinaryJoinerType joiner,
                           UnaryTransform transformer) {
  namespace KE = ::Kokkos::Experimental;
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::transform_reduce_custom_functors_impl(
      label, ex, KE::cbegin(view), KE::cend(view),
      std::move(init_reduction_value), std::move(joiner),
      std::move(transformer));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
