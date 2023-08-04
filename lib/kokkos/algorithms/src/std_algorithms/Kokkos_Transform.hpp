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

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class UnaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator, OutputIterator>::value,
                  OutputIterator>
transform(const ExecutionSpace& ex, InputIterator first1, InputIterator last1,
          OutputIterator d_first, UnaryOperation unary_op) {
  return Impl::transform_impl("Kokkos::transform_iterator_api_default", ex,
                              first1, last1, d_first, std::move(unary_op));
}

template <class ExecutionSpace, class InputIterator, class OutputIterator,
          class UnaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator, OutputIterator>::value,
                  OutputIterator>
transform(const std::string& label, const ExecutionSpace& ex,
          InputIterator first1, InputIterator last1, OutputIterator d_first,
          UnaryOperation unary_op) {
  return Impl::transform_impl(label, ex, first1, last1, d_first,
                              std::move(unary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class UnaryOperation>
auto transform(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source,
               ::Kokkos::View<DataType2, Properties2...>& dest,
               UnaryOperation unary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl("Kokkos::transform_view_api_default", ex,
                              begin(source), end(source), begin(dest),
                              std::move(unary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class UnaryOperation>
auto transform(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source,
               ::Kokkos::View<DataType2, Properties2...>& dest,
               UnaryOperation unary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl(label, ex, begin(source), end(source),
                              begin(dest), std::move(unary_op));
}

template <class ExecutionSpace, class InputIterator1, class InputIterator2,
          class OutputIterator, class BinaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator1, InputIterator2, OutputIterator>::value,
                  OutputIterator>
transform(const ExecutionSpace& ex, InputIterator1 first1, InputIterator1 last1,
          InputIterator2 first2, OutputIterator d_first,
          BinaryOperation binary_op) {
  return Impl::transform_impl("Kokkos::transform_iterator_api_default", ex,
                              first1, last1, first2, d_first,
                              std::move(binary_op));
}

template <class ExecutionSpace, class InputIterator1, class InputIterator2,
          class OutputIterator, class BinaryOperation>
std::enable_if_t< ::Kokkos::Experimental::Impl::are_iterators<
                      InputIterator1, InputIterator2, OutputIterator>::value,
                  OutputIterator>
transform(const std::string& label, const ExecutionSpace& ex,
          InputIterator1 first1, InputIterator1 last1, InputIterator2 first2,
          OutputIterator d_first, BinaryOperation binary_op) {
  return Impl::transform_impl(label, ex, first1, last1, first2, d_first,
                              std::move(binary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class BinaryOperation>
auto transform(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source1,
               const ::Kokkos::View<DataType2, Properties2...>& source2,
               ::Kokkos::View<DataType3, Properties3...>& dest,
               BinaryOperation binary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source2);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl("Kokkos::transform_view_api_default", ex,
                              begin(source1), end(source1), begin(source2),
                              begin(dest), std::move(binary_op));
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class DataType3,
          class... Properties3, class BinaryOperation>
auto transform(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType1, Properties1...>& source1,
               const ::Kokkos::View<DataType2, Properties2...>& source2,
               ::Kokkos::View<DataType3, Properties3...>& dest,
               BinaryOperation binary_op) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source1);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source2);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  return Impl::transform_impl(label, ex, begin(source1), end(source1),
                              begin(source2), begin(dest),
                              std::move(binary_op));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
