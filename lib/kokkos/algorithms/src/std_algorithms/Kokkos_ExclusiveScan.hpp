/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_STD_ALGORITHMS_EXCLUSIVE_SCAN_HPP
#define KOKKOS_STD_ALGORITHMS_EXCLUSIVE_SCAN_HPP

#include "impl/Kokkos_ExclusiveScan.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

// overload set 1
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType>
std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators<
                     InputIteratorType, OutputIteratorType>::value,
                 OutputIteratorType>
exclusive_scan(const ExecutionSpace& ex, InputIteratorType first,
               InputIteratorType last, OutputIteratorType first_dest,
               ValueType init_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::exclusive_scan_default_op_impl(
      "Kokkos::exclusive_scan_default_functors_iterator_api", ex, first, last,
      first_dest, init_value);
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType>
std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators<
                     InputIteratorType, OutputIteratorType>::value,
                 OutputIteratorType>
exclusive_scan(const std::string& label, const ExecutionSpace& ex,
               InputIteratorType first, InputIteratorType last,
               OutputIteratorType first_dest, ValueType init_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::exclusive_scan_default_op_impl(label, ex, first, last,
                                              first_dest, init_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
auto exclusive_scan(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    ValueType init_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::exclusive_scan_default_op_impl(
      "Kokkos::exclusive_scan_default_functors_view_api", ex,
      KE::cbegin(view_from), KE::cend(view_from), KE::begin(view_dest),
      init_value);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType>
auto exclusive_scan(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    ValueType init_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::exclusive_scan_default_op_impl(label, ex, KE::cbegin(view_from),
                                              KE::cend(view_from),
                                              KE::begin(view_dest), init_value);
}

// overload set 2
template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType>
std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators<
                     InputIteratorType, OutputIteratorType>::value,
                 OutputIteratorType>
exclusive_scan(const ExecutionSpace& ex, InputIteratorType first,
               InputIteratorType last, OutputIteratorType first_dest,
               ValueType init_value, BinaryOpType bop) {
  Impl::static_assert_is_not_openmptarget(ex);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::exclusive_scan_custom_op_impl(
      "Kokkos::exclusive_scan_custom_functors_iterator_api", ex, first, last,
      first_dest, init_value, bop);
}

template <class ExecutionSpace, class InputIteratorType,
          class OutputIteratorType, class ValueType, class BinaryOpType>
std::enable_if_t<::Kokkos::Experimental::Impl::are_iterators<
                     InputIteratorType, OutputIteratorType>::value,
                 OutputIteratorType>
exclusive_scan(const std::string& label, const ExecutionSpace& ex,
               InputIteratorType first, InputIteratorType last,
               OutputIteratorType first_dest, ValueType init_value,
               BinaryOpType bop) {
  Impl::static_assert_is_not_openmptarget(ex);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  return Impl::exclusive_scan_custom_op_impl(label, ex, first, last, first_dest,
                                             init_value, bop);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryOpType>
auto exclusive_scan(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    ValueType init_value, BinaryOpType bop) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::exclusive_scan_custom_op_impl(
      "Kokkos::exclusive_scan_custom_functors_view_api", ex,
      KE::cbegin(view_from), KE::cend(view_from), KE::begin(view_dest),
      init_value, bop);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2, class ValueType,
          class BinaryOpType>
auto exclusive_scan(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType1, Properties1...>& view_from,
                    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                    ValueType init_value, BinaryOpType bop) {
  Impl::static_assert_is_not_openmptarget(ex);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");
  namespace KE = ::Kokkos::Experimental;
  return Impl::exclusive_scan_custom_op_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), init_value, bop);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
