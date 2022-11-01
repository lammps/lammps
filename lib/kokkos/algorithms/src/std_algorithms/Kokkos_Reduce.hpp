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

#ifndef KOKKOS_STD_ALGORITHMS_REDUCE_HPP
#define KOKKOS_STD_ALGORITHMS_REDUCE_HPP

#include "impl/Kokkos_Reduce.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set 1
//
template <class ExecutionSpace, class IteratorType>
typename IteratorType::value_type reduce(const ExecutionSpace& ex,
                                         IteratorType first,
                                         IteratorType last) {
  return Impl::reduce_default_functors_impl(
      "Kokkos::reduce_default_functors_iterator_api", ex, first, last,
      typename IteratorType::value_type());
}

template <class ExecutionSpace, class IteratorType>
typename IteratorType::value_type reduce(const std::string& label,
                                         const ExecutionSpace& ex,
                                         IteratorType first,
                                         IteratorType last) {
  return Impl::reduce_default_functors_impl(
      label, ex, first, last, typename IteratorType::value_type());
}

template <class ExecutionSpace, class DataType, class... Properties>
auto reduce(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  using view_type  = ::Kokkos::View<DataType, Properties...>;
  using value_type = typename view_type::value_type;

  return Impl::reduce_default_functors_impl(
      "Kokkos::reduce_default_functors_view_api", ex, KE::cbegin(view),
      KE::cend(view), value_type());
}

template <class ExecutionSpace, class DataType, class... Properties>
auto reduce(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  using view_type  = ::Kokkos::View<DataType, Properties...>;
  using value_type = typename view_type::value_type;

  return Impl::reduce_default_functors_impl(label, ex, KE::cbegin(view),
                                            KE::cend(view), value_type());
}

//
// overload set2:
//
template <class ExecutionSpace, class IteratorType, class ValueType>
ValueType reduce(const ExecutionSpace& ex, IteratorType first,
                 IteratorType last, ValueType init_reduction_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::reduce_default_functors_impl(
      "Kokkos::reduce_default_functors_iterator_api", ex, first, last,
      init_reduction_value);
}

template <class ExecutionSpace, class IteratorType, class ValueType>
ValueType reduce(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last,
                 ValueType init_reduction_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::reduce_default_functors_impl(label, ex, first, last,
                                            init_reduction_value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType>
ValueType reduce(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 ValueType init_reduction_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::reduce_default_functors_impl(
      "Kokkos::reduce_default_functors_view_api", ex, KE::cbegin(view),
      KE::cend(view), init_reduction_value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType>
ValueType reduce(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 ValueType init_reduction_value) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::reduce_default_functors_impl(
      label, ex, KE::cbegin(view), KE::cend(view), init_reduction_value);
}

//
// overload set 3
//
template <class ExecutionSpace, class IteratorType, class ValueType,
          class BinaryOp>
ValueType reduce(const ExecutionSpace& ex, IteratorType first,
                 IteratorType last, ValueType init_reduction_value,
                 BinaryOp joiner) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::reduce_custom_functors_impl(
      "Kokkos::reduce_default_functors_iterator_api", ex, first, last,
      init_reduction_value, joiner);
}

template <class ExecutionSpace, class IteratorType, class ValueType,
          class BinaryOp>
ValueType reduce(const std::string& label, const ExecutionSpace& ex,
                 IteratorType first, IteratorType last,
                 ValueType init_reduction_value, BinaryOp joiner) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  return Impl::reduce_custom_functors_impl(label, ex, first, last,
                                           init_reduction_value, joiner);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType, class BinaryOp>
ValueType reduce(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 ValueType init_reduction_value, BinaryOp joiner) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::reduce_custom_functors_impl(
      "Kokkos::reduce_custom_functors_view_api", ex, KE::cbegin(view),
      KE::cend(view), init_reduction_value, joiner);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType, class BinaryOp>
ValueType reduce(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType, Properties...>& view,
                 ValueType init_reduction_value, BinaryOp joiner) {
  static_assert(std::is_move_constructible<ValueType>::value,
                "ValueType must be move constructible.");

  namespace KE = ::Kokkos::Experimental;
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::reduce_custom_functors_impl(label, ex, KE::cbegin(view),
                                           KE::cend(view), init_reduction_value,
                                           joiner);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
