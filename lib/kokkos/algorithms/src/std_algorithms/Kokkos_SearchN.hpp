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

#ifndef KOKKOS_STD_ALGORITHMS_SEARCH_N_HPP
#define KOKKOS_STD_ALGORITHMS_SEARCH_N_HPP

#include "impl/Kokkos_SearchN.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

// overload set 1: no binary predicate passed
template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType>
IteratorType search_n(const ExecutionSpace& ex, IteratorType first,
                      IteratorType last, SizeType count,
                      const ValueType& value) {
  return Impl::search_n_impl("Kokkos::search_n_iterator_api_default", ex, first,
                             last, count, value);
}

template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType>
IteratorType search_n(const std::string& label, const ExecutionSpace& ex,
                      IteratorType first, IteratorType last, SizeType count,
                      const ValueType& value) {
  return Impl::search_n_impl(label, ex, first, last, count, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class ValueType>
auto search_n(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              SizeType count, const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_n_impl("Kokkos::search_n_view_api_default", ex,
                             KE::begin(view), KE::end(view), count, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class ValueType>
auto search_n(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              SizeType count, const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_n_impl(label, ex, KE::begin(view), KE::end(view), count,
                             value);
}

// overload set 2: binary predicate passed
template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType, class BinaryPredicateType>
IteratorType search_n(const ExecutionSpace& ex, IteratorType first,
                      IteratorType last, SizeType count, const ValueType& value,
                      const BinaryPredicateType& pred) {
  return Impl::search_n_impl("Kokkos::search_n_iterator_api_default", ex, first,
                             last, count, value, pred);
}

template <class ExecutionSpace, class IteratorType, class SizeType,
          class ValueType, class BinaryPredicateType>
IteratorType search_n(const std::string& label, const ExecutionSpace& ex,
                      IteratorType first, IteratorType last, SizeType count,
                      const ValueType& value, const BinaryPredicateType& pred) {
  return Impl::search_n_impl(label, ex, first, last, count, value, pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class ValueType, class BinaryPredicateType>
auto search_n(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              SizeType count, const ValueType& value,
              const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_n_impl("Kokkos::search_n_view_api_default", ex,
                             KE::begin(view), KE::end(view), count, value,
                             pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class SizeType, class ValueType, class BinaryPredicateType>
auto search_n(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              SizeType count, const ValueType& value,
              const BinaryPredicateType& pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::search_n_impl(label, ex, KE::begin(view), KE::end(view), count,
                             value, pred);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
