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

// note: the enable_if below is to avoid "call to ... is ambiguous"
// for example in the unit test when using a variadic function

// overload set1
template <class ExecutionSpace, class IteratorType>
std::enable_if_t<!::Kokkos::is_view<IteratorType>::value, IteratorType> unique(
    const ExecutionSpace& ex, IteratorType first, IteratorType last) {
  return Impl::unique_impl("Kokkos::unique_iterator_api_default", ex, first,
                           last);
}

template <class ExecutionSpace, class IteratorType>
std::enable_if_t<!::Kokkos::is_view<IteratorType>::value, IteratorType> unique(
    const std::string& label, const ExecutionSpace& ex, IteratorType first,
    IteratorType last) {
  return Impl::unique_impl(label, ex, first, last);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto unique(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return ::Kokkos::Experimental::unique("Kokkos::unique_view_api_default", ex,
                                        begin(view), end(view));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto unique(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return ::Kokkos::Experimental::unique(label, ex, begin(view), end(view));
}

// overload set2
template <class ExecutionSpace, class IteratorType, class BinaryPredicate>
IteratorType unique(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last, BinaryPredicate pred) {
  return Impl::unique_impl("Kokkos::unique_iterator_api_default", ex, first,
                           last, pred);
}

template <class ExecutionSpace, class IteratorType, class BinaryPredicate>
IteratorType unique(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last,
                    BinaryPredicate pred) {
  return Impl::unique_impl(label, ex, first, last, pred);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class BinaryPredicate>
auto unique(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            BinaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::unique_impl("Kokkos::unique_view_api_default", ex, begin(view),
                           end(view), std::move(pred));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class BinaryPredicate>
auto unique(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            BinaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::unique_impl(label, ex, begin(view), end(view), std::move(pred));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
