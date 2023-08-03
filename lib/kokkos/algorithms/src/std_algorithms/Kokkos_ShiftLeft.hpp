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

#ifndef KOKKOS_STD_ALGORITHMS_SHIFT_LEFT_HPP
#define KOKKOS_STD_ALGORITHMS_SHIFT_LEFT_HPP

#include "impl/Kokkos_ShiftLeft.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType>
IteratorType shift_left(const ExecutionSpace& ex, IteratorType first,
                        IteratorType last,
                        typename IteratorType::difference_type n) {
  return Impl::shift_left_impl("Kokkos::shift_left_iterator_api_default", ex,
                               first, last, n);
}

template <class ExecutionSpace, class IteratorType>
IteratorType shift_left(const std::string& label, const ExecutionSpace& ex,
                        IteratorType first, IteratorType last,
                        typename IteratorType::difference_type n) {
  return Impl::shift_left_impl(label, ex, first, last, n);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto shift_left(const ExecutionSpace& ex,
                const ::Kokkos::View<DataType, Properties...>& view,
                typename decltype(begin(view))::difference_type n) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::shift_left_impl("Kokkos::shift_left_view_api_default", ex,
                               begin(view), end(view), n);
}

template <class ExecutionSpace, class DataType, class... Properties>
auto shift_left(const std::string& label, const ExecutionSpace& ex,
                const ::Kokkos::View<DataType, Properties...>& view,
                typename decltype(begin(view))::difference_type n) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::shift_left_impl(label, ex, begin(view), end(view), n);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
