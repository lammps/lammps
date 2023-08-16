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

#ifndef KOKKOS_STD_ALGORITHMS_FIND_HPP
#define KOKKOS_STD_ALGORITHMS_FIND_HPP

#include "impl/Kokkos_FindIfOrNot.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class InputIterator, class T>
InputIterator find(const ExecutionSpace& ex, InputIterator first,
                   InputIterator last, const T& value) {
  return Impl::find_impl("Kokkos::find_iterator_api_default", ex, first, last,
                         value);
}

template <class ExecutionSpace, class InputIterator, class T>
InputIterator find(const std::string& label, const ExecutionSpace& ex,
                   InputIterator first, InputIterator last, const T& value) {
  return Impl::find_impl(label, ex, first, last, value);
}

template <class ExecutionSpace, class DataType, class... Properties, class T>
auto find(const ExecutionSpace& ex,
          const ::Kokkos::View<DataType, Properties...>& view, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_impl("Kokkos::find_view_api_default", ex, KE::begin(view),
                         KE::end(view), value);
}

template <class ExecutionSpace, class DataType, class... Properties, class T>
auto find(const std::string& label, const ExecutionSpace& ex,
          const ::Kokkos::View<DataType, Properties...>& view, const T& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  namespace KE = ::Kokkos::Experimental;
  return Impl::find_impl(label, ex, KE::begin(view), KE::end(view), value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
