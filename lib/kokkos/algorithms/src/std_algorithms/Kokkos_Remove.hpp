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

#ifndef KOKKOS_STD_ALGORITHMS_REMOVE_HPP
#define KOKKOS_STD_ALGORITHMS_REMOVE_HPP

#include "impl/Kokkos_RemoveAllVariants.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class Iterator, class ValueType>
Iterator remove(const ExecutionSpace& ex, Iterator first, Iterator last,
                const ValueType& value) {
  return Impl::remove_impl("Kokkos::remove_iterator_api_default", ex, first,
                           last, value);
}

template <class ExecutionSpace, class Iterator, class ValueType>
Iterator remove(const std::string& label, const ExecutionSpace& ex,
                Iterator first, Iterator last, const ValueType& value) {
  return Impl::remove_impl(label, ex, first, last, value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType>
auto remove(const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::remove_impl("Kokkos::remove_iterator_api_default", ex,
                           ::Kokkos::Experimental::begin(view),
                           ::Kokkos::Experimental::end(view), value);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class ValueType>
auto remove(const std::string& label, const ExecutionSpace& ex,
            const ::Kokkos::View<DataType, Properties...>& view,
            const ValueType& value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::remove_impl(label, ex, ::Kokkos::Experimental::begin(view),
                           ::Kokkos::Experimental::end(view), value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
