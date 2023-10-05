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

#ifndef KOKKOS_STD_ALGORITHMS_MINMAX_ELEMENT_HPP
#define KOKKOS_STD_ALGORITHMS_MINMAX_ELEMENT_HPP

#include "impl/Kokkos_MinMaxMinmaxElement.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType>
auto minmax_element(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last) {
  return Impl::minmax_element_impl<MinMaxFirstLastLoc>(
      "Kokkos::minmax_element_iterator_api_default", ex, first, last);
}

template <class ExecutionSpace, class IteratorType>
auto minmax_element(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last) {
  return Impl::minmax_element_impl<MinMaxFirstLastLoc>(label, ex, first, last);
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
auto minmax_element(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last, ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::minmax_element_impl<MinMaxFirstLastLocCustomComparator>(
      "Kokkos::minmax_element_iterator_api_default", ex, first, last,
      std::move(comp));
}

template <class ExecutionSpace, class IteratorType, class ComparatorType>
auto minmax_element(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last,
                    ComparatorType comp) {
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::minmax_element_impl<MinMaxFirstLastLocCustomComparator>(
      label, ex, first, last, std::move(comp));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto minmax_element(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::minmax_element_impl<MinMaxFirstLastLoc>(
      "Kokkos::minmax_element_view_api_default", ex, begin(v), end(v));
}

template <class ExecutionSpace, class DataType, class... Properties>
auto minmax_element(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::minmax_element_impl<MinMaxFirstLastLoc>(label, ex, begin(v),
                                                       end(v));
}

template <class ExecutionSpace, class DataType, class ComparatorType,
          class... Properties>
auto minmax_element(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v,
                    ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::minmax_element_impl<MinMaxFirstLastLocCustomComparator>(
      "Kokkos::minmax_element_view_api_default", ex, begin(v), end(v),
      std::move(comp));
}

template <class ExecutionSpace, class DataType, class ComparatorType,
          class... Properties>
auto minmax_element(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v,
                    ComparatorType comp) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  Impl::static_assert_is_not_openmptarget(ex);

  return Impl::minmax_element_impl<MinMaxFirstLastLocCustomComparator>(
      label, ex, begin(v), end(v), std::move(comp));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
