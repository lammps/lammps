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

#ifndef KOKKOS_STD_ALGORITHMS_SWAP_RANGES_HPP
#define KOKKOS_STD_ALGORITHMS_SWAP_RANGES_HPP

#include "impl/Kokkos_SwapRanges.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 swap_ranges(const ExecutionSpace& ex, IteratorType1 first1,
                          IteratorType1 last1, IteratorType2 first2) {
  return Impl::swap_ranges_impl("Kokkos::swap_ranges_iterator_api_default", ex,
                                first1, last1, first2);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto swap_ranges(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  assert(source.extent(0) == dest.extent(0));
  return Impl::swap_ranges_impl("Kokkos::swap_ranges_view_api_default", ex,
                                begin(source), end(source), begin(dest));
}

template <class ExecutionSpace, class IteratorType1, class IteratorType2>
IteratorType2 swap_ranges(const std::string& label, const ExecutionSpace& ex,
                          IteratorType1 first1, IteratorType1 last1,
                          IteratorType2 first2) {
  return Impl::swap_ranges_impl(label, ex, first1, last1, first2);
}

template <class ExecutionSpace, class DataType1, class... Properties1,
          class DataType2, class... Properties2>
auto swap_ranges(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  assert(source.extent(0) == dest.extent(0));
  return Impl::swap_ranges_impl(label, ex, begin(source), end(source),
                                begin(dest));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
