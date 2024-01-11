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

//
// overload set accepting execution space
//
template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType2 swap_ranges(const ExecutionSpace& ex, IteratorType1 first1,
                          IteratorType1 last1, IteratorType2 first2) {
  return Impl::swap_ranges_exespace_impl(
      "Kokkos::swap_ranges_iterator_api_default", ex, first1, last1, first2);
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
auto swap_ranges(const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  assert(source.extent(0) == dest.extent(0));
  return Impl::swap_ranges_exespace_impl("Kokkos::swap_ranges_view_api_default",
                                         ex, begin(source), end(source),
                                         begin(dest));
}

template <typename ExecutionSpace, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
IteratorType2 swap_ranges(const std::string& label, const ExecutionSpace& ex,
                          IteratorType1 first1, IteratorType1 last1,
                          IteratorType2 first2) {
  return Impl::swap_ranges_exespace_impl(label, ex, first1, last1, first2);
}

template <typename ExecutionSpace, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
auto swap_ranges(const std::string& label, const ExecutionSpace& ex,
                 const ::Kokkos::View<DataType1, Properties1...>& source,
                 ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  assert(source.extent(0) == dest.extent(0));
  return Impl::swap_ranges_exespace_impl(label, ex, begin(source), end(source),
                                         begin(dest));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType1,
          typename IteratorType2,
          std::enable_if_t<is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION IteratorType2 swap_ranges(const TeamHandleType& teamHandle,
                                          IteratorType1 first1,
                                          IteratorType1 last1,
                                          IteratorType2 first2) {
  return Impl::swap_ranges_team_impl(teamHandle, first1, last1, first2);
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2,
          std::enable_if_t<is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto swap_ranges(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& source,
    ::Kokkos::View<DataType2, Properties2...>& dest) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(source);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(dest);

  assert(source.extent(0) == dest.extent(0));
  return Impl::swap_ranges_team_impl(teamHandle, begin(source), end(source),
                                     begin(dest));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
