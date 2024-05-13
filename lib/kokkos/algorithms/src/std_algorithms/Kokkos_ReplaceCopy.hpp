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

#ifndef KOKKOS_STD_ALGORITHMS_REPLACE_COPY_HPP
#define KOKKOS_STD_ALGORITHMS_REPLACE_COPY_HPP

#include "impl/Kokkos_ReplaceCopy.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename InputIterator, typename OutputIterator,
    typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
OutputIterator replace_copy(const ExecutionSpace& ex, InputIterator first_from,
                            InputIterator last_from, OutputIterator first_dest,
                            const ValueType& old_value,
                            const ValueType& new_value) {
  return Impl::replace_copy_exespace_impl("Kokkos::replace_copy_iterator_api",
                                          ex, first_from, last_from, first_dest,
                                          old_value, new_value);
}

template <
    typename ExecutionSpace, typename InputIterator, typename OutputIterator,
    typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
OutputIterator replace_copy(const std::string& label, const ExecutionSpace& ex,
                            InputIterator first_from, InputIterator last_from,
                            OutputIterator first_dest,
                            const ValueType& old_value,
                            const ValueType& new_value) {
  return Impl::replace_copy_exespace_impl(label, ex, first_from, last_from,
                                          first_dest, old_value, new_value);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto replace_copy(const ExecutionSpace& ex,
                  const ::Kokkos::View<DataType1, Properties1...>& view_from,
                  const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                  const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_copy_exespace_impl(
      "Kokkos::replace_copy_view_api", ex, KE::cbegin(view_from),
      KE::cend(view_from), KE::begin(view_dest), old_value, new_value);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename DataType2, typename... Properties2, typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto replace_copy(const std::string& label, const ExecutionSpace& ex,
                  const ::Kokkos::View<DataType1, Properties1...>& view_from,
                  const ::Kokkos::View<DataType2, Properties2...>& view_dest,
                  const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_copy_exespace_impl(
      label, ex, KE::cbegin(view_from), KE::cend(view_from),
      KE::begin(view_dest), old_value, new_value);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename InputIterator,
          typename OutputIterator, typename ValueType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION OutputIterator replace_copy(const TeamHandleType& teamHandle,
                                            InputIterator first_from,
                                            InputIterator last_from,
                                            OutputIterator first_dest,
                                            const ValueType& old_value,
                                            const ValueType& new_value) {
  return Impl::replace_copy_team_impl(teamHandle, first_from, last_from,
                                      first_dest, old_value, new_value);
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename DataType2, typename... Properties2, typename ValueType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto replace_copy(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view_from,
    const ::Kokkos::View<DataType2, Properties2...>& view_dest,
    const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_from);
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view_dest);
  namespace KE = ::Kokkos::Experimental;
  return Impl::replace_copy_team_impl(teamHandle, KE::cbegin(view_from),
                                      KE::cend(view_from), KE::begin(view_dest),
                                      old_value, new_value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
