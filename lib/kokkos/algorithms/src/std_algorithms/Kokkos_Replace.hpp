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

#ifndef KOKKOS_STD_ALGORITHMS_REPLACE_HPP
#define KOKKOS_STD_ALGORITHMS_REPLACE_HPP

#include "impl/Kokkos_Replace.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename Iterator, typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void replace(const ExecutionSpace& ex, Iterator first, Iterator last,
             const ValueType& old_value, const ValueType& new_value) {
  Impl::replace_exespace_impl("Kokkos::replace_iterator_api", ex, first, last,
                              old_value, new_value);
}

template <
    typename ExecutionSpace, typename Iterator, typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void replace(const std::string& label, const ExecutionSpace& ex, Iterator first,
             Iterator last, const ValueType& old_value,
             const ValueType& new_value) {
  Impl::replace_exespace_impl(label, ex, first, last, old_value, new_value);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void replace(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& view,
             const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  Impl::replace_exespace_impl("Kokkos::replace_view_api", ex, KE::begin(view),
                              KE::end(view), old_value, new_value);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void replace(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType1, Properties1...>& view,
             const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  Impl::replace_exespace_impl(label, ex, KE::begin(view), KE::end(view),
                              old_value, new_value);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename Iterator, typename ValueType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION void replace(const TeamHandleType& teamHandle, Iterator first,
                             Iterator last, const ValueType& old_value,
                             const ValueType& new_value) {
  Impl::replace_team_impl(teamHandle, first, last, old_value, new_value);
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename ValueType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION void replace(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view,
    const ValueType& old_value, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  Impl::replace_team_impl(teamHandle, KE::begin(view), KE::end(view), old_value,
                          new_value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
