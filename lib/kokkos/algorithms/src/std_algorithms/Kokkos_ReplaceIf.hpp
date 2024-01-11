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

#ifndef KOKKOS_STD_ALGORITHMS_REPLACE_IF_HPP
#define KOKKOS_STD_ALGORITHMS_REPLACE_IF_HPP

#include "impl/Kokkos_ReplaceIf.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename InputIterator, typename Predicate,
    typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void replace_if(const ExecutionSpace& ex, InputIterator first,
                InputIterator last, Predicate pred,
                const ValueType& new_value) {
  Impl::replace_if_exespace_impl("Kokkos::replace_if_iterator_api", ex, first,
                                 last, pred, new_value);
}

template <
    typename ExecutionSpace, typename InputIterator, typename Predicate,
    typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void replace_if(const std::string& label, const ExecutionSpace& ex,
                InputIterator first, InputIterator last, Predicate pred,
                const ValueType& new_value) {
  Impl::replace_if_exespace_impl(label, ex, first, last, pred, new_value);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename Predicate, typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void replace_if(const ExecutionSpace& ex,
                const ::Kokkos::View<DataType1, Properties1...>& view,
                Predicate pred, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  Impl::replace_if_exespace_impl("Kokkos::replace_if_view_api", ex,
                                 KE::begin(view), KE::end(view), pred,
                                 new_value);
}

template <
    typename ExecutionSpace, typename DataType1, typename... Properties1,
    typename Predicate, typename ValueType,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
void replace_if(const std::string& label, const ExecutionSpace& ex,
                const ::Kokkos::View<DataType1, Properties1...>& view,
                Predicate pred, const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  Impl::replace_if_exespace_impl(label, ex, KE::begin(view), KE::end(view),
                                 pred, new_value);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename InputIterator, typename Predicate,
          typename ValueType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION void replace_if(const TeamHandleType& teamHandle,
                                InputIterator first, InputIterator last,
                                Predicate pred, const ValueType& new_value) {
  Impl::replace_if_team_impl(teamHandle, first, last, pred, new_value);
}

template <typename TeamHandleType, typename DataType1, typename... Properties1,
          typename Predicate, typename ValueType,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION void replace_if(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType1, Properties1...>& view, Predicate pred,
    const ValueType& new_value) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  namespace KE = ::Kokkos::Experimental;
  Impl::replace_if_team_impl(teamHandle, KE::begin(view), KE::end(view), pred,
                             new_value);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
