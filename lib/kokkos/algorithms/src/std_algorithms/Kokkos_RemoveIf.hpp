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

#ifndef KOKKOS_STD_ALGORITHMS_REMOVE_IF_HPP
#define KOKKOS_STD_ALGORITHMS_REMOVE_IF_HPP

#include "impl/Kokkos_RemoveAllVariants.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <
    typename ExecutionSpace, typename Iterator, typename UnaryPredicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
Iterator remove_if(const ExecutionSpace& ex, Iterator first, Iterator last,
                   UnaryPredicate pred) {
  return Impl::remove_if_exespace_impl("Kokkos::remove_if_iterator_api_default",
                                       ex, first, last, pred);
}

template <
    typename ExecutionSpace, typename Iterator, typename UnaryPredicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
Iterator remove_if(const std::string& label, const ExecutionSpace& ex,
                   Iterator first, Iterator last, UnaryPredicate pred) {
  return Impl::remove_if_exespace_impl(label, ex, first, last, pred);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename UnaryPredicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto remove_if(const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view,
               UnaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::remove_if_exespace_impl("Kokkos::remove_if_iterator_api_default",
                                       ex, ::Kokkos::Experimental::begin(view),
                                       ::Kokkos::Experimental::end(view), pred);
}

template <
    typename ExecutionSpace, typename DataType, typename... Properties,
    typename UnaryPredicate,
    std::enable_if_t<::Kokkos::is_execution_space_v<ExecutionSpace>, int> = 0>
auto remove_if(const std::string& label, const ExecutionSpace& ex,
               const ::Kokkos::View<DataType, Properties...>& view,
               UnaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  return Impl::remove_if_exespace_impl(label, ex,
                                       ::Kokkos::Experimental::begin(view),
                                       ::Kokkos::Experimental::end(view), pred);
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename Iterator, typename UnaryPredicate,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION Iterator remove_if(const TeamHandleType& teamHandle,
                                   Iterator first, Iterator last,
                                   UnaryPredicate pred) {
  return Impl::remove_if_team_impl(teamHandle, first, last, pred);
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          typename UnaryPredicate,
          std::enable_if_t<::Kokkos::is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION auto remove_if(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType, Properties...>& view, UnaryPredicate pred) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  return Impl::remove_if_team_impl(teamHandle,
                                   ::Kokkos::Experimental::begin(view),
                                   ::Kokkos::Experimental::end(view), pred);
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
