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

#ifndef KOKKOS_STD_ALGORITHMS_GENERATE_HPP
#define KOKKOS_STD_ALGORITHMS_GENERATE_HPP

#include "impl/Kokkos_GenerateGenerateN.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

//
// overload set accepting execution space
//
template <typename ExecutionSpace, typename IteratorType, typename Generator,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
void generate(const ExecutionSpace& ex, IteratorType first, IteratorType last,
              Generator g) {
  Impl::generate_exespace_impl("Kokkos::generate_iterator_api_default", ex,
                               first, last, std::move(g));
}

template <typename ExecutionSpace, typename IteratorType, typename Generator,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
void generate(const std::string& label, const ExecutionSpace& ex,
              IteratorType first, IteratorType last, Generator g) {
  Impl::generate_exespace_impl(label, ex, first, last, std::move(g));
}

template <typename ExecutionSpace, typename DataType, typename... Properties,
          typename Generator,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
void generate(const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              Generator g) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  Impl::generate_exespace_impl("Kokkos::generate_view_api_default", ex,
                               begin(view), end(view), std::move(g));
}

template <typename ExecutionSpace, typename DataType, typename... Properties,
          typename Generator,
          std::enable_if_t<is_execution_space_v<ExecutionSpace>, int> = 0>
void generate(const std::string& label, const ExecutionSpace& ex,
              const ::Kokkos::View<DataType, Properties...>& view,
              Generator g) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);

  Impl::generate_exespace_impl(label, ex, begin(view), end(view), std::move(g));
}

//
// overload set accepting a team handle
// Note: for now omit the overloads accepting a label
// since they cause issues on device because of the string allocation.
//
template <typename TeamHandleType, typename IteratorType, typename Generator,
          std::enable_if_t<is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION void generate(const TeamHandleType& teamHandle,
                              IteratorType first, IteratorType last,
                              Generator g) {
  Impl::generate_team_impl(teamHandle, first, last, std::move(g));
}

template <typename TeamHandleType, typename DataType, typename... Properties,
          typename Generator,
          std::enable_if_t<is_team_handle_v<TeamHandleType>, int> = 0>
KOKKOS_FUNCTION void generate(
    const TeamHandleType& teamHandle,
    const ::Kokkos::View<DataType, Properties...>& view, Generator g) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(view);
  Impl::generate_team_impl(teamHandle, begin(view), end(view), std::move(g));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
