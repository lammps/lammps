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

#ifndef KOKKOS_STD_ALGORITHMS_IS_PARTITIONED_HPP
#define KOKKOS_STD_ALGORITHMS_IS_PARTITIONED_HPP

#include "impl/Kokkos_IsPartitioned.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType, class PredicateType>
bool is_partitioned(const ExecutionSpace& ex, IteratorType first,
                    IteratorType last, PredicateType p) {
  return Impl::is_partitioned_impl(
      "Kokkos::is_partitioned_iterator_api_default", ex, first, last,
      std::move(p));
}

template <class ExecutionSpace, class IteratorType, class PredicateType>
bool is_partitioned(const std::string& label, const ExecutionSpace& ex,
                    IteratorType first, IteratorType last, PredicateType p) {
  return Impl::is_partitioned_impl(label, ex, first, last, std::move(p));
}

template <class ExecutionSpace, class PredicateType, class DataType,
          class... Properties>
bool is_partitioned(const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v,
                    PredicateType p) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::is_partitioned_impl("Kokkos::is_partitioned_view_api_default",
                                   ex, cbegin(v), cend(v), std::move(p));
}

template <class ExecutionSpace, class PredicateType, class DataType,
          class... Properties>
bool is_partitioned(const std::string& label, const ExecutionSpace& ex,
                    const ::Kokkos::View<DataType, Properties...>& v,
                    PredicateType p) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  return Impl::is_partitioned_impl(label, ex, cbegin(v), cend(v), std::move(p));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
