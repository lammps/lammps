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

#ifndef KOKKOS_STD_ALGORITHMS_PARTITION_POINT_HPP
#define KOKKOS_STD_ALGORITHMS_PARTITION_POINT_HPP

#include "impl/Kokkos_PartitionPoint.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType, class UnaryPredicate>
IteratorType partition_point(const ExecutionSpace& ex, IteratorType first,
                             IteratorType last, UnaryPredicate p) {
  return Impl::partition_point_impl(
      "Kokkos::partitioned_point_iterator_api_default", ex, first, last,
      std::move(p));
}

template <class ExecutionSpace, class IteratorType, class UnaryPredicate>
IteratorType partition_point(const std::string& label, const ExecutionSpace& ex,
                             IteratorType first, IteratorType last,
                             UnaryPredicate p) {
  return Impl::partition_point_impl(label, ex, first, last, std::move(p));
}

template <class ExecutionSpace, class UnaryPredicate, class DataType,
          class... Properties>
auto partition_point(const std::string& label, const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType, Properties...>& v,
                     UnaryPredicate p) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  return Impl::partition_point_impl(label, ex, begin(v), end(v), std::move(p));
}

template <class ExecutionSpace, class UnaryPredicate, class DataType,
          class... Properties>
auto partition_point(const ExecutionSpace& ex,
                     const ::Kokkos::View<DataType, Properties...>& v,
                     UnaryPredicate p) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  return Impl::partition_point_impl("Kokkos::partition_point_view_api_default",
                                    ex, begin(v), end(v), std::move(p));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
