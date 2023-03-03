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

#ifndef KOKKOS_STD_ALGORITHMS_FIND_IF_HPP
#define KOKKOS_STD_ALGORITHMS_FIND_IF_HPP

#include "impl/Kokkos_FindIfOrNot.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType find_if(const ExecutionSpace& ex, IteratorType first,
                     IteratorType last, PredicateType predicate) {
  return Impl::find_if_or_not_impl<true>("Kokkos::find_if_iterator_api_default",
                                         ex, first, last, std::move(predicate));
}

template <class ExecutionSpace, class IteratorType, class PredicateType>
IteratorType find_if(const std::string& label, const ExecutionSpace& ex,
                     IteratorType first, IteratorType last,
                     PredicateType predicate) {
  return Impl::find_if_or_not_impl<true>(label, ex, first, last,
                                         std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
auto find_if(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  namespace KE = ::Kokkos::Experimental;
  return Impl::find_if_or_not_impl<true>("Kokkos::find_if_view_api_default", ex,
                                         KE::begin(v), KE::end(v),
                                         std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
auto find_if(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);
  namespace KE = ::Kokkos::Experimental;
  return Impl::find_if_or_not_impl<true>(label, ex, KE::begin(v), KE::end(v),
                                         std::move(predicate));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
