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

#ifndef KOKKOS_STD_ALGORITHMS_NONE_OF_HPP
#define KOKKOS_STD_ALGORITHMS_NONE_OF_HPP

#include "impl/Kokkos_AllOfAnyOfNoneOf.hpp"
#include "Kokkos_BeginEnd.hpp"

namespace Kokkos {
namespace Experimental {

template <class ExecutionSpace, class IteratorType, class Predicate>
bool none_of(const ExecutionSpace& ex, IteratorType first, IteratorType last,
             Predicate predicate) {
  return Impl::none_of_impl("Kokkos::none_of_iterator_api_default", ex, first,
                            last, predicate);
}

template <class ExecutionSpace, class IteratorType, class Predicate>
bool none_of(const std::string& label, const ExecutionSpace& ex,
             IteratorType first, IteratorType last, Predicate predicate) {
  return Impl::none_of_impl(label, ex, first, last, predicate);
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
bool none_of(const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::none_of_impl("Kokkos::none_of_view_api_default", ex,
                            KE::cbegin(v), KE::cend(v), std::move(predicate));
}

template <class ExecutionSpace, class DataType, class... Properties,
          class Predicate>
bool none_of(const std::string& label, const ExecutionSpace& ex,
             const ::Kokkos::View<DataType, Properties...>& v,
             Predicate predicate) {
  Impl::static_assert_is_admissible_to_kokkos_std_algorithms(v);

  namespace KE = ::Kokkos::Experimental;
  return Impl::none_of_impl(label, ex, KE::cbegin(v), KE::cend(v),
                            std::move(predicate));
}

}  // namespace Experimental
}  // namespace Kokkos

#endif
