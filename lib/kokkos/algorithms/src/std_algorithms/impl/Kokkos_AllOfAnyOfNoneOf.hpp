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

#ifndef KOKKOS_STD_ALGORITHMS_ALL_OF_ANY_OF_NONE_OF_IMPL_HPP
#define KOKKOS_STD_ALGORITHMS_ALL_OF_ANY_OF_NONE_OF_IMPL_HPP

#include "Kokkos_FindIfOrNot.hpp"

namespace Kokkos {
namespace Experimental {
namespace Impl {

//
// exespace impl
//
template <class ExecutionSpace, class InputIterator, class Predicate>
bool all_of_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                          InputIterator first, InputIterator last,
                          Predicate predicate) {
  return (find_if_or_not_exespace_impl<false>(label, ex, first, last,
                                              predicate) == last);
}

template <class ExecutionSpace, class InputIterator, class Predicate>
bool any_of_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                          InputIterator first, InputIterator last,
                          Predicate predicate) {
  return (find_if_or_not_exespace_impl<true>(label, ex, first, last,
                                             predicate) != last);
}

template <class ExecutionSpace, class IteratorType, class Predicate>
bool none_of_exespace_impl(const std::string& label, const ExecutionSpace& ex,
                           IteratorType first, IteratorType last,
                           Predicate predicate) {
  return (find_if_or_not_exespace_impl<true>(label, ex, first, last,
                                             predicate) == last);
}

//
// team impl
//
template <class TeamHandleType, class InputIterator, class Predicate>
KOKKOS_FUNCTION bool all_of_team_impl(const TeamHandleType& teamHandle,
                                      InputIterator first, InputIterator last,
                                      Predicate predicate) {
  return (find_if_or_not_team_impl<false>(teamHandle, first, last, predicate) ==
          last);
}

template <class TeamHandleType, class InputIterator, class Predicate>
KOKKOS_FUNCTION bool any_of_team_impl(const TeamHandleType& teamHandle,
                                      InputIterator first, InputIterator last,
                                      Predicate predicate) {
  return (find_if_or_not_team_impl<true>(teamHandle, first, last, predicate) !=
          last);
}

template <class TeamHandleType, class IteratorType, class Predicate>
KOKKOS_FUNCTION bool none_of_team_impl(const TeamHandleType& teamHandle,
                                       IteratorType first, IteratorType last,
                                       Predicate predicate) {
  return (find_if_or_not_team_impl<true>(teamHandle, first, last, predicate) ==
          last);
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
