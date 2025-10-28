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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_EXPERIMENTAL_MDSPAN_HPP
#define KOKKOS_EXPERIMENTAL_MDSPAN_HPP

// Opt in for Kokkos::pair to submdspan/subview
// submdspan does only take index_pair_like which is derived from tuple_like
// tuple_like is an enumerated list:
// tuple, pair, array, complex, ranges::subrange
// Needs to be defined before including mdspan header

#include <Kokkos_Pair.hpp>

namespace Kokkos {
namespace detail {
template <class IdxT1, class IdxT2>
KOKKOS_INLINE_FUNCTION constexpr auto first_of(
    const pair<IdxT1, IdxT2> &slice) {
  return slice.first;
}
template <class IdxT1, class IdxT2, class Extents, size_t k>
KOKKOS_INLINE_FUNCTION constexpr auto last_of(std::integral_constant<size_t, k>,
                                              const Extents &,
                                              const pair<IdxT1, IdxT2> &slice) {
  return slice.second;
}

template <class T, class IndexType>
struct index_pair_like;

template <class IdxT1, class IdxT2, class IndexType>
struct index_pair_like<Kokkos::pair<IdxT1, IdxT2>, IndexType> {
  static constexpr bool value = std::is_convertible_v<IdxT1, IndexType> &&
                                std::is_convertible_v<IdxT2, IndexType>;
};
}  // namespace detail
}  // namespace Kokkos
#include <mdspan/mdspan.hpp>

#endif  // KOKKOS_EXPERIMENTAL_MDSPAN_HPP
