// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

#ifdef KOKKOS_ENABLE_OPENMPTARGET
#include <cassert>
#include <iostream>
#endif

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

// FIXME_OPENMPTARGET We need to inject our own error handler as the default
// mdspan one cannot be called from device code
#ifdef KOKKOS_ENABLE_OPENMPTARGET
namespace Kokkos::detail {
KOKKOS_INLINE_FUNCTION void openmptarget_precondition_handler(const char *cond,
                                                              const char *file,
                                                              unsigned line) {
  ::printf("%s:%u: precondition failure: `%s`\n", file, line, cond);
  assert(0);
}
}  // namespace Kokkos::detail
#define MDSPAN_IMPL_PRECONDITION_VIOLATION_HANDLER(cond, file, line) \
  Kokkos::detail::openmptarget_precondition_handler(cond, file, line)
#endif

#include <mdspan/mdspan.hpp>

#endif  // KOKKOS_EXPERIMENTAL_MDSPAN_HPP
