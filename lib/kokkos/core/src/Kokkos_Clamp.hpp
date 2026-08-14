// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CLAMP_HPP
#define KOKKOS_CLAMP_HPP

#include <Kokkos_Macros.hpp>

namespace Kokkos {

template <class T>
constexpr KOKKOS_INLINE_FUNCTION const T& clamp(const T& value, const T& lo,
                                                const T& hi) {
  KOKKOS_EXPECTS(!(hi < lo));
  // Capturing the result of std::clamp by reference produces a dangling
  // reference if one of the parameters is a temporary and that parameter is
  // returned.
  // NOLINTNEXTLINE(bugprone-return-const-ref-from-parameter)
  return (value < lo) ? lo : (hi < value) ? hi : value;
}

template <class T, class ComparatorType>
constexpr KOKKOS_INLINE_FUNCTION const T& clamp(const T& value, const T& lo,
                                                const T& hi,
                                                ComparatorType comp) {
  KOKKOS_EXPECTS(!comp(hi, lo));
  // Capturing the result of std::clamp by reference produces a dangling
  // reference if one of the parameters is a temporary and that parameter is
  // returned.
  // NOLINTNEXTLINE(bugprone-return-const-ref-from-parameter)
  return comp(value, lo) ? lo : comp(hi, value) ? hi : value;
}

}  // namespace Kokkos

#endif
