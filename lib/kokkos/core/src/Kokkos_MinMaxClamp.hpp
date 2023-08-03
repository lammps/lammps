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
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_MIN_MAX_CLAMP_HPP
#define KOKKOS_MIN_MAX_CLAMP_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Pair.hpp>

#include <initializer_list>

namespace Kokkos {

// clamp
template <class T>
constexpr KOKKOS_INLINE_FUNCTION const T& clamp(const T& value, const T& lo,
                                                const T& hi) {
  KOKKOS_EXPECTS(!(hi < lo));
  return (value < lo) ? lo : (hi < value) ? hi : value;
}

template <class T, class ComparatorType>
constexpr KOKKOS_INLINE_FUNCTION const T& clamp(const T& value, const T& lo,
                                                const T& hi,
                                                ComparatorType comp) {
  KOKKOS_EXPECTS(!comp(hi, lo));
  return comp(value, lo) ? lo : comp(hi, value) ? hi : value;
}

// max
template <class T>
constexpr KOKKOS_INLINE_FUNCTION const T& max(const T& a, const T& b) {
  return (a < b) ? b : a;
}

template <class T, class ComparatorType>
constexpr KOKKOS_INLINE_FUNCTION const T& max(const T& a, const T& b,
                                              ComparatorType comp) {
  return comp(a, b) ? b : a;
}

template <class T>
KOKKOS_INLINE_FUNCTION constexpr T max(std::initializer_list<T> ilist) {
  auto first      = ilist.begin();
  auto const last = ilist.end();
  auto result     = *first;
  if (first == last) return result;
  while (++first != last) {
    if (result < *first) result = *first;
  }
  return result;
}

template <class T, class Compare>
KOKKOS_INLINE_FUNCTION constexpr T max(std::initializer_list<T> ilist,
                                       Compare comp) {
  auto first      = ilist.begin();
  auto const last = ilist.end();
  auto result     = *first;
  if (first == last) return result;
  while (++first != last) {
    if (comp(result, *first)) result = *first;
  }
  return result;
}

// min
template <class T>
constexpr KOKKOS_INLINE_FUNCTION const T& min(const T& a, const T& b) {
  return (b < a) ? b : a;
}

template <class T, class ComparatorType>
constexpr KOKKOS_INLINE_FUNCTION const T& min(const T& a, const T& b,
                                              ComparatorType comp) {
  return comp(b, a) ? b : a;
}

template <class T>
KOKKOS_INLINE_FUNCTION constexpr T min(std::initializer_list<T> ilist) {
  auto first      = ilist.begin();
  auto const last = ilist.end();
  auto result     = *first;
  if (first == last) return result;
  while (++first != last) {
    if (*first < result) result = *first;
  }
  return result;
}

template <class T, class Compare>
KOKKOS_INLINE_FUNCTION constexpr T min(std::initializer_list<T> ilist,
                                       Compare comp) {
  auto first      = ilist.begin();
  auto const last = ilist.end();
  auto result     = *first;
  if (first == last) return result;
  while (++first != last) {
    if (comp(*first, result)) result = *first;
  }
  return result;
}

// minmax
template <class T>
constexpr KOKKOS_INLINE_FUNCTION auto minmax(const T& a, const T& b) {
  using return_t = ::Kokkos::pair<const T&, const T&>;
  return (b < a) ? return_t{b, a} : return_t{a, b};
}

template <class T, class ComparatorType>
constexpr KOKKOS_INLINE_FUNCTION auto minmax(const T& a, const T& b,
                                             ComparatorType comp) {
  using return_t = ::Kokkos::pair<const T&, const T&>;
  return comp(b, a) ? return_t{b, a} : return_t{a, b};
}

template <class T>
KOKKOS_INLINE_FUNCTION constexpr Kokkos::pair<T, T> minmax(
    std::initializer_list<T> ilist) {
  auto first      = ilist.begin();
  auto const last = ilist.end();
  auto next       = first;
  Kokkos::pair<T, T> result{*first, *first};
  if (first == last || ++next == last) return result;
  if (*next < *first)
    result.first = *next;
  else
    result.second = *next;
  first = next;
  while (++first != last) {
    if (++next == last) {
      if (*first < result.first)
        result.first = *first;
      else if (!(*first < result.second))
        result.second = *first;
      break;
    }
    if (*next < *first) {
      if (*next < result.first) result.first = *next;
      if (!(*first < result.second)) result.second = *first;
    } else {
      if (*first < result.first) result.first = *first;
      if (!(*next < result.second)) result.second = *next;
    }
    first = next;
  }
  return result;
}

template <class T, class Compare>
KOKKOS_INLINE_FUNCTION constexpr Kokkos::pair<T, T> minmax(
    std::initializer_list<T> ilist, Compare comp) {
  auto first      = ilist.begin();
  auto const last = ilist.end();
  auto next       = first;
  Kokkos::pair<T, T> result{*first, *first};
  if (first == last || ++next == last) return result;
  if (comp(*next, *first))
    result.first = *next;
  else
    result.second = *next;
  first = next;
  while (++first != last) {
    if (++next == last) {
      if (comp(*first, result.first))
        result.first = *first;
      else if (!comp(*first, result.second))
        result.second = *first;
      break;
    }
    if (comp(*next, *first)) {
      if (comp(*next, result.first)) result.first = *next;
      if (!comp(*first, result.second)) result.second = *first;
    } else {
      if (comp(*first, result.first)) result.first = *first;
      if (!comp(*next, result.second)) result.second = *next;
    }
    first = next;
  }
  return result;
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
namespace Experimental {
using ::Kokkos::clamp;
using ::Kokkos::max;
using ::Kokkos::min;
using ::Kokkos::minmax;
}  // namespace Experimental
#endif

}  // namespace Kokkos

#endif
