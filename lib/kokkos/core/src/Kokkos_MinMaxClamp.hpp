/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_MIN_MAX_CLAMP_HPP
#define KOKKOS_MIN_MAX_CLAMP_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Pair.hpp>

#include <initializer_list>

namespace Kokkos {
namespace Experimental {

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

}  // namespace Experimental
}  // namespace Kokkos

#endif
