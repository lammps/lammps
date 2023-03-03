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
#ifndef KOKKOS_ATOMIC_WINDOWS_HPP
#define KOKKOS_ATOMIC_WINDOWS_HPP

#ifdef _WIN32

#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <winsock2.h>
#include <windows.h>

namespace Kokkos {
namespace Impl {
#ifdef _MSC_VER
_declspec(align(16))
#endif
    struct cas128_t {
  LONGLONG lower;
  LONGLONG upper;
  KOKKOS_INLINE_FUNCTION
  bool operator!=(const cas128_t& a) const {
    return (lower != a.lower) || upper != a.upper;
  }
}
#if defined(__GNUC__) || defined(__clang__)
__attribute__((aligned(16)))
#endif
;
}  // namespace Impl

#if !defined(__CUDA_ARCH__) || defined(__clang__)
template <typename T>
inline T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    std::enable_if_t<sizeof(T) == sizeof(CHAR), const T&> val) {
  union U {
    CHAR i;
    T t;
    KOKKOS_INLINE_FUNCTION U(){};
  } tmp;

  tmp.i = _InterlockedCompareExchange8((CHAR*)dest, *((CHAR*)&val),
                                       *((CHAR*)&compare));
  return tmp.t;
}

template <typename T>
inline T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    std::enable_if_t<sizeof(T) == sizeof(SHORT), const T&> val) {
  union U {
    SHORT i;
    T t;
    KOKKOS_INLINE_FUNCTION U(){};
  } tmp;

  tmp.i = _InterlockedCompareExchange16((SHORT*)dest, *((SHORT*)&val),
                                        *((SHORT*)&compare));
  return tmp.t;
}

template <typename T>
inline T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    std::enable_if_t<sizeof(T) == sizeof(LONG), const T&> val) {
  union U {
    LONG i;
    T t;
    KOKKOS_INLINE_FUNCTION U() {}
  } tmp;

  tmp.i = _InterlockedCompareExchange((LONG*)dest, *((LONG*)&val),
                                      *((LONG*)&compare));
  return tmp.t;
}

template <typename T>
inline T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    std::enable_if_t<sizeof(T) == sizeof(LONGLONG), const T&> val) {
  union U {
    LONGLONG i;
    T t;
    KOKKOS_INLINE_FUNCTION U() {}
  } tmp;

  tmp.i = _InterlockedCompareExchange64((LONGLONG*)dest, *((LONGLONG*)&val),
                                        *((LONGLONG*)&compare));
  return tmp.t;
}

template <typename T>
inline T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    std::enable_if_t<sizeof(T) == sizeof(Impl::cas128_t), const T&> val) {
  T compare_and_result(compare);
  union U {
    Impl::cas128_t i;
    T t;
    KOKKOS_INLINE_FUNCTION U(){};
  } newval;
  newval.t = val;
  _InterlockedCompareExchange128((LONGLONG*)dest, newval.i.upper,
                                 newval.i.lower,
                                 ((LONGLONG*)&compare_and_result));
  return compare_and_result;
}
#endif

}  // namespace Kokkos
#endif
#endif
