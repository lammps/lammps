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
#ifdef __GNUC__
__attribute__((aligned(16)))
#endif
;
}  // namespace Impl

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    typename Kokkos::Impl::enable_if<sizeof(T) == sizeof(LONG), const T&>::type
        val) {
  union U {
    LONG i;
    T t;
    KOKKOS_INLINE_FUNCTION U(){};
  } tmp;

  tmp.i = _InterlockedCompareExchange((LONG*)dest, *((LONG*)&val),
                                      *((LONG*)&compare));
  return tmp.t;
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    typename Kokkos::Impl::enable_if<sizeof(T) == sizeof(LONGLONG),
                                     const T&>::type val) {
  union U {
    LONGLONG i;
    T t;
    KOKKOS_INLINE_FUNCTION U(){};
  } tmp;

  tmp.i = _InterlockedCompareExchange64((LONGLONG*)dest, *((LONGLONG*)&val),
                                        *((LONGLONG*)&compare));
  return tmp.t;
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_compare_exchange(
    volatile T* const dest, const T& compare,
    typename Kokkos::Impl::enable_if<sizeof(T) == sizeof(Impl::cas128_t),
                                     const T&>::type val) {
  union U {
    Impl::cas128_t i;
    T t;
    KOKKOS_INLINE_FUNCTION U(){};
  } tmp, newval;
  newval.t = val;
  _InterlockedCompareExchange128((LONGLONG*)dest, newval.i.upper,
                                 newval.i.lower, ((LONGLONG*)&compare));
  tmp.t = dest;
  return tmp.t;
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_compare_exchange_strong(volatile T* const dest,
                                                        const T& compare,
                                                        const T& val) {
  return atomic_compare_exchange(dest, compare, val);
}

template <typename T>
T atomic_fetch_or(volatile T* const dest, const T val) {
  T oldval = *dest;
  T assume;
  do {
    assume   = oldval;
    T newval = val | oldval;
    oldval   = atomic_compare_exchange(dest, assume, newval);
  } while (assume != oldval);

  return oldval;
}

template <typename T>
T atomic_fetch_and(volatile T* const dest, const T val) {
  T oldval = *dest;
  T assume;
  do {
    assume   = oldval;
    T newval = val & oldval;
    oldval   = atomic_compare_exchange(dest, assume, newval);
  } while (assume != oldval);

  return oldval;
}

template <typename T>
T atomic_fetch_add(volatile T* const dest, const T val) {
  T oldval = *dest;
  T assume;
  do {
    assume   = oldval;
    T newval = val + oldval;
    oldval   = atomic_compare_exchange(dest, assume, newval);
  } while (assume != oldval);

  return oldval;
}

template <typename T>
T atomic_fetch_sub(volatile T* const dest, const T val) {
  T oldval = *dest;
  T assume;
  do {
    assume   = oldval;
    T newval = val - oldval;
    oldval   = atomic_compare_exchange(dest, assume, newval);
  } while (assume != oldval);

  return oldval;
}

template <typename T>
T atomic_exchange(volatile T* const dest, const T val) {
  T oldval = *dest;
  T assume;
  do {
    assume = oldval;
    oldval = atomic_compare_exchange(dest, assume, val);
  } while (assume != oldval);

  return oldval;
}

template <typename T>
void atomic_or(volatile T* const dest, const T val) {
  atomic_fetch_or(dest, val);
}

template <typename T>
void atomic_and(volatile T* const dest, const T val) {
  atomic_fetch_and(dest, val);
}

template <typename T>
void atomic_add(volatile T* const dest, const T val) {
  atomic_fetch_add(dest, val);
}

template <typename T>
void atomic_sub(volatile T* const dest, const T val) {
  atomic_fetch_sub(dest, val);
}

template <typename T>
void atomic_assign(volatile T* const dest, const T val) {
  atomic_fetch_exchange(dest, val);
}

template <typename T>
T atomic_increment(volatile T* const dest) {
  T oldval = *dest;
  T assume;
  do {
    assume   = oldval;
    T newval = assume++;
    oldval   = atomic_compare_exchange(dest, assume, newval);
  } while (assume != oldval);
}

template <typename T>
T atomic_decrement(volatile T* const dest) {
  T oldval = *dest;
  T assume;
  do {
    assume   = oldval;
    T newval = assume--;
    oldval   = atomic_compare_exchange(dest, assume, newval);
  } while (assume != oldval);
}

}  // namespace Kokkos
#endif
#endif
