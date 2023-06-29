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

#ifndef KOKKOS_BITOPS_HPP
#define KOKKOS_BITOPS_HPP

#include <Kokkos_Macros.hpp>
#include <cstdint>
#include <climits>

#if defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_INTEL_LLVM)
#include <immintrin.h>
#endif

namespace Kokkos {
namespace Impl {

KOKKOS_FORCEINLINE_FUNCTION
int int_log2_fallback(unsigned i) {
  constexpr int shift = sizeof(unsigned) * CHAR_BIT - 1;

  int offset = 0;
  if (i) {
    for (offset = shift; (i & (1 << offset)) == 0; --offset)
      ;
  }
  return offset;
}

KOKKOS_IMPL_DEVICE_FUNCTION
inline int int_log2_device(unsigned i) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  constexpr int shift = sizeof(unsigned) * CHAR_BIT - 1;
  return shift - __clz(i);
#elif defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_INTEL_LLVM)
  return _bit_scan_reverse(i);
#else
  return int_log2_fallback(i);
#endif
}

KOKKOS_IMPL_HOST_FUNCTION
inline int int_log2_host(unsigned i) {
// duplicating shift to avoid unused warning in else branch
#if defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_INTEL_LLVM)
  constexpr int shift = sizeof(unsigned) * CHAR_BIT - 1;
  (void)shift;
  return _bit_scan_reverse(i);
#elif defined(KOKKOS_COMPILER_CRAYC)
  constexpr int shift = sizeof(unsigned) * CHAR_BIT - 1;
  return i ? shift - _leadz32(i) : 0;
#elif defined(__GNUC__) || defined(__GNUG__)
  constexpr int shift = sizeof(unsigned) * CHAR_BIT - 1;
  return shift - __builtin_clz(i);
#else
  return int_log2_fallback(i);
#endif
}

#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma push
#pragma diag_suppress implicit_return_from_non_void_function
#endif
KOKKOS_FORCEINLINE_FUNCTION
int int_log2(unsigned i) {
  KOKKOS_IF_ON_DEVICE((return int_log2_device(i);))
  KOKKOS_IF_ON_HOST((return int_log2_host(i);))
}
#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma pop
#endif

/**\brief  Find first zero bit.
 *
 *  If none then return -1 ;
 */
KOKKOS_FORCEINLINE_FUNCTION
int bit_first_zero_fallback(unsigned i) noexcept {
  constexpr unsigned full = ~0u;

  int offset = -1;
  if (full != i) {
    for (offset = 0; i & (1 << offset); ++offset)
      ;
  }
  return offset;
}

KOKKOS_IMPL_DEVICE_FUNCTION
inline int bit_first_zero_device(unsigned i) noexcept {
  constexpr unsigned full = ~0u;
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  return full != i ? __ffs(~i) - 1 : -1;
#elif defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_INTEL_LLVM)
  return full != i ? _bit_scan_forward(~i) : -1;
#else
  (void)full;
  return bit_first_zero_fallback(i);
#endif
}

KOKKOS_IMPL_HOST_FUNCTION
inline int bit_first_zero_host(unsigned i) noexcept {
  constexpr unsigned full = ~0u;
#if defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_INTEL_LLVM)
  return full != i ? _bit_scan_forward(~i) : -1;
#elif defined(KOKKOS_COMPILER_CRAYC)
  return full != i ? _popcnt(i ^ (i + 1)) - 1 : -1;
#elif defined(KOKKOS_COMPILER_GNU) || defined(__GNUC__) || defined(__GNUG__)
  return full != i ? __builtin_ffs(~i) - 1 : -1;
#else
  (void)full;
  return bit_first_zero_fallback(i);
#endif
}

#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma push
#pragma diag_suppress implicit_return_from_non_void_function
#endif
KOKKOS_FORCEINLINE_FUNCTION
int bit_first_zero(unsigned i) noexcept {
  KOKKOS_IF_ON_DEVICE((return bit_first_zero_device(i);))
  KOKKOS_IF_ON_HOST((return bit_first_zero_host(i);))
}
#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma pop
#endif

KOKKOS_FORCEINLINE_FUNCTION
int bit_scan_forward_fallback(unsigned i) {
  int offset = -1;
  if (i) {
    for (offset = 0; (i & (1 << offset)) == 0; ++offset)
      ;
  }
  return offset;
}

KOKKOS_IMPL_DEVICE_FUNCTION inline int bit_scan_forward_device(unsigned i) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  return __ffs(i) - 1;
#elif defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_INTEL_LLVM)
  return _bit_scan_forward(i);
#else
  return bit_scan_forward_fallback(i);
#endif
}

KOKKOS_IMPL_HOST_FUNCTION inline int bit_scan_forward_host(unsigned i) {
#if defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_INTEL_LLVM)
  return _bit_scan_forward(i);
#elif defined(KOKKOS_COMPILER_CRAYC)
  return i ? _popcnt(~i & (i - 1)) : -1;
#elif defined(KOKKOS_COMPILER_GNU) || defined(__GNUC__) || defined(__GNUG__)
  return __builtin_ffs(i) - 1;
#else
  return bit_scan_forward_fallback(i);
#endif
}

#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma push
#pragma diag_suppress implicit_return_from_non_void_function
#endif
KOKKOS_FORCEINLINE_FUNCTION
int bit_scan_forward(unsigned i) {
  KOKKOS_IF_ON_DEVICE((return bit_scan_forward_device(i);))
  KOKKOS_IF_ON_HOST((return bit_scan_forward_host(i);))
}
#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma pop
#endif

/// Count the number of bits set.
KOKKOS_FORCEINLINE_FUNCTION
int bit_count_fallback(unsigned i) {
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive
  i = i - ((i >> 1) & ~0u / 3u);                           // temp
  i = (i & ~0u / 15u * 3u) + ((i >> 2) & ~0u / 15u * 3u);  // temp
  i = (i + (i >> 4)) & ~0u / 255u * 15u;                   // temp

  // count
  return (int)((i * (~0u / 255u)) >> (sizeof(unsigned) - 1) * CHAR_BIT);
}

KOKKOS_IMPL_DEVICE_FUNCTION inline int bit_count_device(unsigned i) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  return __popc(i);
#elif defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_INTEL_LLVM)
  return _popcnt32(i);
#else
  return bit_count_fallback(i);
#endif
}

KOKKOS_IMPL_HOST_FUNCTION inline int bit_count_host(unsigned i) {
#if defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_INTEL_LLVM)
  return _popcnt32(i);
#elif defined(KOKKOS_COMPILER_CRAYC)
  return _popcnt(i);
#elif defined(__GNUC__) || defined(__GNUG__)
  return __builtin_popcount(i);
#else
  return bit_count_fallback(i);
#endif
}

#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma push
#pragma diag_suppress implicit_return_from_non_void_function
#endif
KOKKOS_FORCEINLINE_FUNCTION
int bit_count(unsigned i) {
  KOKKOS_IF_ON_DEVICE((return bit_count_device(i);))
  KOKKOS_IF_ON_HOST((return bit_count_host(i);))
}
#if defined(__EDG__) && !defined(KOKKOS_COMPILER_INTEL)
#pragma pop
#endif

KOKKOS_INLINE_FUNCTION
unsigned integral_power_of_two_that_contains(const unsigned N) {
  const unsigned i = int_log2(N);
  return ((1u << i) < N) ? i + 1 : i;
}

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_BITOPS_HPP
