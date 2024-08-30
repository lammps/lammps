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

#ifndef KOKKOS_FUNCTIONAL_IMPL_HPP
#define KOKKOS_FUNCTIONAL_IMPL_HPP

#include <Kokkos_Macros.hpp>
#include <cstdint>

namespace Kokkos {
namespace Impl {

// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.
KOKKOS_FORCEINLINE_FUNCTION
uint32_t getblock32(const uint8_t* p, int i) {
  // used to avoid aliasing error which could cause errors with
  // forced inlining
  return ((uint32_t)p[i * 4 + 0]) | ((uint32_t)p[i * 4 + 1] << 8) |
         ((uint32_t)p[i * 4 + 2] << 16) | ((uint32_t)p[i * 4 + 3] << 24);
}

KOKKOS_FORCEINLINE_FUNCTION
uint32_t rotl32(uint32_t x, int8_t r) { return (x << r) | (x >> (32 - r)); }

KOKKOS_FORCEINLINE_FUNCTION
uint32_t fmix32(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

KOKKOS_INLINE_FUNCTION
uint32_t MurmurHash3_x86_32(const void* key, int len, uint32_t seed) {
  const uint8_t* data = static_cast<const uint8_t*>(key);
  const int nblocks   = len / 4;

  uint32_t h1 = seed;

  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;

  //----------
  // body

  for (int i = 0; i < nblocks; ++i) {
    uint32_t k1 = getblock32(data, i);

    k1 *= c1;
    k1 = rotl32(k1, 15);
    k1 *= c2;

    h1 ^= k1;
    h1 = rotl32(h1, 13);
    h1 = h1 * 5 + 0xe6546b64;
  }

  //----------
  // tail

  const uint8_t* tail = (const uint8_t*)(data + nblocks * 4);

  uint32_t k1 = 0;

  switch (len & 3) {
    case 3: k1 ^= tail[2] << 16; [[fallthrough]];
    case 2: k1 ^= tail[1] << 8; [[fallthrough]];
    case 1:
      k1 ^= tail[0];
      k1 *= c1;
      k1 = rotl32(k1, 15);
      k1 *= c2;
      h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len;

  h1 = fmix32(h1);

  return h1;
}

#if defined(__GNUC__) /* GNU C   */ || defined(__GNUG__) /* GNU C++ */ || \
    defined(__clang__)

#define KOKKOS_IMPL_MAY_ALIAS __attribute__((__may_alias__))

#else

#define KOKKOS_IMPL_MAY_ALIAS

#endif

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION bool bitwise_equal(T const* const a_ptr,
                                               T const* const b_ptr) {
  typedef uint64_t KOKKOS_IMPL_MAY_ALIAS T64;  // NOLINT(modernize-use-using)
  typedef uint32_t KOKKOS_IMPL_MAY_ALIAS T32;  // NOLINT(modernize-use-using)
  typedef uint16_t KOKKOS_IMPL_MAY_ALIAS T16;  // NOLINT(modernize-use-using)
  typedef uint8_t KOKKOS_IMPL_MAY_ALIAS T8;    // NOLINT(modernize-use-using)

  enum {
    NUM_8  = sizeof(T),
    NUM_16 = NUM_8 / 2,
    NUM_32 = NUM_8 / 4,
    NUM_64 = NUM_8 / 8
  };

  union {
    T const* const ptr;
    T64 const* const ptr64;
    T32 const* const ptr32;
    T16 const* const ptr16;
    T8 const* const ptr8;
  } a = {a_ptr}, b = {b_ptr};

  bool result = true;

  for (int i = 0; i < NUM_64; ++i) {
    result = result && a.ptr64[i] == b.ptr64[i];
  }

  if (NUM_64 * 2 < NUM_32) {
    result = result && a.ptr32[NUM_64 * 2] == b.ptr32[NUM_64 * 2];
  }

  if (NUM_32 * 2 < NUM_16) {
    result = result && a.ptr16[NUM_32 * 2] == b.ptr16[NUM_32 * 2];
  }

  if (NUM_16 * 2 < NUM_8) {
    result = result && a.ptr8[NUM_16 * 2] == b.ptr8[NUM_16 * 2];
  }

  return result;
}

#undef KOKKOS_IMPL_MAY_ALIAS

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_FUNCTIONAL_IMPL_HPP
