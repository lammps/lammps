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

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_VOLATILE_LOAD_HPP)
#define KOKKOS_VOLATILE_LOAD_HPP

#if defined(__GNUC__) /* GNU C   */ || defined(__GNUG__) /* GNU C++ */ || \
    defined(__clang__)

#define KOKKOS_IMPL_MAY_ALIAS __attribute__((__may_alias__))

#else

#define KOKKOS_IMPL_MAY_ALIAS

#endif

namespace Kokkos {

//----------------------------------------------------------------------------

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION T volatile_load(T const volatile* const src_ptr) {
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
    T const volatile* const ptr;
    T64 const volatile* const ptr64;
    T32 const volatile* const ptr32;
    T16 const volatile* const ptr16;
    T8 const volatile* const ptr8;
  } src = {src_ptr};

  T result;

  union {
    T* const ptr;
    T64* const ptr64;
    T32* const ptr32;
    T16* const ptr16;
    T8* const ptr8;
  } dst = {&result};

  for (int i = 0; i < NUM_64; ++i) {
    dst.ptr64[i] = src.ptr64[i];
  }

  if (NUM_64 * 2 < NUM_32) {
    dst.ptr32[NUM_64 * 2] = src.ptr32[NUM_64 * 2];
  }

  if (NUM_32 * 2 < NUM_16) {
    dst.ptr16[NUM_32 * 2] = src.ptr16[NUM_32 * 2];
  }

  if (NUM_16 * 2 < NUM_8) {
    dst.ptr8[NUM_16 * 2] = src.ptr8[NUM_16 * 2];
  }

  return result;
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION void volatile_store(
    T volatile* const dst_ptr, T const volatile* const src_ptr) {
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
    T const volatile* const ptr;
    T64 const volatile* const ptr64;
    T32 const volatile* const ptr32;
    T16 const volatile* const ptr16;
    T8 const volatile* const ptr8;
  } src = {src_ptr};

  union {
    T volatile* const ptr;
    T64 volatile* const ptr64;
    T32 volatile* const ptr32;
    T16 volatile* const ptr16;
    T8 volatile* const ptr8;
  } dst = {dst_ptr};

  for (int i = 0; i < NUM_64; ++i) {
    dst.ptr64[i] = src.ptr64[i];
  }

  if (NUM_64 * 2 < NUM_32) {
    dst.ptr32[NUM_64 * 2] = src.ptr32[NUM_64 * 2];
  }

  if (NUM_32 * 2 < NUM_16) {
    dst.ptr16[NUM_32 * 2] = src.ptr16[NUM_32 * 2];
  }

  if (NUM_16 * 2 < NUM_8) {
    dst.ptr8[NUM_16 * 2] = src.ptr8[NUM_16 * 2];
  }
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION void volatile_store(T volatile* const dst_ptr,
                                                T const* const src_ptr) {
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
  } src = {src_ptr};

  union {
    T volatile* const ptr;
    T64 volatile* const ptr64;
    T32 volatile* const ptr32;
    T16 volatile* const ptr16;
    T8 volatile* const ptr8;
  } dst = {dst_ptr};

  for (int i = 0; i < NUM_64; ++i) {
    dst.ptr64[i] = src.ptr64[i];
  }

  if (NUM_64 * 2 < NUM_32) {
    dst.ptr32[NUM_64 * 2] = src.ptr32[NUM_64 * 2];
  }

  if (NUM_32 * 2 < NUM_16) {
    dst.ptr16[NUM_32 * 2] = src.ptr16[NUM_32 * 2];
  }

  if (NUM_16 * 2 < NUM_8) {
    dst.ptr8[NUM_16 * 2] = src.ptr8[NUM_16 * 2];
  }
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION void volatile_store(T volatile* dst_ptr,
                                                T const volatile& src) {
  volatile_store(dst_ptr, &src);
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION void volatile_store(T volatile* dst_ptr,
                                                T const& src) {
  volatile_store(dst_ptr, &src);
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION T safe_load(T const* const ptr) {
#if !defined(__MIC__)
  return *ptr;
#else
  return volatile_load(ptr);
#endif
}

}  // namespace Kokkos

#undef KOKKOS_IMPL_MAY_ALIAS

#endif
