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

#ifndef KOKKOS_HOST_EXP_ITERATE_TILE_HPP
#define KOKKOS_HOST_EXP_ITERATE_TILE_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_AGGRESSIVE_VECTORIZATION) && \
    defined(KOKKOS_ENABLE_PRAGMA_IVDEP) && !defined(__CUDA_ARCH__)
#define KOKKOS_MDRANGE_IVDEP
#endif

#ifdef KOKKOS_MDRANGE_IVDEP
#define KOKKOS_ENABLE_IVDEP_MDRANGE _Pragma("ivdep")
#else
#define KOKKOS_ENABLE_IVDEP_MDRANGE
#endif

#include <algorithm>

namespace Kokkos {
namespace Impl {

// Temporary, for testing new loop macros
#define KOKKOS_ENABLE_NEW_LOOP_MACROS 1

#define KOKKOS_IMPL_LOOP_1L(type, tile) \
  KOKKOS_ENABLE_IVDEP_MDRANGE           \
  for (type i0 = 0; i0 < static_cast<type>(tile[0]); ++i0)

#define KOKKOS_IMPL_LOOP_2L(type, tile)                    \
  for (type i1 = 0; i1 < static_cast<type>(tile[1]); ++i1) \
  KOKKOS_IMPL_LOOP_1L(type, tile)

#define KOKKOS_IMPL_LOOP_3L(type, tile)                    \
  for (type i2 = 0; i2 < static_cast<type>(tile[2]); ++i2) \
  KOKKOS_IMPL_LOOP_2L(type, tile)

#define KOKKOS_IMPL_LOOP_4L(type, tile)                    \
  for (type i3 = 0; i3 < static_cast<type>(tile[3]); ++i3) \
  KOKKOS_IMPL_LOOP_3L(type, tile)

#define KOKKOS_IMPL_LOOP_5L(type, tile)                    \
  for (type i4 = 0; i4 < static_cast<type>(tile[4]); ++i4) \
  KOKKOS_IMPL_LOOP_4L(type, tile)

#define KOKKOS_IMPL_LOOP_6L(type, tile)                    \
  for (type i5 = 0; i5 < static_cast<type>(tile[5]); ++i5) \
  KOKKOS_IMPL_LOOP_5L(type, tile)

#define KOKKOS_IMPL_LOOP_7L(type, tile)                    \
  for (type i6 = 0; i6 < static_cast<type>(tile[6]); ++i6) \
  KOKKOS_IMPL_LOOP_6L(type, tile)

#define KOKKOS_IMPL_LOOP_8L(type, tile)                    \
  for (type i7 = 0; i7 < static_cast<type>(tile[7]); ++i7) \
  KOKKOS_IMPL_LOOP_7L(type, tile)

#define KOKKOS_IMPL_LOOP_1R(type, tile) \
  KOKKOS_ENABLE_IVDEP_MDRANGE           \
  for (type i0 = 0; i0 < static_cast<type>(tile[0]); ++i0)

#define KOKKOS_IMPL_LOOP_2R(type, tile) \
  KOKKOS_IMPL_LOOP_1R(type, tile)       \
  for (type i1 = 0; i1 < static_cast<type>(tile[1]); ++i1)

#define KOKKOS_IMPL_LOOP_3R(type, tile) \
  KOKKOS_IMPL_LOOP_2R(type, tile)       \
  for (type i2 = 0; i2 < static_cast<type>(tile[2]); ++i2)

#define KOKKOS_IMPL_LOOP_4R(type, tile) \
  KOKKOS_IMPL_LOOP_3R(type, tile)       \
  for (type i3 = 0; i3 < static_cast<type>(tile[3]); ++i3)

#define KOKKOS_IMPL_LOOP_5R(type, tile) \
  KOKKOS_IMPL_LOOP_4R(type, tile)       \
  for (type i4 = 0; i4 < static_cast<type>(tile[4]); ++i4)

#define KOKKOS_IMPL_LOOP_6R(type, tile) \
  KOKKOS_IMPL_LOOP_5R(type, tile)       \
  for (type i5 = 0; i5 < static_cast<type>(tile[5]); ++i5)

#define KOKKOS_IMPL_LOOP_7R(type, tile) \
  KOKKOS_IMPL_LOOP_6R(type, tile)       \
  for (type i6 = 0; i6 < static_cast<type>(tile[6]); ++i6)

#define KOKKOS_IMPL_LOOP_8R(type, tile) \
  KOKKOS_IMPL_LOOP_7R(type, tile)       \
  for (type i7 = 0; i7 < static_cast<type>(tile[7]); ++i7)

#define KOKKOS_IMPL_LOOP_ARGS_1 i0 + m_offset[0]
#define KOKKOS_IMPL_LOOP_ARGS_2 KOKKOS_IMPL_LOOP_ARGS_1, i1 + m_offset[1]
#define KOKKOS_IMPL_LOOP_ARGS_3 KOKKOS_IMPL_LOOP_ARGS_2, i2 + m_offset[2]
#define KOKKOS_IMPL_LOOP_ARGS_4 KOKKOS_IMPL_LOOP_ARGS_3, i3 + m_offset[3]
#define KOKKOS_IMPL_LOOP_ARGS_5 KOKKOS_IMPL_LOOP_ARGS_4, i4 + m_offset[4]
#define KOKKOS_IMPL_LOOP_ARGS_6 KOKKOS_IMPL_LOOP_ARGS_5, i5 + m_offset[5]
#define KOKKOS_IMPL_LOOP_ARGS_7 KOKKOS_IMPL_LOOP_ARGS_6, i6 + m_offset[6]
#define KOKKOS_IMPL_LOOP_ARGS_8 KOKKOS_IMPL_LOOP_ARGS_7, i7 + m_offset[7]

// New Loop Macros...
// parallel_for, non-tagged
#define KOKKOS_IMPL_APPLY(func, ...) func(__VA_ARGS__);

// LayoutRight
// d = 0 to start
#define KOKKOS_IMPL_LOOP_R_1(func, type, m_offset, extent, d, ...)   \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                        \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) { \
    KOKKOS_IMPL_APPLY(func, __VA_ARGS__, i0 + m_offset[d])           \
  }

#define KOKKOS_IMPL_LOOP_R_2(func, type, m_offset, extent, d, ...)         \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {       \
    KOKKOS_IMPL_LOOP_R_1(func, type, m_offset, extent, d + 1, __VA_ARGS__, \
                         i1 + m_offset[d])                                 \
  }

#define KOKKOS_IMPL_LOOP_R_3(func, type, m_offset, extent, d, ...)         \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {       \
    KOKKOS_IMPL_LOOP_R_2(func, type, m_offset, extent, d + 1, __VA_ARGS__, \
                         i2 + m_offset[d])                                 \
  }

#define KOKKOS_IMPL_LOOP_R_4(func, type, m_offset, extent, d, ...)         \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {       \
    KOKKOS_IMPL_LOOP_R_3(func, type, m_offset, extent, d + 1, __VA_ARGS__, \
                         i3 + m_offset[d])                                 \
  }

#define KOKKOS_IMPL_LOOP_R_5(func, type, m_offset, extent, d, ...)         \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {       \
    KOKKOS_IMPL_LOOP_R_4(func, type, m_offset, extent, d + 1, __VA_ARGS__, \
                         i4 + m_offset[d])                                 \
  }

#define KOKKOS_IMPL_LOOP_R_6(func, type, m_offset, extent, d, ...)         \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {       \
    KOKKOS_IMPL_LOOP_R_5(func, type, m_offset, extent, d + 1, __VA_ARGS__, \
                         i5 + m_offset[d])                                 \
  }

#define KOKKOS_IMPL_LOOP_R_7(func, type, m_offset, extent, d, ...)         \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {       \
    KOKKOS_IMPL_LOOP_R_6(func, type, m_offset, extent, d + 1, __VA_ARGS__, \
                         i6 + m_offset[d])                                 \
  }

#define KOKKOS_IMPL_LOOP_R_8(func, type, m_offset, extent, d, ...)         \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {       \
    KOKKOS_IMPL_LOOP_R_7(func, type, m_offset, extent, d + 1, __VA_ARGS__, \
                         i7 + m_offset[d])                                 \
  }

// LayoutLeft
// d = rank-1 to start
#define KOKKOS_IMPL_LOOP_L_1(func, type, m_offset, extent, d, ...)   \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                        \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) { \
    KOKKOS_IMPL_APPLY(func, i0 + m_offset[d], __VA_ARGS__)           \
  }

#define KOKKOS_IMPL_LOOP_L_2(func, type, m_offset, extent, d, ...)   \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) { \
    KOKKOS_IMPL_LOOP_L_1(func, type, m_offset, extent, d - 1,        \
                         i1 + m_offset[d], __VA_ARGS__)              \
  }

#define KOKKOS_IMPL_LOOP_L_3(func, type, m_offset, extent, d, ...)   \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) { \
    KOKKOS_IMPL_LOOP_L_2(func, type, m_offset, extent, d - 1,        \
                         i2 + m_offset[d], __VA_ARGS__)              \
  }

#define KOKKOS_IMPL_LOOP_L_4(func, type, m_offset, extent, d, ...)   \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) { \
    KOKKOS_IMPL_LOOP_L_3(func, type, m_offset, extent, d - 1,        \
                         i3 + m_offset[d], __VA_ARGS__)              \
  }

#define KOKKOS_IMPL_LOOP_L_5(func, type, m_offset, extent, d, ...)   \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) { \
    KOKKOS_IMPL_LOOP_L_4(func, type, m_offset, extent, d - 1,        \
                         i4 + m_offset[d], __VA_ARGS__)              \
  }

#define KOKKOS_IMPL_LOOP_L_6(func, type, m_offset, extent, d, ...)   \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) { \
    KOKKOS_IMPL_LOOP_L_5(func, type, m_offset, extent, d - 1,        \
                         i5 + m_offset[d], __VA_ARGS__)              \
  }

#define KOKKOS_IMPL_LOOP_L_7(func, type, m_offset, extent, d, ...)   \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) { \
    KOKKOS_IMPL_LOOP_L_6(func, type, m_offset, extent, d - 1,        \
                         i6 + m_offset[d], __VA_ARGS__)              \
  }

#define KOKKOS_IMPL_LOOP_L_8(func, type, m_offset, extent, d, ...)   \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) { \
    KOKKOS_IMPL_LOOP_L_7(func, type, m_offset, extent, d - 1,        \
                         i7 + m_offset[d], __VA_ARGS__)              \
  }

// Left vs Right
// TODO: rank not necessary to pass through, can hardcode the values
#define KOKKOS_IMPL_LOOP_LAYOUT_1(func, type, is_left, m_offset, extent, rank) \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                  \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[0]); ++i0) {           \
    KOKKOS_IMPL_APPLY(func, i0 + m_offset[0])                                  \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_2(func, type, is_left, m_offset, extent, rank) \
  if (is_left) {                                                               \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[rank - 1]); ++i1) {  \
      KOKKOS_IMPL_LOOP_L_1(func, type, m_offset, extent, rank - 2,             \
                           i1 + m_offset[rank - 1])                            \
    }                                                                          \
  } else {                                                                     \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[0]); ++i1) {         \
      KOKKOS_IMPL_LOOP_R_1(func, type, m_offset, extent, 1, i1 + m_offset[0])  \
    }                                                                          \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_3(func, type, is_left, m_offset, extent, rank) \
  if (is_left) {                                                               \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[rank - 1]); ++i2) {  \
      KOKKOS_IMPL_LOOP_L_2(func, type, m_offset, extent, rank - 2,             \
                           i2 + m_offset[rank - 1])                            \
    }                                                                          \
  } else {                                                                     \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[0]); ++i2) {         \
      KOKKOS_IMPL_LOOP_R_2(func, type, m_offset, extent, 1, i2 + m_offset[0])  \
    }                                                                          \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_4(func, type, is_left, m_offset, extent, rank) \
  if (is_left) {                                                               \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[rank - 1]); ++i3) {  \
      KOKKOS_IMPL_LOOP_L_3(func, type, m_offset, extent, rank - 2,             \
                           i3 + m_offset[rank - 1])                            \
    }                                                                          \
  } else {                                                                     \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[0]); ++i3) {         \
      KOKKOS_IMPL_LOOP_R_3(func, type, m_offset, extent, 1, i3 + m_offset[0])  \
    }                                                                          \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_5(func, type, is_left, m_offset, extent, rank) \
  if (is_left) {                                                               \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[rank - 1]); ++i4) {  \
      KOKKOS_IMPL_LOOP_L_4(func, type, m_offset, extent, rank - 2,             \
                           i4 + m_offset[rank - 1])                            \
    }                                                                          \
  } else {                                                                     \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[0]); ++i4) {         \
      KOKKOS_IMPL_LOOP_R_4(func, type, m_offset, extent, 1, i4 + m_offset[0])  \
    }                                                                          \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_6(func, type, is_left, m_offset, extent, rank) \
  if (is_left) {                                                               \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[rank - 1]); ++i5) {  \
      KOKKOS_IMPL_LOOP_L_5(func, type, m_offset, extent, rank - 2,             \
                           i5 + m_offset[rank - 1])                            \
    }                                                                          \
  } else {                                                                     \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[0]); ++i5) {         \
      KOKKOS_IMPL_LOOP_R_5(func, type, m_offset, extent, 1, i5 + m_offset[0])  \
    }                                                                          \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_7(func, type, is_left, m_offset, extent, rank) \
  if (is_left) {                                                               \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[rank - 1]); ++i6) {  \
      KOKKOS_IMPL_LOOP_L_6(func, type, m_offset, extent, rank - 2,             \
                           i6 + m_offset[rank - 1])                            \
    }                                                                          \
  } else {                                                                     \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[0]); ++i6) {         \
      KOKKOS_IMPL_LOOP_R_6(func, type, m_offset, extent, 1, i6 + m_offset[0])  \
    }                                                                          \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_8(func, type, is_left, m_offset, extent, rank) \
  if (is_left) {                                                               \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[rank - 1]); ++i7) {  \
      KOKKOS_IMPL_LOOP_L_7(func, type, m_offset, extent, rank - 2,             \
                           i7 + m_offset[rank - 1])                            \
    }                                                                          \
  } else {                                                                     \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[0]); ++i7) {         \
      KOKKOS_IMPL_LOOP_R_7(func, type, m_offset, extent, 1, i7 + m_offset[0])  \
    }                                                                          \
  }

// Partial vs Full Tile
#define KOKKOS_IMPL_TILE_LOOP_1(func, type, is_left, cond, m_offset,         \
                                extent_full, extent_partial, rank)           \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_1(func, type, is_left, m_offset, extent_full,    \
                              rank)                                          \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_1(func, type, is_left, m_offset, extent_partial, \
                              rank)                                          \
  }

#define KOKKOS_IMPL_TILE_LOOP_2(func, type, is_left, cond, m_offset,         \
                                extent_full, extent_partial, rank)           \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_2(func, type, is_left, m_offset, extent_full,    \
                              rank)                                          \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_2(func, type, is_left, m_offset, extent_partial, \
                              rank)                                          \
  }

#define KOKKOS_IMPL_TILE_LOOP_3(func, type, is_left, cond, m_offset,         \
                                extent_full, extent_partial, rank)           \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_3(func, type, is_left, m_offset, extent_full,    \
                              rank)                                          \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_3(func, type, is_left, m_offset, extent_partial, \
                              rank)                                          \
  }

#define KOKKOS_IMPL_TILE_LOOP_4(func, type, is_left, cond, m_offset,         \
                                extent_full, extent_partial, rank)           \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_4(func, type, is_left, m_offset, extent_full,    \
                              rank)                                          \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_4(func, type, is_left, m_offset, extent_partial, \
                              rank)                                          \
  }

#define KOKKOS_IMPL_TILE_LOOP_5(func, type, is_left, cond, m_offset,         \
                                extent_full, extent_partial, rank)           \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_5(func, type, is_left, m_offset, extent_full,    \
                              rank)                                          \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_5(func, type, is_left, m_offset, extent_partial, \
                              rank)                                          \
  }

#define KOKKOS_IMPL_TILE_LOOP_6(func, type, is_left, cond, m_offset,         \
                                extent_full, extent_partial, rank)           \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_6(func, type, is_left, m_offset, extent_full,    \
                              rank)                                          \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_6(func, type, is_left, m_offset, extent_partial, \
                              rank)                                          \
  }

#define KOKKOS_IMPL_TILE_LOOP_7(func, type, is_left, cond, m_offset,         \
                                extent_full, extent_partial, rank)           \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_7(func, type, is_left, m_offset, extent_full,    \
                              rank)                                          \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_7(func, type, is_left, m_offset, extent_partial, \
                              rank)                                          \
  }

#define KOKKOS_IMPL_TILE_LOOP_8(func, type, is_left, cond, m_offset,         \
                                extent_full, extent_partial, rank)           \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_8(func, type, is_left, m_offset, extent_full,    \
                              rank)                                          \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_8(func, type, is_left, m_offset, extent_partial, \
                              rank)                                          \
  }

// parallel_reduce, non-tagged
// Reduction version
#define KOKKOS_IMPL_APPLY_REDUX(val, func, ...) func(__VA_ARGS__, val);

// LayoutRight
// d = 0 to start
#define KOKKOS_IMPL_LOOP_R_1_REDUX(val, func, type, m_offset, extent, d, ...) \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                 \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {          \
    KOKKOS_IMPL_APPLY_REDUX(val, func, __VA_ARGS__, i0 + m_offset[d])         \
  }

#define KOKKOS_IMPL_LOOP_R_2_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {          \
    KOKKOS_IMPL_LOOP_R_1_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i1 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_3_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {          \
    KOKKOS_IMPL_LOOP_R_2_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i2 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_4_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {          \
    KOKKOS_IMPL_LOOP_R_3_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i3 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_5_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {          \
    KOKKOS_IMPL_LOOP_R_4_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i4 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_6_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {          \
    KOKKOS_IMPL_LOOP_R_5_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i5 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_7_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {          \
    KOKKOS_IMPL_LOOP_R_6_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i6 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_LOOP_R_8_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {          \
    KOKKOS_IMPL_LOOP_R_7_REDUX(val, func, type, m_offset, extent, d + 1,      \
                               __VA_ARGS__, i7 + m_offset[d])                 \
  }

// LayoutLeft
// d = rank-1 to start
#define KOKKOS_IMPL_LOOP_L_1_REDUX(val, func, type, m_offset, extent, d, ...) \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                 \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {          \
    KOKKOS_IMPL_APPLY_REDUX(val, func, i0 + m_offset[d], __VA_ARGS__)         \
  }

#define KOKKOS_IMPL_LOOP_L_2_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {          \
    KOKKOS_IMPL_LOOP_L_1_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i1 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_3_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {          \
    KOKKOS_IMPL_LOOP_L_2_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i2 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_4_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {          \
    KOKKOS_IMPL_LOOP_L_3_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i3 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_5_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {          \
    KOKKOS_IMPL_LOOP_L_4_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i4 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_6_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {          \
    KOKKOS_IMPL_LOOP_L_5_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i5 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_7_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {          \
    KOKKOS_IMPL_LOOP_L_6_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i6 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_LOOP_L_8_REDUX(val, func, type, m_offset, extent, d, ...) \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {          \
    KOKKOS_IMPL_LOOP_L_7_REDUX(val, func, type, m_offset, extent, d - 1,      \
                               i7 + m_offset[d], __VA_ARGS__)                 \
  }

// Left vs Right
#define KOKKOS_IMPL_LOOP_LAYOUT_1_REDUX(val, func, type, is_left, m_offset, \
                                        extent, rank)                       \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                               \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[0]); ++i0) {        \
    KOKKOS_IMPL_APPLY_REDUX(val, func, i0 + m_offset[0])                    \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_2_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[rank - 1]); ++i1) { \
      KOKKOS_IMPL_LOOP_L_1_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i1 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[0]); ++i1) {        \
      KOKKOS_IMPL_LOOP_R_1_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i1 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_3_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[rank - 1]); ++i2) { \
      KOKKOS_IMPL_LOOP_L_2_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i2 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[0]); ++i2) {        \
      KOKKOS_IMPL_LOOP_R_2_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i2 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_4_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[rank - 1]); ++i3) { \
      KOKKOS_IMPL_LOOP_L_3_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i3 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[0]); ++i3) {        \
      KOKKOS_IMPL_LOOP_R_3_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i3 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_5_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[rank - 1]); ++i4) { \
      KOKKOS_IMPL_LOOP_L_4_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i4 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[0]); ++i4) {        \
      KOKKOS_IMPL_LOOP_R_4_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i4 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_6_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[rank - 1]); ++i5) { \
      KOKKOS_IMPL_LOOP_L_5_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i5 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[0]); ++i5) {        \
      KOKKOS_IMPL_LOOP_R_5_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i5 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_7_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[rank - 1]); ++i6) { \
      KOKKOS_IMPL_LOOP_L_6_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i6 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[0]); ++i6) {        \
      KOKKOS_IMPL_LOOP_R_6_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i6 + m_offset[0])                            \
    }                                                                         \
  }

#define KOKKOS_IMPL_LOOP_LAYOUT_8_REDUX(val, func, type, is_left, m_offset,   \
                                        extent, rank)                         \
  if (is_left) {                                                              \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[rank - 1]); ++i7) { \
      KOKKOS_IMPL_LOOP_L_7_REDUX(val, func, type, m_offset, extent, rank - 2, \
                                 i7 + m_offset[rank - 1])                     \
    }                                                                         \
  } else {                                                                    \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[0]); ++i7) {        \
      KOKKOS_IMPL_LOOP_R_7_REDUX(val, func, type, m_offset, extent, 1,        \
                                 i7 + m_offset[0])                            \
    }                                                                         \
  }

// Partial vs Full Tile
#define KOKKOS_IMPL_TILE_LOOP_1_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_1_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_1_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_2_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_2_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_2_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_3_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_3_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_3_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_4_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_4_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_4_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_5_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_5_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_5_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_6_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_6_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_6_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_7_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_7_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_7_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TILE_LOOP_8_REDUX(val, func, type, is_left, cond,        \
                                      m_offset, extent_full, extent_partial, \
                                      rank)                                  \
  if (cond) {                                                                \
    KOKKOS_IMPL_LOOP_LAYOUT_8_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_full, rank)                       \
  } else {                                                                   \
    KOKKOS_IMPL_LOOP_LAYOUT_8_REDUX(val, func, type, is_left, m_offset,      \
                                    extent_partial, rank)                    \
  }
// end New Loop Macros

// tagged macros
#define KOKKOS_IMPL_TAGGED_APPLY(tag, func, ...) func(tag, __VA_ARGS__);

// LayoutRight
// d = 0 to start
#define KOKKOS_IMPL_TAGGED_LOOP_R_1(tag, func, type, m_offset, extent, d, ...) \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                  \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {           \
    KOKKOS_IMPL_TAGGED_APPLY(tag, func, __VA_ARGS__, i0 + m_offset[d])         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_2(tag, func, type, m_offset, extent, d, ...) \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {           \
    KOKKOS_IMPL_TAGGED_LOOP_R_1(tag, func, type, m_offset, extent, d + 1,      \
                                __VA_ARGS__, i1 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_3(tag, func, type, m_offset, extent, d, ...) \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {           \
    KOKKOS_IMPL_TAGGED_LOOP_R_2(tag, func, type, m_offset, extent, d + 1,      \
                                __VA_ARGS__, i2 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_4(tag, func, type, m_offset, extent, d, ...) \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {           \
    KOKKOS_IMPL_TAGGED_LOOP_R_3(tag, func, type, m_offset, extent, d + 1,      \
                                __VA_ARGS__, i3 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_5(tag, func, type, m_offset, extent, d, ...) \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {           \
    KOKKOS_IMPL_TAGGED_LOOP_R_4(tag, func, type, m_offset, extent, d + 1,      \
                                __VA_ARGS__, i4 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_6(tag, func, type, m_offset, extent, d, ...) \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {           \
    KOKKOS_IMPL_TAGGED_LOOP_R_5(tag, func, type, m_offset, extent, d + 1,      \
                                __VA_ARGS__, i5 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_7(tag, func, type, m_offset, extent, d, ...) \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {           \
    KOKKOS_IMPL_TAGGED_LOOP_R_6(tag, func, type, m_offset, extent, d + 1,      \
                                __VA_ARGS__, i6 + m_offset[d])                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_8(tag, func, type, m_offset, extent, d, ...) \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {           \
    KOKKOS_IMPL_TAGGED_LOOP_R_7(tag, func, type, m_offset, extent, d + 1,      \
                                __VA_ARGS__, i7 + m_offset[d])                 \
  }

// LayoutLeft
// d = rank-1 to start
#define KOKKOS_IMPL_TAGGED_LOOP_L_1(tag, func, type, m_offset, extent, d, ...) \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                  \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {           \
    KOKKOS_IMPL_TAGGED_APPLY(tag, func, i0 + m_offset[d], __VA_ARGS__)         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_2(tag, func, type, m_offset, extent, d, ...) \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {           \
    KOKKOS_IMPL_TAGGED_LOOP_L_1(tag, func, type, m_offset, extent, d - 1,      \
                                i1 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_3(tag, func, type, m_offset, extent, d, ...) \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {           \
    KOKKOS_IMPL_TAGGED_LOOP_L_2(tag, func, type, m_offset, extent, d - 1,      \
                                i2 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_4(tag, func, type, m_offset, extent, d, ...) \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {           \
    KOKKOS_IMPL_TAGGED_LOOP_L_3(tag, func, type, m_offset, extent, d - 1,      \
                                i3 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_5(tag, func, type, m_offset, extent, d, ...) \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {           \
    KOKKOS_IMPL_TAGGED_LOOP_L_4(tag, func, type, m_offset, extent, d - 1,      \
                                i4 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_6(tag, func, type, m_offset, extent, d, ...) \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {           \
    KOKKOS_IMPL_TAGGED_LOOP_L_5(tag, func, type, m_offset, extent, d - 1,      \
                                i5 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_7(tag, func, type, m_offset, extent, d, ...) \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {           \
    KOKKOS_IMPL_TAGGED_LOOP_L_6(tag, func, type, m_offset, extent, d - 1,      \
                                i6 + m_offset[d], __VA_ARGS__)                 \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_8(tag, func, type, m_offset, extent, d, ...) \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {           \
    KOKKOS_IMPL_TAGGED_LOOP_L_7(tag, func, type, m_offset, extent, d - 1,      \
                                i7 + m_offset[d], __VA_ARGS__)                 \
  }

// Left vs Right
// TODO: rank not necessary to pass through, can hardcode the values
#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1(tag, func, type, is_left, m_offset, \
                                         extent, rank)                       \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[0]); ++i0) {         \
    KOKKOS_IMPL_TAGGED_APPLY(tag, func, i0 + m_offset[0])                    \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2(tag, func, type, is_left, m_offset,   \
                                         extent, rank)                         \
  if (is_left) {                                                               \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[rank - 1]); ++i1) {  \
      KOKKOS_IMPL_TAGGED_LOOP_L_1(tag, func, type, m_offset, extent, rank - 2, \
                                  i1 + m_offset[rank - 1])                     \
    }                                                                          \
  } else {                                                                     \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[0]); ++i1) {         \
      KOKKOS_IMPL_TAGGED_LOOP_R_1(tag, func, type, m_offset, extent, 1,        \
                                  i1 + m_offset[0])                            \
    }                                                                          \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3(tag, func, type, is_left, m_offset,   \
                                         extent, rank)                         \
  if (is_left) {                                                               \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[rank - 1]); ++i2) {  \
      KOKKOS_IMPL_TAGGED_LOOP_L_2(tag, func, type, m_offset, extent, rank - 2, \
                                  i2 + m_offset[rank - 1])                     \
    }                                                                          \
  } else {                                                                     \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[0]); ++i2) {         \
      KOKKOS_IMPL_TAGGED_LOOP_R_2(tag, func, type, m_offset, extent, 1,        \
                                  i2 + m_offset[0])                            \
    }                                                                          \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4(tag, func, type, is_left, m_offset,   \
                                         extent, rank)                         \
  if (is_left) {                                                               \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[rank - 1]); ++i3) {  \
      KOKKOS_IMPL_TAGGED_LOOP_L_3(tag, func, type, m_offset, extent, rank - 2, \
                                  i3 + m_offset[rank - 1])                     \
    }                                                                          \
  } else {                                                                     \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[0]); ++i3) {         \
      KOKKOS_IMPL_TAGGED_LOOP_R_3(tag, func, type, m_offset, extent, 1,        \
                                  i3 + m_offset[0])                            \
    }                                                                          \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5(tag, func, type, is_left, m_offset,   \
                                         extent, rank)                         \
  if (is_left) {                                                               \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[rank - 1]); ++i4) {  \
      KOKKOS_IMPL_TAGGED_LOOP_L_4(tag, func, type, m_offset, extent, rank - 2, \
                                  i4 + m_offset[rank - 1])                     \
    }                                                                          \
  } else {                                                                     \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[0]); ++i4) {         \
      KOKKOS_IMPL_TAGGED_LOOP_R_4(tag, func, type, m_offset, extent, 1,        \
                                  i4 + m_offset[0])                            \
    }                                                                          \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6(tag, func, type, is_left, m_offset,   \
                                         extent, rank)                         \
  if (is_left) {                                                               \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[rank - 1]); ++i5) {  \
      KOKKOS_IMPL_TAGGED_LOOP_L_5(tag, func, type, m_offset, extent, rank - 2, \
                                  i5 + m_offset[rank - 1])                     \
    }                                                                          \
  } else {                                                                     \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[0]); ++i5) {         \
      KOKKOS_IMPL_TAGGED_LOOP_R_5(tag, func, type, m_offset, extent, 1,        \
                                  i5 + m_offset[0])                            \
    }                                                                          \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7(tag, func, type, is_left, m_offset,   \
                                         extent, rank)                         \
  if (is_left) {                                                               \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[rank - 1]); ++i6) {  \
      KOKKOS_IMPL_TAGGED_LOOP_L_6(tag, func, type, m_offset, extent, rank - 2, \
                                  i6 + m_offset[rank - 1])                     \
    }                                                                          \
  } else {                                                                     \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[0]); ++i6) {         \
      KOKKOS_IMPL_TAGGED_LOOP_R_6(tag, func, type, m_offset, extent, 1,        \
                                  i6 + m_offset[0])                            \
    }                                                                          \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8(tag, func, type, is_left, m_offset,   \
                                         extent, rank)                         \
  if (is_left) {                                                               \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[rank - 1]); ++i7) {  \
      KOKKOS_IMPL_TAGGED_LOOP_L_7(tag, func, type, m_offset, extent, rank - 2, \
                                  i7 + m_offset[rank - 1])                     \
    }                                                                          \
  } else {                                                                     \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[0]); ++i7) {         \
      KOKKOS_IMPL_TAGGED_LOOP_R_7(tag, func, type, m_offset, extent, 1,        \
                                  i7 + m_offset[0])                            \
    }                                                                          \
  }

// Partial vs Full Tile
#define KOKKOS_IMPL_TAGGED_TILE_LOOP_1(tag, func, type, is_left, cond,        \
                                       m_offset, extent_full, extent_partial, \
                                       rank)                                  \
  if (cond) {                                                                 \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1(tag, func, type, is_left, m_offset,      \
                                     extent_full, rank)                       \
  } else {                                                                    \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1(tag, func, type, is_left, m_offset,      \
                                     extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_2(tag, func, type, is_left, cond,        \
                                       m_offset, extent_full, extent_partial, \
                                       rank)                                  \
  if (cond) {                                                                 \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2(tag, func, type, is_left, m_offset,      \
                                     extent_full, rank)                       \
  } else {                                                                    \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2(tag, func, type, is_left, m_offset,      \
                                     extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_3(tag, func, type, is_left, cond,        \
                                       m_offset, extent_full, extent_partial, \
                                       rank)                                  \
  if (cond) {                                                                 \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3(tag, func, type, is_left, m_offset,      \
                                     extent_full, rank)                       \
  } else {                                                                    \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3(tag, func, type, is_left, m_offset,      \
                                     extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_4(tag, func, type, is_left, cond,        \
                                       m_offset, extent_full, extent_partial, \
                                       rank)                                  \
  if (cond) {                                                                 \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4(tag, func, type, is_left, m_offset,      \
                                     extent_full, rank)                       \
  } else {                                                                    \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4(tag, func, type, is_left, m_offset,      \
                                     extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_5(tag, func, type, is_left, cond,        \
                                       m_offset, extent_full, extent_partial, \
                                       rank)                                  \
  if (cond) {                                                                 \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5(tag, func, type, is_left, m_offset,      \
                                     extent_full, rank)                       \
  } else {                                                                    \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5(tag, func, type, is_left, m_offset,      \
                                     extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_6(tag, func, type, is_left, cond,        \
                                       m_offset, extent_full, extent_partial, \
                                       rank)                                  \
  if (cond) {                                                                 \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6(tag, func, type, is_left, m_offset,      \
                                     extent_full, rank)                       \
  } else {                                                                    \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6(tag, func, type, is_left, m_offset,      \
                                     extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_7(tag, func, type, is_left, cond,        \
                                       m_offset, extent_full, extent_partial, \
                                       rank)                                  \
  if (cond) {                                                                 \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7(tag, func, type, is_left, m_offset,      \
                                     extent_full, rank)                       \
  } else {                                                                    \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7(tag, func, type, is_left, m_offset,      \
                                     extent_partial, rank)                    \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_8(tag, func, type, is_left, cond,        \
                                       m_offset, extent_full, extent_partial, \
                                       rank)                                  \
  if (cond) {                                                                 \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8(tag, func, type, is_left, m_offset,      \
                                     extent_full, rank)                       \
  } else {                                                                    \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8(tag, func, type, is_left, m_offset,      \
                                     extent_partial, rank)                    \
  }

// parallel_reduce, tagged
// Reduction version
#define KOKKOS_IMPL_TAGGED_APPLY_REDUX(val, tag, func, ...) \
  func(tag, __VA_ARGS__, val);

// LayoutRight
// d = 0 to start
#define KOKKOS_IMPL_TAGGED_LOOP_R_1_REDUX(val, tag, func, type, m_offset, \
                                          extent, d, ...)                 \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                             \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {      \
    KOKKOS_IMPL_TAGGED_APPLY_REDUX(val, tag, func, __VA_ARGS__,           \
                                   i0 + m_offset[d])                      \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_2_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_1_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i1 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_3_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_2_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i2 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_4_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_3_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i3 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_5_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_4_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i4 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_6_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_5_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i5 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_7_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_6_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i6 + m_offset[d])   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_R_8_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {          \
    KOKKOS_IMPL_TAGGED_LOOP_R_7_REDUX(val, tag, func, type, m_offset, extent, \
                                      d + 1, __VA_ARGS__, i7 + m_offset[d])   \
  }

// LayoutLeft
// d = rank-1 to start
#define KOKKOS_IMPL_TAGGED_LOOP_L_1_REDUX(val, tag, func, type, m_offset, \
                                          extent, d, ...)                 \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                             \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[d]); ++i0) {      \
    KOKKOS_IMPL_TAGGED_APPLY_REDUX(val, tag, func, i0 + m_offset[d],      \
                                   __VA_ARGS__)                           \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_2_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i1 = (type)0; i1 < static_cast<type>(extent[d]); ++i1) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_1_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i1 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_3_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i2 = (type)0; i2 < static_cast<type>(extent[d]); ++i2) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_2_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i2 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_4_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i3 = (type)0; i3 < static_cast<type>(extent[d]); ++i3) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_3_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i3 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_5_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i4 = (type)0; i4 < static_cast<type>(extent[d]); ++i4) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_4_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i4 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_6_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i5 = (type)0; i5 < static_cast<type>(extent[d]); ++i5) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_5_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i5 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_7_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i6 = (type)0; i6 < static_cast<type>(extent[d]); ++i6) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_6_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i6 + m_offset[d], __VA_ARGS__)   \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_L_8_REDUX(val, tag, func, type, m_offset,     \
                                          extent, d, ...)                     \
  for (type i7 = (type)0; i7 < static_cast<type>(extent[d]); ++i7) {          \
    KOKKOS_IMPL_TAGGED_LOOP_L_7_REDUX(val, tag, func, type, m_offset, extent, \
                                      d - 1, i7 + m_offset[d], __VA_ARGS__)   \
  }

// Left vs Right
#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  KOKKOS_ENABLE_IVDEP_MDRANGE                                                 \
  for (type i0 = (type)0; i0 < static_cast<type>(extent[0]); ++i0) {          \
    KOKKOS_IMPL_TAGGED_APPLY_REDUX(val, tag, func, i0 + m_offset[0])          \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[rank - 1]); ++i1) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_1_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i1 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i1 = (type)0; i1 < static_cast<type>(extent[0]); ++i1) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_1_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i1 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[rank - 1]); ++i2) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_2_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i2 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i2 = (type)0; i2 < static_cast<type>(extent[0]); ++i2) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_2_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i2 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[rank - 1]); ++i3) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_3_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i3 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i3 = (type)0; i3 < static_cast<type>(extent[0]); ++i3) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_3_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i3 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[rank - 1]); ++i4) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_4_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i4 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i4 = (type)0; i4 < static_cast<type>(extent[0]); ++i4) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_4_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i4 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[rank - 1]); ++i5) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_5_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i5 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i5 = (type)0; i5 < static_cast<type>(extent[0]); ++i5) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_5_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i5 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[rank - 1]); ++i6) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_6_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i6 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i6 = (type)0; i6 < static_cast<type>(extent[0]); ++i6) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_6_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i6 + m_offset[0])          \
    }                                                                         \
  }

#define KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8_REDUX(val, tag, func, type, is_left, \
                                               m_offset, extent, rank)        \
  if (is_left) {                                                              \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[rank - 1]); ++i7) { \
      KOKKOS_IMPL_TAGGED_LOOP_L_7_REDUX(val, tag, func, type, m_offset,       \
                                        extent, rank - 2,                     \
                                        i7 + m_offset[rank - 1])              \
    }                                                                         \
  } else {                                                                    \
    for (type i7 = (type)0; i7 < static_cast<type>(extent[0]); ++i7) {        \
      KOKKOS_IMPL_TAGGED_LOOP_R_7_REDUX(val, tag, func, type, m_offset,       \
                                        extent, 1, i7 + m_offset[0])          \
    }                                                                         \
  }

// Partial vs Full Tile
#define KOKKOS_IMPL_TAGGED_TILE_LOOP_1_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_2_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_3_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_4_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_5_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_6_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_7_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

#define KOKKOS_IMPL_TAGGED_TILE_LOOP_8_REDUX(val, tag, func, type, is_left, \
                                             cond, m_offset, extent_full,   \
                                             extent_partial, rank)          \
  if (cond) {                                                               \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_full, rank)     \
  } else {                                                                  \
    KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8_REDUX(val, tag, func, type, is_left,   \
                                           m_offset, extent_partial, rank)  \
  }

// end tagged macros

// Structs for calling loops
template <int Rank, bool IsLeft, typename IType, typename Tagged,
          typename Enable = void>
struct Tile_Loop_Type;

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<1, IsLeft, IType, void, void> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_1(func, IType, IsLeft, cond, offset, a, b, 1);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_1_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 1);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<2, IsLeft, IType, void, void> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_2(func, IType, IsLeft, cond, offset, a, b, 2);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_2_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 2);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<3, IsLeft, IType, void, void> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_3(func, IType, IsLeft, cond, offset, a, b, 3);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_3_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 3);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<4, IsLeft, IType, void, void> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_4(func, IType, IsLeft, cond, offset, a, b, 4);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_4_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 4);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<5, IsLeft, IType, void, void> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_5(func, IType, IsLeft, cond, offset, a, b, 5);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_5_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 5);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<6, IsLeft, IType, void, void> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_6(func, IType, IsLeft, cond, offset, a, b, 6);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_6_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 6);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<7, IsLeft, IType, void, void> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_7(func, IType, IsLeft, cond, offset, a, b, 7);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_7_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 7);
  }
};

template <bool IsLeft, typename IType>
struct Tile_Loop_Type<8, IsLeft, IType, void, void> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_8(func, IType, IsLeft, cond, offset, a, b, 8);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TILE_LOOP_8_REDUX(value, func, IType, IsLeft, cond, offset, a,
                                  b, 8);
  }
};

// tagged versions

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<1, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void<Tagged>::value>> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_1(Tagged(), func, IType, IsLeft, cond, offset,
                                   a, b, 1);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_1_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 1);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<2, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void<Tagged>::value>> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_2(Tagged(), func, IType, IsLeft, cond, offset,
                                   a, b, 2);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_2_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 2);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<3, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void<Tagged>::value>> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_3(Tagged(), func, IType, IsLeft, cond, offset,
                                   a, b, 3);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_3_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 3);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<4, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void<Tagged>::value>> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_4(Tagged(), func, IType, IsLeft, cond, offset,
                                   a, b, 4);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_4_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 4);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<5, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void<Tagged>::value>> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_5(Tagged(), func, IType, IsLeft, cond, offset,
                                   a, b, 5);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_5_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 5);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<6, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void<Tagged>::value>> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_6(Tagged(), func, IType, IsLeft, cond, offset,
                                   a, b, 6);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_6_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 6);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<7, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void<Tagged>::value>> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_7(Tagged(), func, IType, IsLeft, cond, offset,
                                   a, b, 7);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_7_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 7);
  }
};

template <bool IsLeft, typename IType, typename Tagged>
struct Tile_Loop_Type<8, IsLeft, IType, Tagged,
                      std::enable_if_t<!std::is_void<Tagged>::value>> {
  template <typename Func, typename Offset, typename ExtentA, typename ExtentB>
  static void apply(Func const& func, bool cond, Offset const& offset,
                    ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_8(Tagged(), func, IType, IsLeft, cond, offset,
                                   a, b, 8);
  }

  template <typename ValType, typename Func, typename Offset, typename ExtentA,
            typename ExtentB>
  static void apply(ValType& value, Func const& func, bool cond,
                    Offset const& offset, ExtentA const& a, ExtentB const& b) {
    KOKKOS_IMPL_TAGGED_TILE_LOOP_8_REDUX(value, Tagged(), func, IType, IsLeft,
                                         cond, offset, a, b, 8);
  }
};
// end Structs for calling loops

template <typename RP, typename Functor, typename Tag = void,
          typename ValueType = void, typename Enable = void>
struct HostIterateTile;

// For ParallelFor
template <typename RP, typename Functor, typename Tag, typename ValueType>
struct HostIterateTile<RP, Functor, Tag, ValueType,
                       std::enable_if_t<std::is_void<ValueType>::value>> {
  using index_type = typename RP::index_type;
  using point_type = typename RP::point_type;

  using value_type = ValueType;

  inline HostIterateTile(RP const& rp, Functor const& func)
      : m_rp(rp), m_func(func) {}

  inline bool check_iteration_bounds(point_type& partial_tile,
                                     point_type& offset) const {
    bool is_full_tile = true;

    for (int i = 0; i < RP::rank; ++i) {
      if ((offset[i] + m_rp.m_tile[i]) <= m_rp.m_upper[i]) {
        partial_tile[i] = m_rp.m_tile[i];
      } else {
        is_full_tile = false;
        partial_tile[i] =
            (m_rp.m_upper[i] - 1 - offset[i]) == 0
                ? 1
                : (m_rp.m_upper[i] - m_rp.m_tile[i]) > 0
                      ? (m_rp.m_upper[i] - offset[i])
                      : (m_rp.m_upper[i] -
                         m_rp.m_lower[i]);  // when single tile encloses range
      }
    }

    return is_full_tile;
  }  // end check bounds

  template <int Rank>
  struct RankTag {
    using type = RankTag<Rank>;
    enum { value = (int)Rank };
  };

#if KOKKOS_ENABLE_NEW_LOOP_MACROS
  template <typename IType>
  inline void operator()(IType tile_idx) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    Tile_Loop_Type<RP::rank, (RP::inner_direction == Iterate::Left), index_type,
                   Tag>::apply(m_func, full_tile, m_offset, m_rp.m_tile,
                               m_tiledims);
  }

#else
  template <typename IType>
  inline void operator()(IType tile_idx) const {
    operator_impl(tile_idx, RankTag<RP::rank>());
  }
  // added due to compiler error when using sfinae to choose operator based on
  // rank w/ cuda+serial

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<2>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2L(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2L(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2R(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2R(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 2

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<3>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3L(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3L(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3R(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3R(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 3

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<4>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4L(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4L(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4R(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4R(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 4

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<5>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5L(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5L(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5R(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5R(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 5

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<6>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6L(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6L(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6R(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6R(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 6

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<7>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7L(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7L(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7R(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7R(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 7

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<8>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8L(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8L(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8R(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8R(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 8
#endif

  template <typename... Args>
  std::enable_if_t<(sizeof...(Args) == RP::rank && std::is_void<Tag>::value),
                   void>
  apply(Args&&... args) const {
    m_func(args...);
  }

  template <typename... Args>
  std::enable_if_t<(sizeof...(Args) == RP::rank && !std::is_void<Tag>::value),
                   void>
  apply(Args&&... args) const {
    m_func(m_tag, args...);
  }

  RP const m_rp;
  Functor const m_func;
  std::conditional_t<std::is_void<Tag>::value, int, Tag> m_tag;
};

// For ParallelReduce
// ValueType - scalar: For reductions
template <typename RP, typename Functor, typename Tag, typename ValueType>
struct HostIterateTile<RP, Functor, Tag, ValueType,
                       std::enable_if_t<!std::is_void<ValueType>::value &&
                                        !std::is_array<ValueType>::value>> {
  using index_type = typename RP::index_type;
  using point_type = typename RP::point_type;

  using value_type = ValueType;

  inline HostIterateTile(RP const& rp, Functor const& func)
      : m_rp(rp)  // Cuda 7.0 does not like braces...
        ,
        m_func(func) {
    // Errors due to braces rather than parenthesis for init (with cuda 7.0)
    //      /home/ndellin/kokkos/core/src/impl/KokkosExp_Host_IterateTile.hpp:1216:98:
    //      error: too many braces around initializer for ‘int’ [-fpermissive]
    //      /home/ndellin/kokkos/core/src/impl/KokkosExp_Host_IterateTile.hpp:1216:98:
    //      error: aggregate value used where an integer was expected
  }

  inline bool check_iteration_bounds(point_type& partial_tile,
                                     point_type& offset) const {
    bool is_full_tile = true;

    for (int i = 0; i < RP::rank; ++i) {
      if ((offset[i] + m_rp.m_tile[i]) <= m_rp.m_upper[i]) {
        partial_tile[i] = m_rp.m_tile[i];
      } else {
        is_full_tile = false;
        partial_tile[i] =
            (m_rp.m_upper[i] - 1 - offset[i]) == 0
                ? 1
                : (m_rp.m_upper[i] - m_rp.m_tile[i]) > 0
                      ? (m_rp.m_upper[i] - offset[i])
                      : (m_rp.m_upper[i] -
                         m_rp.m_lower[i]);  // when single tile encloses range
      }
    }

    return is_full_tile;
  }  // end check bounds

  template <int Rank>
  struct RankTag {
    using type = RankTag<Rank>;
    enum { value = (int)Rank };
  };

#if KOKKOS_ENABLE_NEW_LOOP_MACROS
  template <typename IType>
  inline void operator()(IType tile_idx, value_type& val) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    Tile_Loop_Type<RP::rank, (RP::inner_direction == Iterate::Left), index_type,
                   Tag>::apply(val, m_func.get_functor(), full_tile, m_offset,
                               m_rp.m_tile, m_tiledims);
  }

#else
  template <typename IType>
  inline void operator()(IType tile_idx) const {
    operator_impl(tile_idx, RankTag<RP::rank>());
  }
  // added due to compiler error when using sfinae to choose operator based on
  // rank

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<2>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2L(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2L(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2R(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2R(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 2

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<3>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3L(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3L(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3R(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3R(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 3

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<4>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4L(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4L(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4R(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4R(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 4

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<5>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5L(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5L(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5R(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5R(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 5

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<6>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6L(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6L(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6R(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6R(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 6

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<7>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7L(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7L(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7R(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7R(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 7

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<8>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8L(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8L(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8R(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8R(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 8

  template <typename... Args>
  std::enable_if_t<(sizeof...(Args) == RP::rank && std::is_void<Tag>::value),
                   void>
  apply(Args&&... args) const {
    m_func(args..., m_v);
  }

  template <typename... Args>
  std::enable_if_t<(sizeof...(Args) == RP::rank && !std::is_void<Tag>::value),
                   void>
  apply(Args&&... args) const {
    m_func(m_tag, args..., m_v);
  }
#endif

  RP const m_rp;
  Functor const m_func;
  std::conditional_t<std::is_void<Tag>::value, int, Tag> m_tag;
};

// For ParallelReduce
// Extra specialization for array reductions
// ValueType[]: For array reductions
template <typename RP, typename Functor, typename Tag, typename ValueType>
struct HostIterateTile<RP, Functor, Tag, ValueType,
                       std::enable_if_t<!std::is_void<ValueType>::value &&
                                        std::is_array<ValueType>::value>> {
  using index_type = typename RP::index_type;
  using point_type = typename RP::point_type;

  using value_type =
      std::remove_extent_t<ValueType>;  // strip away the
                                        // 'array-ness' [], only
                                        // underlying type remains

  inline HostIterateTile(RP const& rp, Functor const& func)
      : m_rp(rp)  // Cuda 7.0 does not like braces...
        ,
        m_func(func) {}

  inline bool check_iteration_bounds(point_type& partial_tile,
                                     point_type& offset) const {
    bool is_full_tile = true;

    for (int i = 0; i < RP::rank; ++i) {
      if ((offset[i] + m_rp.m_tile[i]) <= m_rp.m_upper[i]) {
        partial_tile[i] = m_rp.m_tile[i];
      } else {
        is_full_tile = false;
        partial_tile[i] =
            (m_rp.m_upper[i] - 1 - offset[i]) == 0
                ? 1
                : (m_rp.m_upper[i] - m_rp.m_tile[i]) > 0
                      ? (m_rp.m_upper[i] - offset[i])
                      : (m_rp.m_upper[i] -
                         m_rp.m_lower[i]);  // when single tile encloses range
      }
    }

    return is_full_tile;
  }  // end check bounds

  template <int Rank>
  struct RankTag {
    using type = RankTag<Rank>;
    enum { value = (int)Rank };
  };

#if KOKKOS_ENABLE_NEW_LOOP_MACROS
  template <typename IType>
  inline void operator()(IType tile_idx, value_type* val) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    Tile_Loop_Type<RP::rank, (RP::inner_direction == Iterate::Left), index_type,
                   Tag>::apply(val, m_func, full_tile, m_offset, m_rp.m_tile,
                               m_tiledims);
  }

#else
  template <typename IType>
  inline void operator()(IType tile_idx) const {
    operator_impl(tile_idx, RankTag<RP::rank>());
  }
  // added due to compiler error when using sfinae to choose operator based on
  // rank

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<2>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2L(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2L(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2R(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_2R(index_type, m_tiledims) { apply(LOOP_ARGS_2); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 2

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<3>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3L(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3L(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3R(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_3R(index_type, m_tiledims) { apply(LOOP_ARGS_3); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 3

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<4>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4L(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4L(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4R(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_4R(index_type, m_tiledims) { apply(LOOP_ARGS_4); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 4

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<5>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5L(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5L(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5R(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_5R(index_type, m_tiledims) { apply(LOOP_ARGS_5); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 5

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<6>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6L(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6L(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6R(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_6R(index_type, m_tiledims) { apply(LOOP_ARGS_6); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 6

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<7>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7L(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7L(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7R(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_7R(index_type, m_tiledims) { apply(LOOP_ARGS_7); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 7

  template <typename IType>
  inline void operator_impl(IType tile_idx, const RankTag<8>) const {
    point_type m_offset;
    point_type m_tiledims;

    if (RP::outer_direction == Iterate::Left) {
      for (int i = 0; i < RP::rank; ++i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    } else {
      for (int i = RP::rank - 1; i >= 0; --i) {
        m_offset[i] =
            (tile_idx % m_rp.m_tile_end[i]) * m_rp.m_tile[i] + m_rp.m_lower[i];
        tile_idx /= m_rp.m_tile_end[i];
      }
    }

    // Check if offset+tiledim in bounds - if not, replace tile dims with the
    // partial tile dims
    const bool full_tile = check_iteration_bounds(m_tiledims, m_offset);

    if (RP::inner_direction == Iterate::Left) {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8L(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8L(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      }
    }  // end Iterate::Left
    else {
      if (full_tile) {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8R(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      } else {
        //      #pragma simd
        KOKKOS_IMPL_LOOP_8R(index_type, m_tiledims) { apply(LOOP_ARGS_8); }
      }
    }  // end Iterate::Right

  }  // end op() rank == 8
  template <typename... Args>
  std::enable_if_t<(sizeof...(Args) == RP::rank && std::is_void<Tag>::value),
                   void>
  apply(Args&&... args) const {
    m_func(args..., m_v);
  }

  template <typename... Args>
  std::enable_if_t<(sizeof...(Args) == RP::rank && !std::is_void<Tag>::value),
                   void>
  apply(Args&&... args) const {
    m_func(m_tag, args..., m_v);
  }
#endif

  RP const m_rp;
  Functor const m_func;
  std::conditional_t<std::is_void<Tag>::value, int, Tag> m_tag;
};

// ------------------------------------------------------------------ //

#undef KOKKOS_ENABLE_NEW_LOOP_MACROS
#undef KOKKOS_IMPL_LOOP_1L
#undef KOKKOS_IMPL_LOOP_2L
#undef KOKKOS_IMPL_LOOP_3L
#undef KOKKOS_IMPL_LOOP_4L
#undef KOKKOS_IMPL_LOOP_5L
#undef KOKKOS_IMPL_LOOP_6L
#undef KOKKOS_IMPL_LOOP_7L
#undef KOKKOS_IMPL_LOOP_8L
#undef KOKKOS_IMPL_LOOP_1R
#undef KOKKOS_IMPL_LOOP_2R
#undef KOKKOS_IMPL_LOOP_3R
#undef KOKKOS_IMPL_LOOP_4R
#undef KOKKOS_IMPL_LOOP_5R
#undef KOKKOS_IMPL_LOOP_6R
#undef KOKKOS_IMPL_LOOP_7R
#undef KOKKOS_IMPL_LOOP_8R
#undef KOKKOS_IMPL_LOOP_ARGS_1
#undef KOKKOS_IMPL_LOOP_ARGS_2
#undef KOKKOS_IMPL_LOOP_ARGS_3
#undef KOKKOS_IMPL_LOOP_ARGS_4
#undef KOKKOS_IMPL_LOOP_ARGS_5
#undef KOKKOS_IMPL_LOOP_ARGS_6
#undef KOKKOS_IMPL_LOOP_ARGS_7
#undef KOKKOS_IMPL_LOOP_ARGS_8
#undef KOKKOS_IMPL_APPLY
#undef KOKKOS_IMPL_LOOP_R_1
#undef KOKKOS_IMPL_LOOP_R_2
#undef KOKKOS_IMPL_LOOP_R_3
#undef KOKKOS_IMPL_LOOP_R_4
#undef KOKKOS_IMPL_LOOP_R_5
#undef KOKKOS_IMPL_LOOP_R_6
#undef KOKKOS_IMPL_LOOP_R_7
#undef KOKKOS_IMPL_LOOP_R_8
#undef KOKKOS_IMPL_LOOP_L_1
#undef KOKKOS_IMPL_LOOP_L_2
#undef KOKKOS_IMPL_LOOP_L_3
#undef KOKKOS_IMPL_LOOP_L_4
#undef KOKKOS_IMPL_LOOP_L_5
#undef KOKKOS_IMPL_LOOP_L_6
#undef KOKKOS_IMPL_LOOP_L_7
#undef KOKKOS_IMPL_LOOP_L_8
#undef KOKKOS_IMPL_LOOP_LAYOUT_1
#undef KOKKOS_IMPL_LOOP_LAYOUT_2
#undef KOKKOS_IMPL_LOOP_LAYOUT_3
#undef KOKKOS_IMPL_LOOP_LAYOUT_4
#undef KOKKOS_IMPL_LOOP_LAYOUT_5
#undef KOKKOS_IMPL_LOOP_LAYOUT_6
#undef KOKKOS_IMPL_LOOP_LAYOUT_7
#undef KOKKOS_IMPL_LOOP_LAYOUT_8
#undef KOKKOS_IMPL_TILE_LOOP_1
#undef KOKKOS_IMPL_TILE_LOOP_2
#undef KOKKOS_IMPL_TILE_LOOP_3
#undef KOKKOS_IMPL_TILE_LOOP_4
#undef KOKKOS_IMPL_TILE_LOOP_5
#undef KOKKOS_IMPL_TILE_LOOP_6
#undef KOKKOS_IMPL_TILE_LOOP_7
#undef KOKKOS_IMPL_TILE_LOOP_8
#undef KOKKOS_IMPL_APPLY_REDUX
#undef KOKKOS_IMPL_LOOP_R_1_REDUX
#undef KOKKOS_IMPL_LOOP_R_2_REDUX
#undef KOKKOS_IMPL_LOOP_R_3_REDUX
#undef KOKKOS_IMPL_LOOP_R_4_REDUX
#undef KOKKOS_IMPL_LOOP_R_5_REDUX
#undef KOKKOS_IMPL_LOOP_R_6_REDUX
#undef KOKKOS_IMPL_LOOP_R_7_REDUX
#undef KOKKOS_IMPL_LOOP_R_8_REDUX
#undef KOKKOS_IMPL_LOOP_L_1_REDUX
#undef KOKKOS_IMPL_LOOP_L_2_REDUX
#undef KOKKOS_IMPL_LOOP_L_3_REDUX
#undef KOKKOS_IMPL_LOOP_L_4_REDUX
#undef KOKKOS_IMPL_LOOP_L_5_REDUX
#undef KOKKOS_IMPL_LOOP_L_6_REDUX
#undef KOKKOS_IMPL_LOOP_L_7_REDUX
#undef KOKKOS_IMPL_LOOP_L_8_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_1_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_2_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_3_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_4_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_5_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_6_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_7_REDUX
#undef KOKKOS_IMPL_LOOP_LAYOUT_8_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_1_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_2_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_3_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_4_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_5_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_6_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_7_REDUX
#undef KOKKOS_IMPL_TILE_LOOP_8_REDUX
#undef KOKKOS_IMPL_TAGGED_APPLY
#undef KOKKOS_IMPL_TAGGED_LOOP_R_1
#undef KOKKOS_IMPL_TAGGED_LOOP_R_2
#undef KOKKOS_IMPL_TAGGED_LOOP_R_3
#undef KOKKOS_IMPL_TAGGED_LOOP_R_4
#undef KOKKOS_IMPL_TAGGED_LOOP_R_5
#undef KOKKOS_IMPL_TAGGED_LOOP_R_6
#undef KOKKOS_IMPL_TAGGED_LOOP_R_7
#undef KOKKOS_IMPL_TAGGED_LOOP_R_8
#undef KOKKOS_IMPL_TAGGED_LOOP_L_1
#undef KOKKOS_IMPL_TAGGED_LOOP_L_2
#undef KOKKOS_IMPL_TAGGED_LOOP_L_3
#undef KOKKOS_IMPL_TAGGED_LOOP_L_4
#undef KOKKOS_IMPL_TAGGED_LOOP_L_5
#undef KOKKOS_IMPL_TAGGED_LOOP_L_6
#undef KOKKOS_IMPL_TAGGED_LOOP_L_7
#undef KOKKOS_IMPL_TAGGED_LOOP_L_8
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_1
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_2
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_3
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_4
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_5
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_6
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_7
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_8
#undef KOKKOS_IMPL_TAGGED_APPLY_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_1_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_2_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_3_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_4_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_5_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_6_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_7_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_R_8_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_1_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_2_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_3_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_4_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_5_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_6_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_7_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_L_8_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_1_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_2_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_3_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_4_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_5_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_6_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_7_REDUX
#undef KOKKOS_IMPL_TAGGED_LOOP_LAYOUT_8_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_1_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_2_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_3_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_4_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_5_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_6_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_7_REDUX
#undef KOKKOS_IMPL_TAGGED_TILE_LOOP_8_REDUX

}  // namespace Impl
}  // namespace Kokkos

#endif
