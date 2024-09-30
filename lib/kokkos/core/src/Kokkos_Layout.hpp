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

/// \file Kokkos_Layout.hpp
/// \brief Declaration of various \c MemoryLayout options.

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_LAYOUT_HPP
#define KOKKOS_LAYOUT_HPP

#include <cstddef>
#include <impl/Kokkos_Traits.hpp>

namespace Kokkos {

enum { ARRAY_LAYOUT_MAX_RANK = 8 };

//----------------------------------------------------------------------------
/// \struct LayoutLeft
/// \brief Memory layout tag indicating left-to-right (Fortran scheme)
///   striding of multi-indices.
///
/// This is an example of a \c MemoryLayout template parameter of
/// View.  The memory layout describes how View maps from a
/// multi-index (i0, i1, ..., ik) to a memory location.
///
/// "Layout left" indicates a mapping where the leftmost index i0
/// refers to contiguous access, and strides increase for dimensions
/// going right from there (i1, i2, ...).  This layout imitates how
/// Fortran stores multi-dimensional arrays.  For the special case of
/// a two-dimensional array, "layout left" is also called "column
/// major."
struct LayoutLeft {
  //! Tag this class as a kokkos array layout
  using array_layout = LayoutLeft;

  size_t dimension[ARRAY_LAYOUT_MAX_RANK];

  enum : bool { is_extent_constructible = true };

  LayoutLeft(LayoutLeft const&) = default;
  LayoutLeft(LayoutLeft&&)      = default;
  LayoutLeft& operator=(LayoutLeft const&) = default;
  LayoutLeft& operator=(LayoutLeft&&) = default;

  KOKKOS_INLINE_FUNCTION
  explicit constexpr LayoutLeft(size_t N0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                size_t N1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                size_t N2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                size_t N3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                size_t N4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                size_t N5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                size_t N6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                size_t N7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
      : dimension{N0, N1, N2, N3, N4, N5, N6, N7} {}

  friend bool operator==(const LayoutLeft& left, const LayoutLeft& right) {
    for (unsigned int rank = 0; rank < ARRAY_LAYOUT_MAX_RANK; ++rank)
      if (left.dimension[rank] != right.dimension[rank]) return false;
    return true;
  }

  friend bool operator!=(const LayoutLeft& left, const LayoutLeft& right) {
    return !(left == right);
  }
};

//----------------------------------------------------------------------------
/// \struct LayoutRight
/// \brief Memory layout tag indicating right-to-left (C or
///   lexigraphical scheme) striding of multi-indices.
///
/// This is an example of a \c MemoryLayout template parameter of
/// View.  The memory layout describes how View maps from a
/// multi-index (i0, i1, ..., ik) to a memory location.
///
/// "Right layout" indicates a mapping where the rightmost index ik
/// refers to contiguous access, and strides increase for dimensions
/// going left from there.  This layout imitates how C stores
/// multi-dimensional arrays.  For the special case of a
/// two-dimensional array, "layout right" is also called "row major."
struct LayoutRight {
  //! Tag this class as a kokkos array layout
  using array_layout = LayoutRight;

  size_t dimension[ARRAY_LAYOUT_MAX_RANK];

  enum : bool { is_extent_constructible = true };

  LayoutRight(LayoutRight const&) = default;
  LayoutRight(LayoutRight&&)      = default;
  LayoutRight& operator=(LayoutRight const&) = default;
  LayoutRight& operator=(LayoutRight&&) = default;

  KOKKOS_INLINE_FUNCTION
  explicit constexpr LayoutRight(size_t N0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                 size_t N1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                 size_t N2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                 size_t N3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                 size_t N4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                 size_t N5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                 size_t N6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                                 size_t N7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
      : dimension{N0, N1, N2, N3, N4, N5, N6, N7} {}

  friend bool operator==(const LayoutRight& left, const LayoutRight& right) {
    for (unsigned int rank = 0; rank < ARRAY_LAYOUT_MAX_RANK; ++rank)
      if (left.dimension[rank] != right.dimension[rank]) return false;
    return true;
  }

  friend bool operator!=(const LayoutRight& left, const LayoutRight& right) {
    return !(left == right);
  }
};

//----------------------------------------------------------------------------
/// \struct LayoutStride
/// \brief  Memory layout tag indicated arbitrarily strided
///         multi-index mapping into contiguous memory.
struct LayoutStride {
  //! Tag this class as a kokkos array layout
  using array_layout = LayoutStride;

  size_t dimension[ARRAY_LAYOUT_MAX_RANK];
  size_t stride[ARRAY_LAYOUT_MAX_RANK];

  enum : bool { is_extent_constructible = false };

  LayoutStride(LayoutStride const&) = default;
  LayoutStride(LayoutStride&&)      = default;
  LayoutStride& operator=(LayoutStride const&) = default;
  LayoutStride& operator=(LayoutStride&&) = default;

  /** \brief  Compute strides from ordered dimensions.
   *
   *  Values of order uniquely form the set [0..rank)
   *  and specify ordering of the dimensions.
   *  Order = {0,1,2,...} is LayoutLeft
   *  Order = {...,2,1,0} is LayoutRight
   */
  template <typename iTypeOrder, typename iTypeDimen>
  KOKKOS_INLINE_FUNCTION static LayoutStride order_dimensions(
      int const rank, iTypeOrder const* const order,
      iTypeDimen const* const dimen) {
    LayoutStride tmp;
    // Verify valid rank order:
    int check_input = ARRAY_LAYOUT_MAX_RANK < rank ? 0 : int(1 << rank) - 1;
    for (int r = 0; r < ARRAY_LAYOUT_MAX_RANK; ++r) {
      tmp.dimension[r] = KOKKOS_IMPL_CTOR_DEFAULT_ARG;
      tmp.stride[r]    = 0;
    }
    for (int r = 0; r < rank; ++r) {
      check_input &= ~int(1 << order[r]);
    }
    if (0 == check_input) {
      size_t n = 1;
      for (int r = 0; r < rank; ++r) {
        tmp.stride[order[r]] = n;
        n *= (dimen[order[r]]);
        tmp.dimension[r] = dimen[r];
      }
    }
    return tmp;
  }

  KOKKOS_INLINE_FUNCTION
  explicit constexpr LayoutStride(
      size_t N0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG, size_t S0 = 0,
      size_t N1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG, size_t S1 = 0,
      size_t N2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG, size_t S2 = 0,
      size_t N3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG, size_t S3 = 0,
      size_t N4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG, size_t S4 = 0,
      size_t N5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG, size_t S5 = 0,
      size_t N6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG, size_t S6 = 0,
      size_t N7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG, size_t S7 = 0)
      : dimension{N0, N1, N2, N3, N4, N5, N6, N7}, stride{S0, S1, S2, S3,
                                                          S4, S5, S6, S7} {}

  friend bool operator==(const LayoutStride& left, const LayoutStride& right) {
    for (unsigned int rank = 0; rank < ARRAY_LAYOUT_MAX_RANK; ++rank)
      if (left.dimension[rank] != right.dimension[rank] ||
          left.stride[rank] != right.stride[rank])
        return false;
    return true;
  }

  friend bool operator!=(const LayoutStride& left, const LayoutStride& right) {
    return !(left == right);
  }
};

// ===================================================================================

//////////////////////////////////////////////////////////////////////////////////////

enum class Iterate {
  Default,
  Left,  // Left indices stride fastest
  Right  // Right indices stride fastest
};

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template <typename Layout, class Enable = void>
struct KOKKOS_DEPRECATED is_layouttiled : std::false_type {};
#endif

namespace Impl {
// For use with view_copy
template <typename... Layout>
struct layout_iterate_type_selector {
  static const Kokkos::Iterate outer_iteration_pattern =
      Kokkos::Iterate::Default;
  static const Kokkos::Iterate inner_iteration_pattern =
      Kokkos::Iterate::Default;
};

template <>
struct layout_iterate_type_selector<Kokkos::LayoutRight> {
  static const Kokkos::Iterate outer_iteration_pattern = Kokkos::Iterate::Right;
  static const Kokkos::Iterate inner_iteration_pattern = Kokkos::Iterate::Right;
};

template <>
struct layout_iterate_type_selector<Kokkos::LayoutLeft> {
  static const Kokkos::Iterate outer_iteration_pattern = Kokkos::Iterate::Left;
  static const Kokkos::Iterate inner_iteration_pattern = Kokkos::Iterate::Left;
};

template <>
struct layout_iterate_type_selector<Kokkos::LayoutStride> {
  static const Kokkos::Iterate outer_iteration_pattern =
      Kokkos::Iterate::Default;
  static const Kokkos::Iterate inner_iteration_pattern =
      Kokkos::Iterate::Default;
};
}  // namespace Impl

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template <typename... Layout>
using layout_iterate_type_selector KOKKOS_DEPRECATED =
    Impl::layout_iterate_type_selector<Layout...>;
#endif

}  // namespace Kokkos

#endif  // #ifndef KOKKOS_LAYOUT_HPP
