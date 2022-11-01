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

/// \file Kokkos_Layout.hpp
/// \brief Declaration of various \c MemoryLayout options.

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_3
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#else
KOKKOS_IMPL_WARNING("Including non-public Kokkos header files is not allowed.")
#endif
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

// To check for LayoutTiled
// This is to hide extra compile-time 'identifier' info within the LayoutTiled
// class by not relying on template specialization to include the ArgN*'s
template <typename LayoutTiledCheck, class Enable = void>
struct is_layouttiled : std::false_type {};

template <typename LayoutTiledCheck>
struct is_layouttiled<LayoutTiledCheck,
                      std::enable_if_t<LayoutTiledCheck::is_array_layout_tiled>>
    : std::true_type {};

namespace Experimental {

/// LayoutTiled
// Must have Rank >= 2
template <
    Kokkos::Iterate OuterP, Kokkos::Iterate InnerP, unsigned ArgN0,
    unsigned ArgN1, unsigned ArgN2 = 0, unsigned ArgN3 = 0, unsigned ArgN4 = 0,
    unsigned ArgN5 = 0, unsigned ArgN6 = 0, unsigned ArgN7 = 0,
    bool IsPowerOfTwo =
        (Kokkos::Impl::is_integral_power_of_two(ArgN0) &&
         Kokkos::Impl::is_integral_power_of_two(ArgN1) &&
         (Kokkos::Impl::is_integral_power_of_two(ArgN2) || (ArgN2 == 0)) &&
         (Kokkos::Impl::is_integral_power_of_two(ArgN3) || (ArgN3 == 0)) &&
         (Kokkos::Impl::is_integral_power_of_two(ArgN4) || (ArgN4 == 0)) &&
         (Kokkos::Impl::is_integral_power_of_two(ArgN5) || (ArgN5 == 0)) &&
         (Kokkos::Impl::is_integral_power_of_two(ArgN6) || (ArgN6 == 0)) &&
         (Kokkos::Impl::is_integral_power_of_two(ArgN7) || (ArgN7 == 0)))>
struct LayoutTiled {
  static_assert(IsPowerOfTwo,
                "LayoutTiled must be given power-of-two tile dimensions");

  using array_layout = LayoutTiled<OuterP, InnerP, ArgN0, ArgN1, ArgN2, ArgN3,
                                   ArgN4, ArgN5, ArgN6, ArgN7, IsPowerOfTwo>;
  static constexpr Iterate outer_pattern = OuterP;
  static constexpr Iterate inner_pattern = InnerP;

  enum { N0 = ArgN0 };
  enum { N1 = ArgN1 };
  enum { N2 = ArgN2 };
  enum { N3 = ArgN3 };
  enum { N4 = ArgN4 };
  enum { N5 = ArgN5 };
  enum { N6 = ArgN6 };
  enum { N7 = ArgN7 };

  size_t dimension[ARRAY_LAYOUT_MAX_RANK];

  enum : bool { is_extent_constructible = true };

  LayoutTiled(LayoutTiled const&) = default;
  LayoutTiled(LayoutTiled&&)      = default;
  LayoutTiled& operator=(LayoutTiled const&) = default;
  LayoutTiled& operator=(LayoutTiled&&) = default;

  KOKKOS_INLINE_FUNCTION
  explicit constexpr LayoutTiled(size_t argN0 = 0, size_t argN1 = 0,
                                 size_t argN2 = 0, size_t argN3 = 0,
                                 size_t argN4 = 0, size_t argN5 = 0,
                                 size_t argN6 = 0, size_t argN7 = 0)
      : dimension{argN0, argN1, argN2, argN3, argN4, argN5, argN6, argN7} {}

  friend bool operator==(const LayoutTiled& left, const LayoutTiled& right) {
    for (unsigned int rank = 0; rank < ARRAY_LAYOUT_MAX_RANK; ++rank)
      if (left.dimension[rank] != right.dimension[rank]) return false;
    return true;
  }

  friend bool operator!=(const LayoutTiled& left, const LayoutTiled& right) {
    return !(left == right);
  }
};

}  // namespace Experimental

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

template <unsigned ArgN0, unsigned ArgN1, unsigned ArgN2, unsigned ArgN3,
          unsigned ArgN4, unsigned ArgN5, unsigned ArgN6, unsigned ArgN7>
struct layout_iterate_type_selector<Kokkos::Experimental::LayoutTiled<
    Kokkos::Iterate::Left, Kokkos::Iterate::Left, ArgN0, ArgN1, ArgN2, ArgN3,
    ArgN4, ArgN5, ArgN6, ArgN7, true>> {
  static const Kokkos::Iterate outer_iteration_pattern = Kokkos::Iterate::Left;
  static const Kokkos::Iterate inner_iteration_pattern = Kokkos::Iterate::Left;
};

template <unsigned ArgN0, unsigned ArgN1, unsigned ArgN2, unsigned ArgN3,
          unsigned ArgN4, unsigned ArgN5, unsigned ArgN6, unsigned ArgN7>
struct layout_iterate_type_selector<Kokkos::Experimental::LayoutTiled<
    Kokkos::Iterate::Right, Kokkos::Iterate::Left, ArgN0, ArgN1, ArgN2, ArgN3,
    ArgN4, ArgN5, ArgN6, ArgN7, true>> {
  static const Kokkos::Iterate outer_iteration_pattern = Kokkos::Iterate::Right;
  static const Kokkos::Iterate inner_iteration_pattern = Kokkos::Iterate::Left;
};

template <unsigned ArgN0, unsigned ArgN1, unsigned ArgN2, unsigned ArgN3,
          unsigned ArgN4, unsigned ArgN5, unsigned ArgN6, unsigned ArgN7>
struct layout_iterate_type_selector<Kokkos::Experimental::LayoutTiled<
    Kokkos::Iterate::Left, Kokkos::Iterate::Right, ArgN0, ArgN1, ArgN2, ArgN3,
    ArgN4, ArgN5, ArgN6, ArgN7, true>> {
  static const Kokkos::Iterate outer_iteration_pattern = Kokkos::Iterate::Left;
  static const Kokkos::Iterate inner_iteration_pattern = Kokkos::Iterate::Right;
};

template <unsigned ArgN0, unsigned ArgN1, unsigned ArgN2, unsigned ArgN3,
          unsigned ArgN4, unsigned ArgN5, unsigned ArgN6, unsigned ArgN7>
struct layout_iterate_type_selector<Kokkos::Experimental::LayoutTiled<
    Kokkos::Iterate::Right, Kokkos::Iterate::Right, ArgN0, ArgN1, ArgN2, ArgN3,
    ArgN4, ArgN5, ArgN6, ArgN7, true>> {
  static const Kokkos::Iterate outer_iteration_pattern = Kokkos::Iterate::Right;
  static const Kokkos::Iterate inner_iteration_pattern = Kokkos::Iterate::Right;
};

}  // namespace Kokkos

#endif  // #ifndef KOKKOS_LAYOUT_HPP
