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
#pragma once

#include "macros.hpp"
#include "trait_backports.hpp"

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {
namespace detail {

//==============================================================================

template <class T, size_t Disambiguator = 0, class Enable = void>
struct no_unique_address_emulation {
  using stored_type = T;
  T m_v;
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T const &ref() const noexcept {
    return m_v;
  }
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T &ref() noexcept {
    return m_v;
  }
};

// Empty case
// This doesn't work if T is final, of course, but we're not using anything
// like that currently. That kind of thing could be added pretty easily though
template <class T, size_t Disambiguator>
struct no_unique_address_emulation<
    T, Disambiguator,
    std::enable_if_t<MDSPAN_IMPL_TRAIT(std::is_empty, T) &&
                // If the type isn't trivially destructible, its destructor
                // won't be called at the right time, so don't use this
                // specialization
                MDSPAN_IMPL_TRAIT(std::is_trivially_destructible, T)>> :
#ifdef MDSPAN_IMPL_COMPILER_MSVC
    // MSVC doesn't allow you to access public static member functions of a type
    // when you *happen* to privately inherit from that type.
    protected
#else
    // But we still want this to be private if possible so that we don't accidentally
    // access members of T directly rather than calling ref() first, which wouldn't
    // work if T happens to be stateful and thus we're using the unspecialized definition
    // of no_unique_address_emulation above.
    private
#endif
    T {
  using stored_type = T;
  MDSPAN_FORCE_INLINE_FUNCTION constexpr T const &ref() const noexcept {
    return *static_cast<T const *>(this);
  }
  MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 T &ref() noexcept {
    return *static_cast<T *>(this);
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr no_unique_address_emulation() noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr no_unique_address_emulation(
      no_unique_address_emulation const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr no_unique_address_emulation(
      no_unique_address_emulation &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED no_unique_address_emulation &
  operator=(no_unique_address_emulation const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED no_unique_address_emulation &
  operator=(no_unique_address_emulation &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~no_unique_address_emulation() noexcept = default;

  // Explicitly make this not a reference so that the copy or move
  // constructor still gets called.
  MDSPAN_INLINE_FUNCTION
  explicit constexpr no_unique_address_emulation(T const& v) noexcept : T(v) {}
  MDSPAN_INLINE_FUNCTION
  explicit constexpr no_unique_address_emulation(T&& v) noexcept : T(::std::move(v)) {}
};

//==============================================================================

} // end namespace detail
} // end namespace MDSPAN_IMPL_STANDARD_NAMESPACE
