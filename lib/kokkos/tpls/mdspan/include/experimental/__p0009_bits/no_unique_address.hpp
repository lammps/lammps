/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2019) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
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

#pragma once

#include "macros.hpp"
#include "trait_backports.hpp"

namespace std {
namespace experimental {
namespace detail {

//==============================================================================

template <class _T, size_t _Disambiguator = 0, class _Enable = void>
struct __no_unique_address_emulation {
  using __stored_type = _T;
  _T __v;
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _T const &__ref() const noexcept {
    return __v;
  }
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _T &__ref() noexcept {
    return __v;
  }
};

// Empty case
// This doesn't work if _T is final, of course, but we're not using anything
// like that currently. That kind of thing could be added pretty easily though
template <class _T, size_t _Disambiguator>
struct __no_unique_address_emulation<
    _T, _Disambiguator,
    enable_if_t<_MDSPAN_TRAIT(is_empty, _T) &&
                // If the type isn't trivially destructible, its destructor
                // won't be called at the right time, so don't use this
                // specialization
                _MDSPAN_TRAIT(is_trivially_destructible, _T)>> : 
#ifdef _MDSPAN_COMPILER_MSVC
    // MSVC doesn't allow you to access public static member functions of a type
    // when you *happen* to privately inherit from that type.
    protected
#else
    // But we still want this to be private if possible so that we don't accidentally 
    // access members of _T directly rather than calling __ref() first, which wouldn't
    // work if _T happens to be stateful and thus we're using the unspecialized definition
    // of __no_unique_address_emulation above.
    private
#endif
    _T {
  using __stored_type = _T;
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _T const &__ref() const noexcept {
    return *static_cast<_T const *>(this);
  }
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _T &__ref() noexcept {
    return *static_cast<_T *>(this);
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __no_unique_address_emulation() noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __no_unique_address_emulation(
      __no_unique_address_emulation const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __no_unique_address_emulation(
      __no_unique_address_emulation &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __no_unique_address_emulation &
  operator=(__no_unique_address_emulation const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __no_unique_address_emulation &
  operator=(__no_unique_address_emulation &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~__no_unique_address_emulation() noexcept = default;

  // Explicitly make this not a reference so that the copy or move
  // constructor still gets called.
  MDSPAN_INLINE_FUNCTION
  explicit constexpr __no_unique_address_emulation(_T const& __v) noexcept : _T(__v) {}
  MDSPAN_INLINE_FUNCTION
  explicit constexpr __no_unique_address_emulation(_T&& __v) noexcept : _T(::std::move(__v)) {}
};

//==============================================================================

} // end namespace detail
} // end namespace experimental
} // end namespace std
