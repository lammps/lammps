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

#if !defined(_MDSPAN_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
#  include "no_unique_address.hpp"
#endif

namespace std {
namespace experimental {
namespace detail {

// For no unique address emulation, this is the case taken when neither are empty.
// For real `[[no_unique_address]]`, this case is always taken.
template <class _T, class _U, class _Enable = void> struct __compressed_pair {
  _MDSPAN_NO_UNIQUE_ADDRESS _T __t_val;
  _MDSPAN_NO_UNIQUE_ADDRESS _U __u_val;
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _T &__first() noexcept { return __t_val; }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _T const &__first() const noexcept {
    return __t_val;
  }
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _U &__second() noexcept { return __u_val; }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _U const &__second() const noexcept {
    return __u_val;
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair() noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair(__compressed_pair const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair(__compressed_pair &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __compressed_pair &
  operator=(__compressed_pair const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __compressed_pair &
  operator=(__compressed_pair &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~__compressed_pair() noexcept = default;
  template <class _TLike, class _ULike>
  MDSPAN_INLINE_FUNCTION constexpr __compressed_pair(_TLike &&__t, _ULike &&__u)
      : __t_val((_TLike &&) __t), __u_val((_ULike &&) __u) {}
};

#if !defined(_MDSPAN_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)

// First empty.
template <class _T, class _U>
struct __compressed_pair<
    _T, _U,
    enable_if_t<_MDSPAN_TRAIT(is_empty, _T) && !_MDSPAN_TRAIT(is_empty, _U)>>
    : private _T {
  _U __u_val;
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _T &__first() noexcept {
    return *static_cast<_T *>(this);
  }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _T const &__first() const noexcept {
    return *static_cast<_T const *>(this);
  }
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _U &__second() noexcept { return __u_val; }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _U const &__second() const noexcept {
    return __u_val;
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair() noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair(__compressed_pair const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair(__compressed_pair &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __compressed_pair &
  operator=(__compressed_pair const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __compressed_pair &
  operator=(__compressed_pair &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~__compressed_pair() noexcept = default;
  template <class _TLike, class _ULike>
  MDSPAN_INLINE_FUNCTION constexpr __compressed_pair(_TLike &&__t, _ULike &&__u)
      : _T((_TLike &&) __t), __u_val((_ULike &&) __u) {}
};

// Second empty.
template <class _T, class _U>
struct __compressed_pair<
    _T, _U,
    enable_if_t<!_MDSPAN_TRAIT(is_empty, _T) && _MDSPAN_TRAIT(is_empty, _U)>>
    : private _U {
  _T __t_val;
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _T &__first() noexcept { return __t_val; }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _T const &__first() const noexcept {
    return __t_val;
  }
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _U &__second() noexcept {
    return *static_cast<_U *>(this);
  }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _U const &__second() const noexcept {
    return *static_cast<_U const *>(this);
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair() noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair(__compressed_pair const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair(__compressed_pair &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __compressed_pair &
  operator=(__compressed_pair const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __compressed_pair &
  operator=(__compressed_pair &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~__compressed_pair() noexcept = default;

  template <class _TLike, class _ULike>
  MDSPAN_INLINE_FUNCTION constexpr __compressed_pair(_TLike &&__t, _ULike &&__u)
      : _U((_ULike &&) __u), __t_val((_TLike &&) __t) {}
};

// Both empty.
template <class _T, class _U>
struct __compressed_pair<
    _T, _U,
    enable_if_t<_MDSPAN_TRAIT(is_empty, _T) && _MDSPAN_TRAIT(is_empty, _U)>>
    // We need to use the __no_unique_address_emulation wrapper here to avoid
    // base class ambiguities.
#ifdef _MDSPAN_COMPILER_MSVC
// MSVC doesn't allow you to access public static member functions of a type
// when you *happen* to privately inherit from that type.
    : protected __no_unique_address_emulation<_T, 0>,
      protected __no_unique_address_emulation<_U, 1>
#else
    : private __no_unique_address_emulation<_T, 0>,
      private __no_unique_address_emulation<_U, 1>
#endif
{
  using __first_base_t = __no_unique_address_emulation<_T, 0>;
  using __second_base_t = __no_unique_address_emulation<_U, 1>;

  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _T &__first() noexcept {
    return this->__first_base_t::__ref();
  }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _T const &__first() const noexcept {
    return this->__first_base_t::__ref();
  }
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 _U &__second() noexcept {
    return this->__second_base_t::__ref();
  }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _U const &__second() const noexcept {
    return this->__second_base_t::__ref();
  }

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair() noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair(__compressed_pair const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __compressed_pair(__compressed_pair &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __compressed_pair &
  operator=(__compressed_pair const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __compressed_pair &
  operator=(__compressed_pair &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~__compressed_pair() noexcept = default;
  template <class _TLike, class _ULike>
  MDSPAN_INLINE_FUNCTION constexpr __compressed_pair(_TLike &&__t, _ULike &&__u) noexcept
    : __first_base_t(_T((_TLike &&) __t)),
      __second_base_t(_U((_ULike &&) __u))
  { }
};

#endif // !defined(_MDSPAN_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)

} // end namespace detail
} // end namespace experimental
} // end namespace std
