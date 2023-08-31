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

#include "dynamic_extent.hpp"
#include "trait_backports.hpp"
#include "maybe_static_value.hpp"
#include "standard_layout_static_array.hpp"
#include "type_list.hpp"

// Needs to be after the includes above to work with the single header generator
#if !_MDSPAN_PRESERVE_STANDARD_LAYOUT
#include <cstddef> // size_t
#include <utility> // integer_sequence
#include <array>

namespace std {
namespace experimental {
namespace detail {

//==============================================================================

template <class _T, _T _Val, bool _Mask> struct __mask_element {};

template <class _T, _T... _Result>
struct __mask_sequence_assign_op {
  template <_T _V>
  __mask_sequence_assign_op<_T, _Result..., _V>
  operator=(__mask_element<_T, _V, true>&&);
  template <_T _V>
  __mask_sequence_assign_op<_T, _Result...>
  operator=(__mask_element<_T, _V, false>&&);
  using __result = integer_sequence<_T, _Result...>;
};

template <class _Seq, class _Mask>
struct __mask_sequence;

template <class _T, _T... _Vals, bool... _Masks>
struct __mask_sequence<integer_sequence<_T, _Vals...>, integer_sequence<bool, _Masks...>>
{
  using type = typename decltype(
    _MDSPAN_FOLD_ASSIGN_LEFT(
      __mask_sequence_assign_op<_T>{}, /* = ... = */ __mask_element<_T, _Vals, _Masks>{}
    )
  )::__result;
};

//==============================================================================

template <class _T, class _static_t, class _Vals, _static_t __sentinal,
          class _Idxs, class _IdxsDynamic, class _IdxsDynamicIdxs>
class __partially_static_array_impl;

template <
  class _T, class _static_t,
  _static_t... __values_or_sentinals, _static_t __sentinal,
  size_t... _Idxs,
  size_t... _IdxsDynamic,
  size_t... _IdxsDynamicIdxs
>
class __partially_static_array_impl<
  _T,
  _static_t,
  integer_sequence<_static_t, __values_or_sentinals...>,
  __sentinal,
  integer_sequence<size_t, _Idxs...>,
  integer_sequence<size_t, _IdxsDynamic...>,
  integer_sequence<size_t, _IdxsDynamicIdxs...>
>
    : private __maybe_static_value<_T, _static_t, __values_or_sentinals, __sentinal,
                                   _Idxs>... {
private:

  template <size_t _N>
  using __base_n = typename __type_at<_N,
    __type_list<__maybe_static_value<_T, _static_t, __values_or_sentinals, __sentinal, _Idxs>...>
  >::type;

public:

  static constexpr auto __size = sizeof...(_Idxs);
  static constexpr auto __size_dynamic =
    _MDSPAN_FOLD_PLUS_RIGHT(static_cast<int>((__values_or_sentinals == __sentinal)), /* + ... + */ 0);

  //--------------------------------------------------------------------------

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __partially_static_array_impl() = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __partially_static_array_impl(
      __partially_static_array_impl const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __partially_static_array_impl(
      __partially_static_array_impl &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __partially_static_array_impl &
  operator=(__partially_static_array_impl const &) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __partially_static_array_impl &
  operator=(__partially_static_array_impl &&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~__partially_static_array_impl() noexcept = default;

  MDSPAN_INLINE_FUNCTION
  constexpr __partially_static_array_impl(
      __construct_psa_from_all_exts_values_tag_t,
      __repeated_with_idxs<_Idxs, _T> const &... __vals) noexcept
      : __base_n<_Idxs>(__base_n<_Idxs>{{__vals}})... {}

  MDSPAN_INLINE_FUNCTION
  constexpr __partially_static_array_impl(
      __construct_psa_from_dynamic_exts_values_tag_t,
      __repeated_with_idxs<_IdxsDynamicIdxs, _T> const &... __vals) noexcept
      : __base_n<_IdxsDynamic>(__base_n<_IdxsDynamic>{{__vals}})... {}

  MDSPAN_INLINE_FUNCTION constexpr explicit __partially_static_array_impl(
    array<_T, sizeof...(_Idxs)> const& __vals) noexcept
    : __partially_static_array_impl(
        __construct_psa_from_all_exts_values_tag,
        ::std::get<_Idxs>(__vals)...) {}

  // clang-format off
  MDSPAN_FUNCTION_REQUIRES(
    (MDSPAN_INLINE_FUNCTION constexpr explicit),
    __partially_static_array_impl,
    (array<_T, __size_dynamic> const &__vals), noexcept,
    /* requires */
      (sizeof...(_Idxs) != __size_dynamic)
  ): __partially_static_array_impl(
       __construct_psa_from_dynamic_exts_values_tag,
       ::std::get<_IdxsDynamicIdxs>(__vals)...) {}
  // clang-format on

  template <class _U, class _static_u, class _UValsSeq, _static_u __u_sentinal, class _UIdxsSeq,
            class _UIdxsDynamicSeq, class _UIdxsDynamicIdxsSeq>
  MDSPAN_INLINE_FUNCTION constexpr __partially_static_array_impl(
    __partially_static_array_impl<
      _U, _static_u, _UValsSeq, __u_sentinal, _UIdxsSeq,
     _UIdxsDynamicSeq, _UIdxsDynamicIdxsSeq> const &__rhs) noexcept
    : __partially_static_array_impl(
        __construct_psa_from_all_exts_values_tag,
        __rhs.template __get_n<_Idxs>()...) {}

  //--------------------------------------------------------------------------

  // See comment in the previous partial specialization for why this is
  // necessary.  Or just trust me that it's messy.
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr __partially_static_array_impl const &__enable_psa_conversion() const
  noexcept {
      return *this;
  }

  template <size_t _I>
  MDSPAN_FORCE_INLINE_FUNCTION constexpr _T __get_n() const noexcept {
    return static_cast<__base_n<_I> const*>(this)->__value();
  }

  template <class _U, size_t _I>
  MDSPAN_FORCE_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 void __set_n(_U&& __rhs) noexcept {
    static_cast<__base_n<_I>*>(this)->__set_value((_U&&)__rhs);
  }

  template <size_t _I, _static_t __default = __sentinal>
  MDSPAN_FORCE_INLINE_FUNCTION static constexpr _static_t
  __get_static_n() noexcept {
    return __base_n<_I>::__static_value == __sentinal ?
      __default : __base_n<_I>::__static_value;
  }

  MDSPAN_FORCE_INLINE_FUNCTION constexpr _T
  __get(size_t __n) const noexcept {
    return _MDSPAN_FOLD_PLUS_RIGHT(
      (_T(_Idxs == __n) * __get_n<_Idxs>()), /* + ... + */ _T(0)
    );
  }

};

//==============================================================================

template <class _T, class _static_t, class _ValSeq, _static_t __sentinal, class _Idxs = make_index_sequence<_ValSeq::size()>>
struct __partially_static_array_impl_maker;

template <
  class _T, class _static_t,  _static_t... _Vals, _static_t __sentinal, size_t... _Idxs
>
struct __partially_static_array_impl_maker<
  _T, _static_t, integer_sequence<_static_t, _Vals...>, __sentinal, integer_sequence<size_t, _Idxs...>
>
{
  using __dynamic_idxs = typename __mask_sequence<
    integer_sequence<size_t, _Idxs...>,
    integer_sequence<bool, (_Vals == __sentinal)...>
  >::type;
  using __impl_base =
    __partially_static_array_impl<_T, _static_t,
      integer_sequence<_static_t, _Vals...>,
      __sentinal, integer_sequence<size_t, _Idxs...>,
      __dynamic_idxs,
      make_index_sequence<__dynamic_idxs::size()>
    >;
};

template <class _T, class _static_t, class _ValsSeq, _static_t __sentinal = dynamic_extent>
class __partially_static_array_with_sentinal
  : public __partially_static_array_impl_maker<_T, _static_t, _ValsSeq, __sentinal>::__impl_base
{
private:
  using __base_t = typename __partially_static_array_impl_maker<_T, _static_t, _ValsSeq, __sentinal>::__impl_base;
public:
  using __base_t::__base_t;
};

//==============================================================================

template <class T, class _static_t, _static_t... __values_or_sentinals>
struct __partially_static_sizes :
  __partially_static_array_with_sentinal<
    T, _static_t, ::std::integer_sequence<_static_t, __values_or_sentinals...>>
{
private:
  using __base_t = __partially_static_array_with_sentinal<
    T, _static_t, ::std::integer_sequence<_static_t, __values_or_sentinals...>>;
public:
  using __base_t::__base_t;
  template <class _UTag>
  MDSPAN_FORCE_INLINE_FUNCTION constexpr __partially_static_sizes<T, _static_t, __values_or_sentinals...>
  __with_tag() const noexcept {
    return *this;
  }
};

// Tags are needed for the standard layout version, but not here
template <class T, class _static_t, _static_t... __values_or_sentinals>
using __partially_static_sizes_tagged = __partially_static_sizes<T, _static_t, __values_or_sentinals...>;

} // end namespace detail
} // end namespace experimental
} // end namespace std

#endif // !_MDSPAN_PRESERVE_STANDARD_LAYOUT
