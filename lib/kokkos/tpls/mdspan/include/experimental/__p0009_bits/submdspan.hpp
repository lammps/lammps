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

#include "mdspan.hpp"
#include "full_extent_t.hpp"
#include "dynamic_extent.hpp"
#include "layout_left.hpp"
#include "layout_right.hpp"
#include "layout_stride.hpp"
#include "macros.hpp"
#include "trait_backports.hpp"

#include <tuple> // std::apply
#include <utility> // std::pair

namespace std {
namespace experimental {

namespace detail {

template <size_t OldExtent, size_t OldStaticStride, class T>
struct __slice_wrap {
  T slice;
  size_t old_extent;
  size_t old_stride;
};

//--------------------------------------------------------------------------------

template <size_t OldExtent, size_t OldStaticStride>
MDSPAN_INLINE_FUNCTION constexpr
__slice_wrap<OldExtent, OldStaticStride, size_t>
__wrap_slice(size_t val, size_t ext, size_t stride) { return { val, ext, stride }; }

template <size_t OldExtent, size_t OldStaticStride, class IntegerType, IntegerType Value0>
MDSPAN_INLINE_FUNCTION constexpr
__slice_wrap<OldExtent, OldStaticStride, std::integral_constant<IntegerType, Value0>>
__wrap_slice(size_t val, size_t ext, std::integral_constant<IntegerType, Value0> stride)
{
#if MDSPAN_HAS_CXX_17
  if constexpr (std::is_signed_v<IntegerType>) {
    static_assert(Value0 >= IntegerType(0), "Invalid slice specifier");
  }
#endif // MDSPAN_HAS_CXX_17

  return { val, ext, stride };
}

template <size_t OldExtent, size_t OldStaticStride>
MDSPAN_INLINE_FUNCTION constexpr
__slice_wrap<OldExtent, OldStaticStride, full_extent_t>
__wrap_slice(full_extent_t val, size_t ext, size_t stride) { return { val, ext, stride }; }

// TODO generalize this to anything that works with std::get<0> and std::get<1>
template <size_t OldExtent, size_t OldStaticStride>
MDSPAN_INLINE_FUNCTION constexpr
__slice_wrap<OldExtent, OldStaticStride, std::tuple<size_t, size_t>>
__wrap_slice(std::tuple<size_t, size_t> const& val, size_t ext, size_t stride)
{
  return { val, ext, stride };
}

template <size_t OldExtent, size_t OldStaticStride,
	  class IntegerType0, IntegerType0 Value0,
	  class IntegerType1, IntegerType1 Value1>
MDSPAN_INLINE_FUNCTION constexpr
  __slice_wrap<OldExtent, OldStaticStride,
               std::tuple<std::integral_constant<IntegerType0, Value0>,
                          std::integral_constant<IntegerType1, Value1>>>
__wrap_slice(std::tuple<std::integral_constant<IntegerType0, Value0>, std::integral_constant<IntegerType1, Value1>> const& val, size_t ext, size_t stride)
{
  static_assert(Value1 >= Value0, "Invalid slice tuple");
  return { val, ext, stride };
}

//--------------------------------------------------------------------------------


// a layout right remains a layout right if it is indexed by 0 or more scalars,
// then optionally a pair and finally 0 or more all
template <
  // what we encountered until now preserves the layout right
  bool result=true,
  // we only encountered 0 or more scalars, no pair or all
  bool encountered_only_scalar=true
>
struct preserve_layout_right_analysis : integral_constant<bool, result> {
  using layout_type_if_preserved = layout_right;
  using encounter_pair = preserve_layout_right_analysis<
    // if we encounter a pair, the layout remains a layout right only if it was one before
    // and that only scalars were encountered until now
    result && encountered_only_scalar,
    // if we encounter a pair, we didn't encounter scalars only
    false
  >;
  using encounter_all = preserve_layout_right_analysis<
    // if we encounter a all, the layout remains a layout right if it was one before
    result,
    // if we encounter a all, we didn't encounter scalars only
    false
  >;
  using encounter_scalar = preserve_layout_right_analysis<
    // if we encounter a scalar, the layout remains a layout right only if it was one before
    // and that only scalars were encountered until now
    result && encountered_only_scalar,
    // if we encounter a scalar, the fact that we encountered scalars only doesn't change
    encountered_only_scalar
  >;
};

// a layout left remains a layout left if it is indexed by 0 or more all,
// then optionally a pair and finally 0 or more scalars
template <
  bool result=true,
  bool encountered_only_all=true
>
struct preserve_layout_left_analysis : integral_constant<bool, result> {
  using layout_type_if_preserved = layout_left;
  using encounter_pair = preserve_layout_left_analysis<
    // if we encounter a pair, the layout remains a layout left only if it was one before
    // and that only all were encountered until now
    result && encountered_only_all,
    // if we encounter a pair, we didn't encounter all only
    false
  >;
  using encounter_all = preserve_layout_left_analysis<
    // if we encounter a all, the layout remains a layout left only if it was one before
    // and that only all were encountered until now
    result && encountered_only_all,
    // if we encounter a all, the fact that we encountered scalars all doesn't change
    encountered_only_all
  >;
  using encounter_scalar = preserve_layout_left_analysis<
    // if we encounter a scalar, the layout remains a layout left if it was one before
    result,
    // if we encounter a scalar, we didn't encounter scalars only
    false
  >;
};

struct ignore_layout_preservation : std::integral_constant<bool, false> {
  using layout_type_if_preserved = void;
  using encounter_pair = ignore_layout_preservation;
  using encounter_all = ignore_layout_preservation;
  using encounter_scalar = ignore_layout_preservation;
};

template <class Layout>
struct preserve_layout_analysis
  : ignore_layout_preservation { };
template <>
struct preserve_layout_analysis<layout_right>
  : preserve_layout_right_analysis<> { };
template <>
struct preserve_layout_analysis<layout_left>
  : preserve_layout_left_analysis<> { };

//--------------------------------------------------------------------------------

template <
  class _IndexT,
  class _PreserveLayoutAnalysis,
  class _OffsetsArray=__partially_static_sizes<_IndexT, size_t>,
  class _ExtsArray=__partially_static_sizes<_IndexT, size_t>,
  class _StridesArray=__partially_static_sizes<_IndexT, size_t>,
  class = make_index_sequence<_OffsetsArray::__size>,
  class = make_index_sequence<_ExtsArray::__size>,
  class = make_index_sequence<_StridesArray::__size>
>
struct __assign_op_slice_handler;

/* clang-format: off */
template <
  class _IndexT,
  class _PreserveLayoutAnalysis,
  size_t... _Offsets,
  size_t... _Exts,
  size_t... _Strides,
  size_t... _OffsetIdxs,
  size_t... _ExtIdxs,
  size_t... _StrideIdxs>
struct __assign_op_slice_handler<
  _IndexT,
  _PreserveLayoutAnalysis,
  __partially_static_sizes<_IndexT, size_t, _Offsets...>,
  __partially_static_sizes<_IndexT, size_t, _Exts...>,
  __partially_static_sizes<_IndexT, size_t, _Strides...>,
  integer_sequence<size_t, _OffsetIdxs...>,
  integer_sequence<size_t, _ExtIdxs...>,
  integer_sequence<size_t, _StrideIdxs...>>
{
  // TODO remove this for better compiler performance
  static_assert(
    _MDSPAN_FOLD_AND((_Strides == dynamic_extent || _Strides > 0) /* && ... */),
    " "
  );
  static_assert(
    _MDSPAN_FOLD_AND((_Offsets == dynamic_extent || _Offsets >= 0) /* && ... */),
    " "
  );

  using __offsets_storage_t = __partially_static_sizes<_IndexT, size_t, _Offsets...>;
  using __extents_storage_t = __partially_static_sizes<_IndexT, size_t, _Exts...>;
  using __strides_storage_t = __partially_static_sizes<_IndexT, size_t, _Strides...>;
  __offsets_storage_t __offsets;
  __extents_storage_t __exts;
  __strides_storage_t __strides;

#ifdef __INTEL_COMPILER
#if __INTEL_COMPILER <= 1800
  MDSPAN_INLINE_FUNCTION constexpr __assign_op_slice_handler(__assign_op_slice_handler&& __other) noexcept
    : __offsets(::std::move(__other.__offsets)), __exts(::std::move(__other.__exts)), __strides(::std::move(__other.__strides))
  { }
  MDSPAN_INLINE_FUNCTION constexpr __assign_op_slice_handler(
    __offsets_storage_t&& __o,
    __extents_storage_t&& __e,
    __strides_storage_t&& __s
  ) noexcept
    : __offsets(::std::move(__o)), __exts(::std::move(__e)), __strides(::std::move(__s))
  { }
#endif
#endif

// Don't define this unless we need it; they have a cost to compile
#ifndef _MDSPAN_USE_RETURN_TYPE_DEDUCTION
  using __extents_type = ::std::experimental::extents<_IndexT, _Exts...>;
#endif

  // For size_t slice, skip the extent and stride, but add an offset corresponding to the value
  template <size_t _OldStaticExtent, size_t _OldStaticStride>
  MDSPAN_FORCE_INLINE_FUNCTION // NOLINT (misc-unconventional-assign-operator)
  _MDSPAN_CONSTEXPR_14 auto
  operator=(__slice_wrap<_OldStaticExtent, _OldStaticStride, size_t>&& __slice) noexcept
    -> __assign_op_slice_handler<
         _IndexT,
         typename _PreserveLayoutAnalysis::encounter_scalar,
         __partially_static_sizes<_IndexT, size_t, _Offsets..., dynamic_extent>,
         __partially_static_sizes<_IndexT, size_t, _Exts...>,
         __partially_static_sizes<_IndexT, size_t, _Strides...>/* intentional space here to work around ICC bug*/> {
    return {
      __partially_static_sizes<_IndexT, size_t, _Offsets..., dynamic_extent>(
        __construct_psa_from_all_exts_values_tag,
        __offsets.template __get_n<_OffsetIdxs>()..., __slice.slice),
      ::std::move(__exts),
      ::std::move(__strides)
    };
  }

  // Treat integral_constant slice like size_t slice, but with a compile-time offset.
  // The result's extents_type can't take advantage of that,
  // but it might help for specialized layouts.
  template <size_t _OldStaticExtent, size_t _OldStaticStride, class IntegerType, IntegerType Value0>
  MDSPAN_FORCE_INLINE_FUNCTION // NOLINT (misc-unconventional-assign-operator)
  _MDSPAN_CONSTEXPR_14 auto
  operator=(__slice_wrap<_OldStaticExtent, _OldStaticStride, std::integral_constant<IntegerType, Value0>>&&) noexcept
    -> __assign_op_slice_handler<
         _IndexT,
         typename _PreserveLayoutAnalysis::encounter_scalar,
         __partially_static_sizes<_IndexT, size_t, _Offsets..., Value0>,
         __partially_static_sizes<_IndexT, size_t, _Exts...>,
         __partially_static_sizes<_IndexT, size_t, _Strides...>/* intentional space here to work around ICC bug*/> {
#if MDSPAN_HAS_CXX_17
    if constexpr (std::is_signed_v<IntegerType>) {
      static_assert(Value0 >= IntegerType(0), "Invalid slice specifier");
    }
#endif // MDSPAN_HAS_CXX_17
    return {
      __partially_static_sizes<_IndexT, size_t, _Offsets..., Value0>(
        __construct_psa_from_all_exts_values_tag,
        __offsets.template __get_n<_OffsetIdxs>()..., size_t(Value0)),
      ::std::move(__exts),
      ::std::move(__strides)
    };
  }

  // For a std::full_extent, offset 0 and old extent
  template <size_t _OldStaticExtent, size_t _OldStaticStride>
  MDSPAN_FORCE_INLINE_FUNCTION // NOLINT (misc-unconventional-assign-operator)
  _MDSPAN_CONSTEXPR_14 auto
  operator=(__slice_wrap<_OldStaticExtent, _OldStaticStride, full_extent_t>&& __slice) noexcept
    -> __assign_op_slice_handler<
         _IndexT,
         typename _PreserveLayoutAnalysis::encounter_all,
         __partially_static_sizes<_IndexT, size_t, _Offsets..., 0>,
         __partially_static_sizes<_IndexT, size_t, _Exts..., _OldStaticExtent>,
         __partially_static_sizes<_IndexT, size_t, _Strides..., _OldStaticStride>/* intentional space here to work around ICC bug*/> {
    return {
      __partially_static_sizes<_IndexT, size_t, _Offsets..., 0>(
        __construct_psa_from_all_exts_values_tag,
        __offsets.template __get_n<_OffsetIdxs>()..., size_t(0)),
      __partially_static_sizes<_IndexT, size_t, _Exts..., _OldStaticExtent>(
        __construct_psa_from_all_exts_values_tag,
        __exts.template __get_n<_ExtIdxs>()..., __slice.old_extent),
      __partially_static_sizes<_IndexT, size_t, _Strides..., _OldStaticStride>(
        __construct_psa_from_all_exts_values_tag,
        __strides.template __get_n<_StrideIdxs>()..., __slice.old_stride)
    };
  }

  // For a std::tuple, add an offset and add a new dynamic extent (strides still preserved)
  template <size_t _OldStaticExtent, size_t _OldStaticStride>
  MDSPAN_FORCE_INLINE_FUNCTION // NOLINT (misc-unconventional-assign-operator)
  _MDSPAN_CONSTEXPR_14 auto
  operator=(__slice_wrap<_OldStaticExtent, _OldStaticStride, tuple<size_t, size_t>>&& __slice) noexcept
    -> __assign_op_slice_handler<
         _IndexT,
         typename _PreserveLayoutAnalysis::encounter_pair,
         __partially_static_sizes<_IndexT, size_t, _Offsets..., dynamic_extent>,
         __partially_static_sizes<_IndexT, size_t, _Exts..., dynamic_extent>,
         __partially_static_sizes<_IndexT, size_t, _Strides..., _OldStaticStride>/* intentional space here to work around ICC bug*/> {
    return {
      __partially_static_sizes<_IndexT, size_t, _Offsets..., dynamic_extent>(
        __construct_psa_from_all_exts_values_tag,
        __offsets.template __get_n<_OffsetIdxs>()..., ::std::get<0>(__slice.slice)),
      __partially_static_sizes<_IndexT, size_t, _Exts..., dynamic_extent>(
        __construct_psa_from_all_exts_values_tag,
        __exts.template __get_n<_ExtIdxs>()..., ::std::get<1>(__slice.slice) - ::std::get<0>(__slice.slice)),
      __partially_static_sizes<_IndexT, size_t, _Strides..., _OldStaticStride>(
        __construct_psa_from_all_exts_values_tag,
        __strides.template __get_n<_StrideIdxs>()..., __slice.old_stride)
    };
  }

  // For a std::tuple of two std::integral_constant, do something like
  // we did above for a tuple of two size_t, but make sure the
  // result's extents type make the values compile-time constants.
  template <size_t _OldStaticExtent, size_t _OldStaticStride,
	    class IntegerType0, IntegerType0 Value0,
	    class IntegerType1, IntegerType1 Value1>
  MDSPAN_FORCE_INLINE_FUNCTION // NOLINT (misc-unconventional-assign-operator)
  _MDSPAN_CONSTEXPR_14 auto
  operator=(__slice_wrap<_OldStaticExtent, _OldStaticStride, tuple<std::integral_constant<IntegerType0, Value0>, std::integral_constant<IntegerType1, Value1>>>&& __slice) noexcept
    -> __assign_op_slice_handler<
         _IndexT,
         typename _PreserveLayoutAnalysis::encounter_pair,
         __partially_static_sizes<_IndexT, size_t, _Offsets..., size_t(Value0)>,
         __partially_static_sizes<_IndexT, size_t, _Exts..., size_t(Value1 - Value0)>,
         __partially_static_sizes<_IndexT, size_t, _Strides..., _OldStaticStride>/* intentional space here to work around ICC bug*/> {
    static_assert(Value1 >= Value0, "Invalid slice specifier");
    return {
      // We're still turning the template parameters Value0 and Value1
      // into (constexpr) run-time values here.
      __partially_static_sizes<_IndexT, size_t, _Offsets..., size_t(Value0) > (
        __construct_psa_from_all_exts_values_tag,
        __offsets.template __get_n<_OffsetIdxs>()..., Value0),
      __partially_static_sizes<_IndexT, size_t, _Exts..., size_t(Value1 - Value0) > (
        __construct_psa_from_all_exts_values_tag,
        __exts.template __get_n<_ExtIdxs>()..., Value1 - Value0),
      __partially_static_sizes<_IndexT, size_t, _Strides..., _OldStaticStride>(
        __construct_psa_from_all_exts_values_tag,
        __strides.template __get_n<_StrideIdxs>()..., __slice.old_stride)
    };
  }

   // TODO defer instantiation of this?
  using layout_type = typename conditional<
    _PreserveLayoutAnalysis::value,
    typename _PreserveLayoutAnalysis::layout_type_if_preserved,
    layout_stride
  >::type;

  // TODO noexcept specification
  template <class NewLayout>
  MDSPAN_INLINE_FUNCTION
  _MDSPAN_DEDUCE_RETURN_TYPE_SINGLE_LINE(
    (
      _MDSPAN_CONSTEXPR_14 /* auto */
      _make_layout_mapping_impl(NewLayout) noexcept
    ),
    (
      /* not layout stride, so don't pass dynamic_strides */
      /* return */ typename NewLayout::template mapping<::std::experimental::extents<_IndexT, _Exts...>>(
        experimental::extents<_IndexT, _Exts...>::__make_extents_impl(::std::move(__exts))
      ) /* ; */
    )
  )

  MDSPAN_INLINE_FUNCTION
  _MDSPAN_DEDUCE_RETURN_TYPE_SINGLE_LINE(
    (
      _MDSPAN_CONSTEXPR_14 /* auto */
      _make_layout_mapping_impl(layout_stride) noexcept
    ),
    (
      /* return */ layout_stride::template mapping<::std::experimental::extents<_IndexT, _Exts...>>
        ::__make_mapping(::std::move(__exts), ::std::move(__strides)) /* ; */
    )
  )

  template <class OldLayoutMapping> // mostly for deferred instantiation, but maybe we'll use this in the future
  MDSPAN_INLINE_FUNCTION
  _MDSPAN_DEDUCE_RETURN_TYPE_SINGLE_LINE(
    (
      _MDSPAN_CONSTEXPR_14 /* auto */
      make_layout_mapping(OldLayoutMapping const&) noexcept
    ),
    (
      /* return */ this->_make_layout_mapping_impl(layout_type{}) /* ; */
    )
  )
};

//==============================================================================

#if _MDSPAN_USE_RETURN_TYPE_DEDUCTION
// Forking this because the C++11 version will be *completely* unreadable
template <class ET, class ST, size_t... Exts, class LP, class AP, class... SliceSpecs, size_t... Idxs>
MDSPAN_INLINE_FUNCTION
constexpr auto _submdspan_impl(
  integer_sequence<size_t, Idxs...>,
  mdspan<ET, std::experimental::extents<ST, Exts...>, LP, AP> const& src,
  SliceSpecs&&... slices
) noexcept
{
  using _IndexT = ST;
  auto _handled =
    _MDSPAN_FOLD_ASSIGN_LEFT(
      (
        detail::__assign_op_slice_handler<
          _IndexT,
          detail::preserve_layout_analysis<LP>
        >{
          __partially_static_sizes<_IndexT, size_t>{},
          __partially_static_sizes<_IndexT, size_t>{},
          __partially_static_sizes<_IndexT, size_t>{}
        }
      ),
        /* = ... = */
      detail::__wrap_slice<
        Exts, dynamic_extent
      >(
        slices, src.extents().template __extent<Idxs>(),
        src.mapping().stride(Idxs)
      )
    );

  size_t offset_size = src.mapping()(_handled.__offsets.template __get_n<Idxs>()...);
  auto offset_ptr = src.accessor().offset(src.data_handle(), offset_size);
  auto map = _handled.make_layout_mapping(src.mapping());
  auto acc_pol = typename AP::offset_policy(src.accessor());
  return mdspan<
    ET, remove_const_t<remove_reference_t<decltype(map.extents())>>,
        typename decltype(_handled)::layout_type, remove_const_t<remove_reference_t<decltype(acc_pol)>>
  >(
    std::move(offset_ptr), std::move(map), std::move(acc_pol)
  );
}
#else

template <class ET, class AP, class Src, class Handled, size_t... Idxs>
auto _submdspan_impl_helper(Src&& src, Handled&& h, std::integer_sequence<size_t, Idxs...>)
  -> mdspan<
       ET, typename Handled::__extents_type, typename Handled::layout_type, typename AP::offset_policy
     >
{
  return {
    src.accessor().offset(src.data_handle(), src.mapping()(h.__offsets.template __get_n<Idxs>()...)),
    h.make_layout_mapping(src.mapping()),
    typename AP::offset_policy(src.accessor())
  };
}

template <class ET, class ST, size_t... Exts, class LP, class AP, class... SliceSpecs, size_t... Idxs>
MDSPAN_INLINE_FUNCTION
_MDSPAN_DEDUCE_RETURN_TYPE_SINGLE_LINE(
  (
    constexpr /* auto */ _submdspan_impl(
      std::integer_sequence<size_t, Idxs...> seq,
      mdspan<ET, std::experimental::extents<ST, Exts...>, LP, AP> const& src,
      SliceSpecs&&... slices
    ) noexcept
  ),
  (
    /* return */ _submdspan_impl_helper<ET, AP>(
      src,
      _MDSPAN_FOLD_ASSIGN_LEFT(
        (
          detail::__assign_op_slice_handler<
            size_t,
            detail::preserve_layout_analysis<LP>
          >{
            __partially_static_sizes<ST, size_t>{},
            __partially_static_sizes<ST, size_t>{},
            __partially_static_sizes<ST, size_t>{}
          }
        ),
        /* = ... = */
        detail::__wrap_slice<
          Exts, dynamic_extent
        >(
          slices, src.extents().template __extent<Idxs>(), src.mapping().stride(Idxs)
        )
      ),
      seq
    ) /* ; */
  )
)

#endif

template <class T> struct _is_layout_stride : std::false_type { };
template<>
struct _is_layout_stride<
  layout_stride
> : std::true_type
{ };

} // namespace detail

//==============================================================================

MDSPAN_TEMPLATE_REQUIRES(
  class ET, class EXT, class LP, class AP, class... SliceSpecs,
  /* requires */ (
    (
      _MDSPAN_TRAIT(is_same, LP, layout_left)
        || _MDSPAN_TRAIT(is_same, LP, layout_right)
        || detail::_is_layout_stride<LP>::value
    ) &&
    _MDSPAN_FOLD_AND((
      _MDSPAN_TRAIT(is_convertible, SliceSpecs, size_t)
        || _MDSPAN_TRAIT(is_convertible, SliceSpecs, tuple<size_t, size_t>)
        || _MDSPAN_TRAIT(is_convertible, SliceSpecs, full_extent_t)
    ) /* && ... */) &&
    sizeof...(SliceSpecs) == EXT::rank()
  )
)
MDSPAN_INLINE_FUNCTION
_MDSPAN_DEDUCE_RETURN_TYPE_SINGLE_LINE(
  (
    constexpr submdspan(
      mdspan<ET, EXT, LP, AP> const& src, SliceSpecs... slices
    ) noexcept
  ),
  (
    /* return */
      detail::_submdspan_impl(std::make_index_sequence<sizeof...(SliceSpecs)>{}, src, slices...) /*;*/
  )
)
/* clang-format: on */

} // end namespace experimental
} // namespace std
