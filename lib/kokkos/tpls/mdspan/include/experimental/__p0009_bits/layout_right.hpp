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
#include "extents.hpp"
#include "layout_stride.hpp"
#include "utility.hpp"
#if MDSPAN_HAS_CXX_17
#include "../__p2642_bits/layout_padded_fwd.hpp"
#endif

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {

//==============================================================================
template <class Extents>
class layout_right::mapping {
  public:
    using extents_type = Extents;
    using index_type = typename extents_type::index_type;
    using size_type = typename extents_type::size_type;
    using rank_type = typename extents_type::rank_type;
    using layout_type = layout_right;
  private:

    static_assert(detail::__is_extents_v<extents_type>,
                  MDSPAN_IMPL_STANDARD_NAMESPACE_STRING "::layout_right::mapping must be instantiated with a specialization of " MDSPAN_IMPL_STANDARD_NAMESPACE_STRING "::extents.");

    template <class>
    friend class mapping;

    // i0+(i1 + E(1)*(i2 + E(2)*i3))
    template <size_t r, size_t Rank>
    struct __rank_count {};

    template <size_t r, size_t Rank, class I, class... Indices>
    _MDSPAN_HOST_DEVICE
    constexpr index_type __compute_offset(
      index_type offset, __rank_count<r,Rank>, const I& i, Indices... idx) const {
      return __compute_offset(offset * __extents.extent(r) + i,__rank_count<r+1,Rank>(),  idx...);
    }

    template<class I, class ... Indices>
    _MDSPAN_HOST_DEVICE
    constexpr index_type __compute_offset(
      __rank_count<0,extents_type::rank()>, const I& i, Indices... idx) const {
      return __compute_offset(i,__rank_count<1,extents_type::rank()>(),idx...);
    }

    _MDSPAN_HOST_DEVICE
    constexpr index_type __compute_offset(size_t offset, __rank_count<extents_type::rank(), extents_type::rank()>) const {
      return static_cast<index_type>(offset);
    }

    _MDSPAN_HOST_DEVICE
    constexpr index_type __compute_offset(__rank_count<0,0>) const { return 0; }

  public:

    //--------------------------------------------------------------------------------

    MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mapping() noexcept = default;
    MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mapping(mapping const&) noexcept = default;

    _MDSPAN_HOST_DEVICE
    constexpr mapping(extents_type const& __exts) noexcept
      :__extents(__exts)
    { }

    MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (
        _MDSPAN_TRAIT(std::is_constructible, extents_type, OtherExtents)
      )
    )
    MDSPAN_CONDITIONAL_EXPLICIT((!std::is_convertible<OtherExtents, extents_type>::value)) // needs two () due to comma
    MDSPAN_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14
    mapping(mapping<OtherExtents> const& other) noexcept // NOLINT(google-explicit-constructor)
      :__extents(other.extents())
    {
       /*
        * TODO: check precondition
        * other.required_span_size() is a representable value of type index_type
        */
    }

    MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (
        _MDSPAN_TRAIT(std::is_constructible, extents_type, OtherExtents) &&
        (extents_type::rank() <= 1)
      )
    )
    MDSPAN_CONDITIONAL_EXPLICIT((!std::is_convertible<OtherExtents, extents_type>::value)) // needs two () due to comma
    MDSPAN_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14
    mapping(layout_left::mapping<OtherExtents> const& other) noexcept // NOLINT(google-explicit-constructor)
      :__extents(other.extents())
    {
       /*
        * TODO: check precondition
        * other.required_span_size() is a representable value of type index_type
        */
    }

    /**
     * Converting constructor from `layout_right_padded::mapping`.
     *
     * This overload participates in overload resolution only if _Mapping is a layout_right_padded mapping and
     * extents_type is constructible from _Mapping::extents_type.
     *
     * \note There is currently a difference from p2642r2, where this function is specified as taking
     * `layout_right_padded< padding_value >::mapping< Extents>`. However, this makes `padding_value` non-deducible.
     */
#if MDSPAN_HAS_CXX_17
    MDSPAN_TEMPLATE_REQUIRES(
        class _Mapping,
        /* requires */ (
        MDSPAN_IMPL_PROPOSED_NAMESPACE::detail::is_layout_right_padded_mapping<_Mapping>::value
        && std::is_constructible_v<extents_type, typename _Mapping::extents_type>))
    MDSPAN_CONDITIONAL_EXPLICIT((!std::is_convertible_v<typename _Mapping::extents_type, extents_type>))
    mapping(const _Mapping &__other) noexcept
        : __extents(__other.extents())
    {
      MDSPAN_IMPL_PROPOSED_NAMESPACE::detail::
          check_padded_layout_converting_constructor_mandates<
            extents_type, _Mapping>(detail::with_rank<extents_type::rank()>{});
      MDSPAN_IMPL_PROPOSED_NAMESPACE::detail::
          check_padded_layout_converting_constructor_preconditions<
            extents_type>(detail::with_rank<extents_type::rank()>{}, __other);
    }
#endif

    MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (
        _MDSPAN_TRAIT(std::is_constructible, extents_type, OtherExtents)
      )
    )
    MDSPAN_CONDITIONAL_EXPLICIT((extents_type::rank() > 0))
    MDSPAN_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14
    mapping(layout_stride::mapping<OtherExtents> const& other) noexcept // NOLINT(google-explicit-constructor)
      :__extents(other.extents())
    {
       /*
        * TODO: check precondition
        * other.required_span_size() is a representable value of type index_type
        */
       detail::validate_strides(detail::with_rank<extents_type::rank()>{}, layout_right{}, __extents, other);
    }

    MDSPAN_INLINE_FUNCTION_DEFAULTED _MDSPAN_CONSTEXPR_14_DEFAULTED mapping& operator=(mapping const&) noexcept = default;

    MDSPAN_INLINE_FUNCTION
    constexpr const extents_type& extents() const noexcept {
      return __extents;
    }

    MDSPAN_INLINE_FUNCTION
    constexpr index_type required_span_size() const noexcept {
      index_type value = 1;
      for(rank_type r=0; r != extents_type::rank(); ++r) value*=__extents.extent(r);
      return value;
    }

    //--------------------------------------------------------------------------------

    MDSPAN_TEMPLATE_REQUIRES(
      class ... Indices,
      /* requires */ (
      (sizeof...(Indices) == extents_type::rank()) &&
      (detail::are_valid_indices<index_type, Indices...>())
      )
    )
    _MDSPAN_HOST_DEVICE
    constexpr index_type operator()(Indices... idxs) const noexcept {
#if ! defined(NDEBUG)
      detail::check_all_indices(this->extents(), idxs...);
#endif // ! NDEBUG
      return __compute_offset(__rank_count<0, extents_type::rank()>(), static_cast<index_type>(idxs)...);
    }

    MDSPAN_INLINE_FUNCTION static constexpr bool is_always_unique() noexcept { return true; }
    MDSPAN_INLINE_FUNCTION static constexpr bool is_always_exhaustive() noexcept { return true; }
    MDSPAN_INLINE_FUNCTION static constexpr bool is_always_strided() noexcept { return true; }
    MDSPAN_INLINE_FUNCTION static constexpr bool is_unique() noexcept { return true; }
    MDSPAN_INLINE_FUNCTION static constexpr bool is_exhaustive() noexcept { return true; }
    MDSPAN_INLINE_FUNCTION static constexpr bool is_strided() noexcept { return true; }

    MDSPAN_INLINE_FUNCTION
    constexpr index_type stride(rank_type i) const noexcept
#if MDSPAN_HAS_CXX_20
      requires ( Extents::rank() > 0 )
#endif
    {
      index_type value = 1;
      for(rank_type r=extents_type::rank()-1; r>i; r--) value*=__extents.extent(r);
      return value;
    }

    MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ ( Extents::rank() == OtherExtents::rank())
    )
    MDSPAN_INLINE_FUNCTION
    friend constexpr bool operator==(mapping const& lhs, mapping<OtherExtents> const& rhs) noexcept {
      return lhs.extents() == rhs.extents();
    }

    // In C++ 20 the not equal exists if equal is found
#if !(MDSPAN_HAS_CXX_20)
    MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (Extents::rank() == OtherExtents::rank())
    )
    MDSPAN_INLINE_FUNCTION
    friend constexpr bool operator!=(mapping const& lhs, mapping<OtherExtents> const& rhs) noexcept {
      return lhs.extents() != rhs.extents();
    }
#endif

    // Not really public, but currently needed to implement fully constexpr useable submdspan:
    template<size_t N, class SizeType, size_t ... E, size_t ... Idx>
    constexpr index_type __get_stride(MDSPAN_IMPL_STANDARD_NAMESPACE::extents<SizeType, E...>,std::integer_sequence<size_t, Idx...>) const {
      return _MDSPAN_FOLD_TIMES_RIGHT((Idx>N? __extents.template __extent<Idx>():1),1);
    }
    template<size_t N>
    constexpr index_type __stride() const noexcept {
      return __get_stride<N>(__extents, std::make_index_sequence<extents_type::rank()>());
    }

private:
   _MDSPAN_NO_UNIQUE_ADDRESS extents_type __extents{};

   // [mdspan.submdspan.mapping], submdspan mapping specialization
   template<class... SliceSpecifiers>
   MDSPAN_INLINE_FUNCTION
   constexpr auto submdspan_mapping_impl(
       SliceSpecifiers... slices) const;

   template<class... SliceSpecifiers>
     friend constexpr auto submdspan_mapping(
       const mapping& src, SliceSpecifiers... slices) {
         return src.submdspan_mapping_impl(slices...);
     }
};

} // end namespace MDSPAN_IMPL_STANDARD_NAMESPACE

