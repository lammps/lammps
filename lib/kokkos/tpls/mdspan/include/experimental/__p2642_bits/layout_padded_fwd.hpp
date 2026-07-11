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
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#pragma once

#include <cassert>
#include "../__p0009_bits/dynamic_extent.hpp"
#include "../__p0009_bits/utility.hpp"

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {
namespace MDSPAN_IMPL_PROPOSED_NAMESPACE {

template <size_t padding_value = dynamic_extent>
struct layout_left_padded {
  template <class Extents>
  class mapping;
};

template <size_t padding_value = dynamic_extent>
struct layout_right_padded {
  template <class Extents>
  class mapping;
};

namespace detail {
// The layout_padded_constants structs are only useful if rank > 1, otherwise they may wrap
template <class Layout, class ExtentsType>
struct layout_padded_constants;

template <class ExtentsType, size_t PaddingStride>
struct layout_padded_constants<layout_left_padded<PaddingStride>, ExtentsType>
{
  using rank_type = typename ExtentsType::rank_type;
  static constexpr rank_type padded_stride_idx = 1;
  static constexpr rank_type extent_to_pad_idx = 0;
};

template <class ExtentsType, size_t PaddingStride>
struct layout_padded_constants<layout_right_padded<PaddingStride>, ExtentsType>
{
  using rank_type = typename ExtentsType::rank_type;
  static constexpr rank_type padded_stride_idx = ExtentsType::rank() - 2;
  static constexpr rank_type extent_to_pad_idx = ExtentsType::rank() - 1;
};

template <class Layout>
struct is_layout_left_padded : std::false_type {};

template <size_t PaddingStride>
struct is_layout_left_padded<layout_left_padded<PaddingStride>> : std::true_type {};

template <class Mapping, class Enabled = void>
struct is_layout_left_padded_mapping : std::false_type {};

template <class Mapping>
struct is_layout_left_padded_mapping<Mapping,
  std::enable_if_t<std::is_same<Mapping, typename layout_left_padded<Mapping::padding_value>::template mapping<typename Mapping::extents_type>>::value>>
    : std::true_type {};

template <class Layout>
struct is_layout_right_padded : std::false_type {};

template <size_t PaddingStride>
struct is_layout_right_padded<layout_right_padded<PaddingStride>> : std::true_type {};

template <class Mapping, class Enabled = void>
struct is_layout_right_padded_mapping : std::false_type {};

template <class Mapping>
struct is_layout_right_padded_mapping<Mapping,
  std::enable_if_t<std::is_same<Mapping, typename layout_right_padded<Mapping::padding_value>::template mapping<typename Mapping::extents_type>>::value>>
    : std::true_type {};


template <class LayoutExtentsType, class PaddedLayoutMappingType>
MDSPAN_INLINE_FUNCTION
constexpr void check_padded_layout_converting_constructor_mandates(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::with_rank<0>) {}

template <class LayoutExtentsType, class PaddedLayoutMappingType>
MDSPAN_INLINE_FUNCTION
constexpr void check_padded_layout_converting_constructor_mandates(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::with_rank<1>) {}

template <class LayoutExtentsType, class PaddedLayoutMappingType, std::size_t N>
MDSPAN_INLINE_FUNCTION
constexpr void check_padded_layout_converting_constructor_mandates(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::with_rank<N>)
{
  using extents_type = typename PaddedLayoutMappingType::extents_type;
  constexpr auto padding_value = PaddedLayoutMappingType::padding_value;
  constexpr auto idx = layout_padded_constants<typename PaddedLayoutMappingType::layout_type, LayoutExtentsType >::extent_to_pad_idx;

  constexpr auto statically_determinable =
    (LayoutExtentsType::static_extent(idx) != dynamic_extent) &&
    (extents_type::static_extent(idx) != dynamic_extent) &&
    (padding_value != dynamic_extent);

  static_assert(!statically_determinable ||
                (padding_value == 0
                 ? LayoutExtentsType::static_extent(idx) == 0
                 : LayoutExtentsType::static_extent(idx) % padding_value == 0),
                "");
}

template <typename ExtentsType, typename OtherMapping>
MDSPAN_INLINE_FUNCTION
constexpr void check_padded_layout_converting_constructor_preconditions(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::with_rank<0>,
                                                                        const OtherMapping&) {}
template <typename ExtentsType, typename OtherMapping>
MDSPAN_INLINE_FUNCTION
constexpr void check_padded_layout_converting_constructor_preconditions(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::with_rank<1>,
                                                                        const OtherMapping&) {}
template <typename ExtentsType, typename OtherMapping, std::size_t N>
MDSPAN_INLINE_FUNCTION
constexpr void check_padded_layout_converting_constructor_preconditions(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::with_rank<N>,
                                                                        const OtherMapping &other_mapping) {
  constexpr auto padded_stride_idx =
    layout_padded_constants<typename OtherMapping::layout_type,
                            ExtentsType>::padded_stride_idx;
  constexpr auto extent_to_pad_idx = layout_padded_constants<typename OtherMapping::layout_type, ExtentsType>::extent_to_pad_idx;
  MDSPAN_IMPL_PRECONDITION(other_mapping.stride(padded_stride_idx) == other_mapping.extents().extent(extent_to_pad_idx));
}


}
}
}
