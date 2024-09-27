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

#include <array>
#include <tuple>
#include <type_traits>
#include <utility> // index_sequence

// Suppress spurious warning with NVCC about no return statement.
// This is a known issue in NVCC and NVC++
// Depending on the CUDA and GCC version we need both the builtin
// and the diagnostic push. I tried really hard to find something shorter
// but no luck ...
#if defined __NVCC__
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress = implicit_return_from_non_void_function
#else
#ifdef __CUDA_ARCH__
#pragma diagnostic push
#pragma diag_suppress implicit_return_from_non_void_function
#endif
#endif
#elif defined __NVCOMPILER
#pragma diagnostic push
#pragma diag_suppress = implicit_return_from_non_void_function
#endif

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {
//******************************************
// Return type of submdspan_mapping overloads
//******************************************
template <class LayoutMapping> struct submdspan_mapping_result {
  _MDSPAN_NO_UNIQUE_ADDRESS LayoutMapping mapping{};
  size_t offset;
};

namespace detail {
// We use const Slice& and not Slice&& because the various
// submdspan_mapping_impl overloads use their slices arguments
// multiple times.  This makes perfect forwarding not useful, but we
// still don't want to pass those (possibly of size 64 x 3 bits)
// objects by value.
template <class IndexType, class Slice>
MDSPAN_INLINE_FUNCTION constexpr bool
one_slice_out_of_bounds(const IndexType &ext, const Slice &slice) {
  using common_t =
      std::common_type_t<decltype(detail::first_of(slice)), IndexType>;
  return static_cast<common_t>(detail::first_of(slice)) ==
         static_cast<common_t>(ext);
}

template <size_t... RankIndices, class IndexType, size_t... Exts,
          class... Slices>
MDSPAN_INLINE_FUNCTION constexpr bool
any_slice_out_of_bounds_helper(std::index_sequence<RankIndices...>,
                               const extents<IndexType, Exts...> &exts,
                               const Slices &... slices) {
  return _MDSPAN_FOLD_OR(
      (one_slice_out_of_bounds(exts.extent(RankIndices), slices)));
}

template <class IndexType, size_t... Exts, class... Slices>
MDSPAN_INLINE_FUNCTION constexpr bool
any_slice_out_of_bounds(const extents<IndexType, Exts...> &exts,
                        const Slices &... slices) {
  return any_slice_out_of_bounds_helper(
      std::make_index_sequence<sizeof...(Slices)>(), exts, slices...);
}

// constructs sub strides
template <class SrcMapping, class... slice_strides, size_t... InvMapIdxs>
MDSPAN_INLINE_FUNCTION constexpr auto construct_sub_strides(
    const SrcMapping &src_mapping, std::index_sequence<InvMapIdxs...>,
    const std::tuple<slice_strides...> &slices_stride_factor) {
  using index_type = typename SrcMapping::index_type;
  return std::array<typename SrcMapping::index_type, sizeof...(InvMapIdxs)>{
      (static_cast<index_type>(src_mapping.stride(InvMapIdxs)) *
       static_cast<index_type>(std::get<InvMapIdxs>(slices_stride_factor)))...};
}

template<class SliceSpecifier, class IndexType>
struct is_range_slice {
  constexpr static bool value =
    std::is_same_v<SliceSpecifier, full_extent_t> ||
    std::is_convertible_v<SliceSpecifier,
                          std::tuple<IndexType, IndexType>>;
};

template<class SliceSpecifier, class IndexType>
constexpr bool is_range_slice_v = is_range_slice<SliceSpecifier, IndexType>::value;

template<class SliceSpecifier, class IndexType>
struct is_index_slice {
  constexpr static bool value = std::is_convertible_v<SliceSpecifier, IndexType>;
};

template<class SliceSpecifier, class IndexType>
constexpr bool is_index_slice_v = is_index_slice<SliceSpecifier, IndexType>::value;

} // namespace detail

//**********************************
// layout_left submdspan_mapping
//*********************************
namespace detail {

// Figure out whether to preserve layout_left
template <class IndexType, size_t SubRank, class IndexSequence,
          class... SliceSpecifiers>
struct deduce_layout_left_submapping;

template <class IndexType, size_t SubRank, size_t... Idx,
          class... SliceSpecifiers>
struct deduce_layout_left_submapping<
    IndexType, SubRank, std::index_sequence<Idx...>, SliceSpecifiers...> {

  using count_range = index_sequence_scan_impl<
      0u, (is_index_slice_v<SliceSpecifiers, IndexType> ? 0u : 1u)...>;

  constexpr static int gap_len =
      (((Idx > 0 && count_range::get(Idx) == 1 &&
         is_index_slice_v<SliceSpecifiers, IndexType>)
            ? 1
            : 0) +
       ... + 0);

  MDSPAN_INLINE_FUNCTION
  constexpr static bool layout_left_value() {
    // Use layout_left for rank 0
    if constexpr (SubRank == 0) {
      return true;
    // Use layout_left for rank 1 result if leftmost slice specifier is range like
    } else if constexpr (SubRank == 1) {
      return ((Idx > 0 || is_range_slice_v<SliceSpecifiers, IndexType>)&&...);
    } else {
      // Preserve if leftmost SubRank-1 slices are full_extent_t and
      // the slice at idx Subrank - 1 is a range and
      // for idx > SubRank the slice is an index
      return ((((Idx <  SubRank - 1) && std::is_same_v<SliceSpecifiers, full_extent_t>) ||
               ((Idx == SubRank - 1) && is_range_slice_v<SliceSpecifiers, IndexType>) ||
               ((Idx >  SubRank - 1) && is_index_slice_v<SliceSpecifiers, IndexType>)) && ...);
    }
#if defined(__NVCC__) && !defined(__CUDA_ARCH__) && defined(__GNUC__)
    __builtin_unreachable();
#endif
  }

  MDSPAN_INLINE_FUNCTION
  constexpr static bool layout_left_padded_value() {
    // Technically could also keep layout_left_padded for SubRank==0
    // and SubRank==1 with leftmost slice specifier being a contiguous range
    // but we intercept these cases separately

    // In all other cases:
    // leftmost slice must be range
    // then there can be a gap with index slices
    // then SubRank - 2 full_extent slices
    // then another range slice
    // then more index slices
    // e.g. R I I I F F F R I I for obtaining a rank-5 from a rank-10
    return ((((Idx == 0)                                       && is_range_slice_v<SliceSpecifiers, IndexType>) ||
             ((Idx > 0 && Idx <= gap_len)                     && is_index_slice_v<SliceSpecifiers, IndexType>) ||
             ((Idx > gap_len && Idx < gap_len + SubRank - 1) && std::is_same_v<SliceSpecifiers, full_extent_t>) || 
             ((Idx == gap_len + SubRank - 1)                  && is_range_slice_v<SliceSpecifiers, IndexType>) ||
             ((Idx >  gap_len + SubRank - 1)                  && is_index_slice_v<SliceSpecifiers, IndexType>)) && ... );
  }
};

// We are reusing the same thing for layout_left and layout_left_padded
// For layout_left as source StaticStride is static_extent(0)
template<class Extents, size_t NumGaps, size_t StaticStride>
struct compute_s_static_layout_left {
  // Neither StaticStride nor any of the provided extents can be zero.
  // StaticStride can never be zero, the static_extents we are looking at are associated with
  // integral slice specifiers - which wouldn't be valid for zero extent
  template<size_t ... Idx>
  MDSPAN_INLINE_FUNCTION
  static constexpr size_t value(std::index_sequence<Idx...>) {
    size_t val = ((Idx>0 && Idx<=NumGaps ? (Extents::static_extent(Idx) == dynamic_extent?0:Extents::static_extent(Idx)) : 1) * ... * (StaticStride == dynamic_extent?0:StaticStride));
    return val == 0?dynamic_extent:val;
  }
};

} // namespace detail

// Actual submdspan mapping call
template <class Extents>
template <class... SliceSpecifiers>
MDSPAN_INLINE_FUNCTION constexpr auto
layout_left::mapping<Extents>::submdspan_mapping_impl(
    SliceSpecifiers... slices) const {

  // compute sub extents
  using src_ext_t = Extents;
  auto dst_ext = submdspan_extents(extents(), slices...);
  using dst_ext_t = decltype(dst_ext);

  // figure out sub layout type
  using deduce_layout = detail::deduce_layout_left_submapping<
      typename dst_ext_t::index_type, dst_ext_t::rank(),
      std::make_index_sequence<src_ext_t::rank()>,
      SliceSpecifiers...>;

  // Figure out if any slice's lower bound equals the corresponding extent.
  // If so, bypass evaluating the layout mapping.  This fixes LWG Issue 4060.
  const bool out_of_bounds =
      detail::any_slice_out_of_bounds(this->extents(), slices...);
  auto offset = static_cast<size_t>(
      out_of_bounds ? this->required_span_size()
                    : this->operator()(detail::first_of(slices)...));

  if constexpr (deduce_layout::layout_left_value()) {
    // layout_left case
    using dst_mapping_t = typename layout_left::template mapping<dst_ext_t>;
    return submdspan_mapping_result<dst_mapping_t>{dst_mapping_t(dst_ext),
                                                   offset};
  } else if constexpr (deduce_layout::layout_left_padded_value()) {
    constexpr size_t S_static = MDSPAN_IMPL_STANDARD_NAMESPACE::detail::compute_s_static_layout_left<Extents, deduce_layout::gap_len, Extents::static_extent(0)>::value(std::make_index_sequence<Extents::rank()>());
    using dst_mapping_t = typename MDSPAN_IMPL_PROPOSED_NAMESPACE::layout_left_padded<S_static>::template mapping<dst_ext_t>;
    return submdspan_mapping_result<dst_mapping_t>{
        dst_mapping_t(dst_ext, stride(1 + deduce_layout::gap_len)), offset};
  } else {
    // layout_stride case
    using dst_mapping_t = typename layout_stride::mapping<dst_ext_t>;
    auto inv_map = detail::inv_map_rank(std::integral_constant<size_t, 0>(),
                                        std::index_sequence<>(), slices...);
    return submdspan_mapping_result<dst_mapping_t> {
      dst_mapping_t(dst_ext,
                    detail::construct_sub_strides(
                        *this, inv_map,
// HIP needs deduction guides to have markups so we need to be explicit
// NVCC 11.0 has a bug with deduction guide here, tested that 11.2 does not have
// the issue But Clang-CUDA also doesn't accept the use of deduction guide so
// disable it for CUDA altogether
#if defined(_MDSPAN_HAS_HIP) || defined(_MDSPAN_HAS_CUDA)
                        std::tuple<decltype(detail::stride_of(slices))...>{
                            detail::stride_of(slices)...})),
#else
                        std::tuple{detail::stride_of(slices)...})),
#endif
          offset
    };
  }
#if defined(__NVCC__) && !defined(__CUDA_ARCH__) && defined(__GNUC__)
  __builtin_unreachable();
#endif
}

template <size_t PaddingValue>
template <class Extents>
template <class... SliceSpecifiers>
MDSPAN_INLINE_FUNCTION constexpr auto
MDSPAN_IMPL_PROPOSED_NAMESPACE::layout_left_padded<PaddingValue>::mapping<Extents>::submdspan_mapping_impl(
    SliceSpecifiers... slices) const {

  // compute sub extents
  using src_ext_t = Extents;
  auto dst_ext = submdspan_extents(extents(), slices...);
  using dst_ext_t = decltype(dst_ext);

  if constexpr (Extents::rank() == 0) { // rank-0 case
    using dst_mapping_t = typename MDSPAN_IMPL_PROPOSED_NAMESPACE::layout_left_padded<PaddingValue>::template mapping<Extents>;
    return submdspan_mapping_result<dst_mapping_t>{*this, 0};
  } else {
    const bool out_of_bounds =
        MDSPAN_IMPL_STANDARD_NAMESPACE::detail::any_slice_out_of_bounds(this->extents(), slices...);
    auto offset = static_cast<size_t>(
        out_of_bounds ? this->required_span_size()
                    : this->operator()(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::first_of(slices)...));
    if constexpr (dst_ext_t::rank() == 0) { // result rank-0
      // The following for some reasons leads to compiler error later, while not using a typedef works:
      // Compilers: CUDA 11.2 with GCC 9.1
      //
      // using dst_mapping_t = typename layout_left::template mapping<dst_ext_t>;
      // return submdspan_mapping_result<dst_mapping_t>{dst_mapping_t{dst_ext}, offset};
      //
      // Error: submdspan_mapping.hpp:299:23: error: 'dst_mapping_t' does not name a type
      //         299 |         using dst_mapping_t = typename layout_left::template mapping<dst_ext_t>;
      // The same error is given (about dst_mapping_t not naming type) when a different name is used in 299:
      //        using dst_mapping_t2 = typename layout_left::template mapping<dst_ext_t>;

      return submdspan_mapping_result<typename layout_left::template mapping<dst_ext_t>>
             {typename layout_left::template mapping<dst_ext_t>{dst_ext}, offset};
    } else { // general case
      // Figure out if any slice's lower bound equals the corresponding extent.
      // If so, bypass evaluating the layout mapping.  This fixes LWG Issue 4060.
      // figure out sub layout type
      using deduce_layout = MDSPAN_IMPL_STANDARD_NAMESPACE::detail::deduce_layout_left_submapping<
        typename dst_ext_t::index_type, dst_ext_t::rank(),
        decltype(std::make_index_sequence<src_ext_t::rank()>()),
        SliceSpecifiers...>;

      if constexpr (deduce_layout::layout_left_value() && dst_ext_t::rank() == 1) { // getting rank-1 from leftmost
        using dst_mapping_t = typename layout_left::template mapping<dst_ext_t>;
        return submdspan_mapping_result<dst_mapping_t>{dst_mapping_t{dst_ext}, offset};
      } else if constexpr (deduce_layout::layout_left_padded_value()) { // can keep layout_left_padded
        constexpr size_t S_static = MDSPAN_IMPL_STANDARD_NAMESPACE::detail::compute_s_static_layout_left<Extents, deduce_layout::gap_len, static_padding_stride>::value(std::make_index_sequence<Extents::rank()>());
        using dst_mapping_t = typename MDSPAN_IMPL_PROPOSED_NAMESPACE::layout_left_padded<S_static>::template mapping<dst_ext_t>;
        return submdspan_mapping_result<dst_mapping_t>{
        dst_mapping_t(dst_ext, stride(1 + deduce_layout::gap_len)), offset};
      } else { // layout_stride
    auto inv_map = MDSPAN_IMPL_STANDARD_NAMESPACE::detail::inv_map_rank(std::integral_constant<size_t, 0>(),
                                        std::index_sequence<>(), slices...);
      using dst_mapping_t = typename layout_stride::template mapping<dst_ext_t>;
    return submdspan_mapping_result<dst_mapping_t> {
      dst_mapping_t(dst_ext,
                    MDSPAN_IMPL_STANDARD_NAMESPACE::detail::construct_sub_strides(
                        *this, inv_map,
// HIP needs deduction guides to have markups so we need to be explicit
// NVCC 11.0 has a bug with deduction guide here, tested that 11.2 does not have
// the issue But Clang-CUDA also doesn't accept the use of deduction guide so
// disable it for CUDA alltogether
#if defined(_MDSPAN_HAS_HIP) || defined(_MDSPAN_HAS_CUDA)
                        std::tuple<decltype(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::stride_of(slices))...>{
                            MDSPAN_IMPL_STANDARD_NAMESPACE::detail::stride_of(slices)...})),
#else
                        std::tuple{MDSPAN_IMPL_STANDARD_NAMESPACE::detail::stride_of(slices)...})),
#endif
          offset
    };
      }
    }
  }


#if defined(__NVCC__) && !defined(__CUDA_ARCH__) && defined(__GNUC__)
  __builtin_unreachable();
#endif
}

//**********************************
// layout_right submdspan_mapping
//*********************************
namespace detail {

// Figure out whether to preserve layout_right
template <class IndexType, size_t SubRank, class IndexSequence,
          class... SliceSpecifiers>
struct deduce_layout_right_submapping;

template <class IndexType, size_t SubRank, size_t... Idx,
          class... SliceSpecifiers>
struct deduce_layout_right_submapping<
    IndexType, SubRank, std::index_sequence<Idx...>, SliceSpecifiers...> {

  static constexpr size_t Rank = sizeof...(Idx);
  using count_range = index_sequence_scan_impl<
      0u, (std::is_convertible_v<SliceSpecifiers, IndexType> ? 0u : 1u)...>;
  //__static_partial_sums<!std::is_convertible_v<SliceSpecifiers,
  // IndexType>...>;
  constexpr static int gap_len =
      (((Idx < Rank - 1 && count_range::get(Idx) == SubRank - 1 &&
         std::is_convertible_v<SliceSpecifiers, IndexType>)
            ? 1
            : 0) +
       ... + 0);

  MDSPAN_INLINE_FUNCTION
  constexpr static bool layout_right_value() {
    // Use layout_right for rank 0
    if constexpr (SubRank == 0) {
      return true;
    // Use layout_right for rank 1 result if rightmost slice specifier is range like
    } else if constexpr (SubRank == 1) {
      return ((Idx < Rank - 1 || is_range_slice_v<SliceSpecifiers, IndexType>)&&...);
    } else {
      // Preserve if rightmost SubRank-1 slices are full_extent_t and
      // the slice at idx Rank-Subrank is a range and
      // for idx < Rank - SubRank the slice is an index
      return ((((Idx >= Rank - SubRank) && std::is_same_v<SliceSpecifiers, full_extent_t>) ||
               ((Idx == Rank - SubRank) && is_range_slice_v<SliceSpecifiers, IndexType>) ||
               ((Idx <  Rank - SubRank) && is_index_slice_v<SliceSpecifiers, IndexType>)) && ...);
    }
#if defined(__NVCC__) && !defined(__CUDA_ARCH__) && defined(__GNUC__)
    __builtin_unreachable();
#endif
  }

  MDSPAN_INLINE_FUNCTION
  constexpr static bool layout_right_padded_value() {
    // Technically could also keep layout_right_padded for SubRank==0
    // and SubRank==1 with rightmost slice specifier being a contiguous range
    // but we intercept these cases separately

    // In all other cases:
    // rightmost slice must be range
    // then there can be a gap with index slices
    // then SubRank - 2 full_extent slices
    // then another range slice
    // then more index slices
    // e.g. I I R F F F I I I R for obtaining a rank-5 from a rank-10
    return ((((Idx == Rank - 1)                                               && is_range_slice_v<SliceSpecifiers, IndexType>) ||
             ((Idx >= Rank - gap_len - 1 && Idx < Rank - 1)                  && is_index_slice_v<SliceSpecifiers, IndexType>) ||
             ((Idx >  Rank - gap_len - SubRank && Idx < Rank - gap_len - 1) && std::is_same_v<SliceSpecifiers, full_extent_t>) ||
             ((Idx == Rank - gap_len - SubRank)                              && is_range_slice_v<SliceSpecifiers, IndexType>) ||
             ((Idx <  Rank - gap_len - SubRank)                              && is_index_slice_v<SliceSpecifiers, IndexType>)) && ... );
  }
};

// We are reusing the same thing for layout_right and layout_right_padded
// For layout_right as source StaticStride is static_extent(Rank-1)
template<class Extents, size_t NumGaps, size_t StaticStride>
struct compute_s_static_layout_right {
  // Neither StaticStride nor any of the provided extents can be zero.
  // StaticStride can never be zero, the static_extents we are looking at are associated with
  // integral slice specifiers - which wouldn't be valid for zero extent
  template<size_t ... Idx>
  MDSPAN_INLINE_FUNCTION
  static constexpr size_t value(std::index_sequence<Idx...>) {
    size_t val = ((Idx >= Extents::rank() - 1 - NumGaps && Idx < Extents::rank() - 1 ? (Extents::static_extent(Idx) == dynamic_extent?0:Extents::static_extent(Idx)) : 1) * ... * (StaticStride == dynamic_extent?0:StaticStride));
    return val == 0?dynamic_extent:val;
  }
};

} // namespace detail

// Actual submdspan mapping call
template <class Extents>
template <class... SliceSpecifiers>
MDSPAN_INLINE_FUNCTION constexpr auto
layout_right::mapping<Extents>::submdspan_mapping_impl(
    SliceSpecifiers... slices) const {

  // compute sub extents
  using src_ext_t = Extents;
  auto dst_ext = submdspan_extents(extents(), slices...);
  using dst_ext_t = decltype(dst_ext);

  // figure out sub layout type
  using deduce_layout = detail::deduce_layout_right_submapping<
      typename dst_ext_t::index_type, dst_ext_t::rank(),
      std::make_index_sequence<src_ext_t::rank()>,
      SliceSpecifiers...>;

  // Figure out if any slice's lower bound equals the corresponding extent.
  // If so, bypass evaluating the layout mapping.  This fixes LWG Issue 4060.
  const bool out_of_bounds =
      detail::any_slice_out_of_bounds(this->extents(), slices...);
  auto offset = static_cast<size_t>(
      out_of_bounds ? this->required_span_size()
                    : this->operator()(detail::first_of(slices)...));

  if constexpr (deduce_layout::layout_right_value()) {
    // layout_right case
    using dst_mapping_t = typename layout_right::mapping<dst_ext_t>;
    return submdspan_mapping_result<dst_mapping_t>{dst_mapping_t(dst_ext),
                                                   offset};
  } else if constexpr (deduce_layout::layout_right_padded_value()) {
    constexpr size_t S_static = MDSPAN_IMPL_STANDARD_NAMESPACE::detail::compute_s_static_layout_left<Extents, deduce_layout::gap_len, Extents::static_extent(Extents::rank() - 1)>::value(std::make_index_sequence<Extents::rank()>());
    using dst_mapping_t = typename MDSPAN_IMPL_PROPOSED_NAMESPACE::layout_right_padded<S_static>::template mapping<dst_ext_t>;
    return submdspan_mapping_result<dst_mapping_t>{
        dst_mapping_t(dst_ext,
                      stride(src_ext_t::rank() - 2 - deduce_layout::gap_len)),
        offset};
  } else {
    // layout_stride case
    using dst_mapping_t = typename layout_stride::mapping<dst_ext_t>;
    auto inv_map = detail::inv_map_rank(std::integral_constant<size_t, 0>(),
                                        std::index_sequence<>(), slices...);
    return submdspan_mapping_result<dst_mapping_t> {
      dst_mapping_t(dst_ext,
                    detail::construct_sub_strides(
                        *this, inv_map,
// HIP needs deduction guides to have markups so we need to be explicit
// NVCC 11.0 has a bug with deduction guide here, tested that 11.2 does not have
// the issue But Clang-CUDA also doesn't accept the use of deduction guide so
// disable it for CUDA altogether
#if defined(_MDSPAN_HAS_HIP) || defined(_MDSPAN_HAS_CUDA)
                        std::tuple<decltype(detail::stride_of(slices))...>{
                            detail::stride_of(slices)...})),
#else
                        std::tuple{detail::stride_of(slices)...})),
#endif
          offset
    };
  }
#if defined(__NVCC__) && !defined(__CUDA_ARCH__) && defined(__GNUC__)
  __builtin_unreachable();
#endif
}

template <size_t PaddingValue>
template <class Extents>
template <class... SliceSpecifiers>
MDSPAN_INLINE_FUNCTION constexpr auto
MDSPAN_IMPL_PROPOSED_NAMESPACE::layout_right_padded<PaddingValue>::mapping<Extents>::submdspan_mapping_impl(
    SliceSpecifiers... slices) const {

  // compute sub extents
  using src_ext_t = Extents;
  auto dst_ext = submdspan_extents(extents(), slices...);
  using dst_ext_t = decltype(dst_ext);

  if constexpr (Extents::rank() == 0) { // rank-0 case
    using dst_mapping_t = typename MDSPAN_IMPL_PROPOSED_NAMESPACE::layout_right_padded<PaddingValue>::template mapping<Extents>;
    return submdspan_mapping_result<dst_mapping_t>{*this, 0};
  } else {
    // Figure out if any slice's lower bound equals the corresponding extent.
    // If so, bypass evaluating the layout mapping.  This fixes LWG Issue 4060.
    // figure out sub layout type
    const bool out_of_bounds =
        MDSPAN_IMPL_STANDARD_NAMESPACE::detail::any_slice_out_of_bounds(this->extents(), slices...);
    auto offset = static_cast<size_t>(
        out_of_bounds ? this->required_span_size()
                    : this->operator()(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::first_of(slices)...));
    if constexpr (dst_ext_t::rank() == 0) { // result rank-0
      // Same issue as in layout_left_padded: see comment there
      // using dst_mapping_t = typename layout_right::template mapping<dst_ext_t>;
      // return submdspan_mapping_result<dst_mapping_t>{dst_mapping_t{dst_ext}, offset};
      return submdspan_mapping_result<typename layout_right::template mapping<dst_ext_t>>
        {typename layout_right::template mapping<dst_ext_t>{dst_ext}, offset};
    } else { // general case
      using deduce_layout = MDSPAN_IMPL_STANDARD_NAMESPACE::detail::deduce_layout_right_submapping<
        typename dst_ext_t::index_type, dst_ext_t::rank(),
        decltype(std::make_index_sequence<src_ext_t::rank()>()),
        SliceSpecifiers...>;

      if constexpr (deduce_layout::layout_right_value() && dst_ext_t::rank() == 1) { // getting rank-1 from rightmost
        using dst_mapping_t = typename layout_right::template mapping<dst_ext_t>;
        return submdspan_mapping_result<dst_mapping_t>{dst_mapping_t{dst_ext}, offset};
      } else if constexpr (deduce_layout::layout_right_padded_value()) { // can keep layout_right_padded
        constexpr size_t S_static = MDSPAN_IMPL_STANDARD_NAMESPACE::detail::compute_s_static_layout_right<Extents, deduce_layout::gap_len, static_padding_stride>::value(std::make_index_sequence<Extents::rank()>());
        using dst_mapping_t = typename MDSPAN_IMPL_PROPOSED_NAMESPACE::layout_right_padded<S_static>::template mapping<dst_ext_t>;
        return submdspan_mapping_result<dst_mapping_t>{
        dst_mapping_t(dst_ext, stride(Extents::rank() - 2 - deduce_layout::gap_len)), offset};
      } else { // layout_stride
    auto inv_map = MDSPAN_IMPL_STANDARD_NAMESPACE::detail::inv_map_rank(std::integral_constant<size_t, 0>(),
                                        std::index_sequence<>(), slices...);
      using dst_mapping_t = typename layout_stride::template mapping<dst_ext_t>;
    return submdspan_mapping_result<dst_mapping_t> {
      dst_mapping_t(dst_ext,
                    MDSPAN_IMPL_STANDARD_NAMESPACE::detail::construct_sub_strides(
                        *this, inv_map,
// HIP needs deduction guides to have markups so we need to be explicit
// NVCC 11.0 has a bug with deduction guide here, tested that 11.2 does not have
// the issue But Clang-CUDA also doesn't accept the use of deduction guide so
// disable it for CUDA alltogether
#if defined(_MDSPAN_HAS_HIP) || defined(_MDSPAN_HAS_CUDA)
                        std::tuple<decltype(MDSPAN_IMPL_STANDARD_NAMESPACE::detail::stride_of(slices))...>{
                            MDSPAN_IMPL_STANDARD_NAMESPACE::detail::stride_of(slices)...})),
#else
                        std::tuple{MDSPAN_IMPL_STANDARD_NAMESPACE::detail::stride_of(slices)...})),
#endif
          offset
    };
      }
    }
  }


#if defined(__NVCC__) && !defined(__CUDA_ARCH__) && defined(__GNUC__)
  __builtin_unreachable();
#endif
}

//**********************************
// layout_stride submdspan_mapping
//*********************************
template <class Extents>
template <class... SliceSpecifiers>
MDSPAN_INLINE_FUNCTION constexpr auto
layout_stride::mapping<Extents>::submdspan_mapping_impl(
    SliceSpecifiers... slices) const {
  auto dst_ext = submdspan_extents(extents(), slices...);
  using dst_ext_t = decltype(dst_ext);
  auto inv_map = detail::inv_map_rank(std::integral_constant<size_t, 0>(),
                                      std::index_sequence<>(), slices...);
  using dst_mapping_t = typename layout_stride::template mapping<dst_ext_t>;

  // Figure out if any slice's lower bound equals the corresponding extent.
  // If so, bypass evaluating the layout mapping.  This fixes LWG Issue 4060.
  const bool out_of_bounds =
      detail::any_slice_out_of_bounds(this->extents(), slices...);
  auto offset = static_cast<size_t>(
      out_of_bounds ? this->required_span_size()
                    : this->operator()(detail::first_of(slices)...));

  return submdspan_mapping_result<dst_mapping_t> {
    dst_mapping_t(dst_ext,
                  detail::construct_sub_strides(
                      *this, inv_map,
// HIP needs deduction guides to have markups so we need to be explicit
// NVCC 11.0 has a bug with deduction guide here, tested that 11.2 does not have
// the issue
#if defined(_MDSPAN_HAS_HIP) ||                                                \
    (defined(__NVCC__) &&                                                      \
     (__CUDACC_VER_MAJOR__ * 100 + __CUDACC_VER_MINOR__ * 10) < 1120)
                      std::tuple<decltype(detail::stride_of(slices))...>(
                          detail::stride_of(slices)...))),
#else
                      std::tuple(detail::stride_of(slices)...))),
#endif
        offset
  };
}

} // namespace MDSPAN_IMPL_STANDARD_NAMESPACE

#if defined __NVCC__
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#else
#ifdef __CUDA_ARCH__
#pragma diagnostic pop
#endif
#endif
#elif defined __NVCOMPILER
#pragma diagnostic pop
#endif
