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
#include "layout_padded_fwd.hpp"
#include "../__p0009_bits/dynamic_extent.hpp"
#include "../__p0009_bits/extents.hpp"
#include "../__p0009_bits/mdspan.hpp"
#include "../__p0009_bits/layout_left.hpp"
#include "../__p0009_bits/layout_right.hpp"
#include "../__p0009_bits/layout_stride.hpp"
#include "../__p0009_bits/utility.hpp"

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {
namespace MDSPAN_IMPL_PROPOSED_NAMESPACE {

namespace detail {
template<class T>
MDSPAN_INLINE_FUNCTION
constexpr T
find_next_multiple(T alignment, T offset)
{
  if ( alignment == 0 ) {
    return T(0);
  } else {
    return ( ( offset + alignment - 1 ) / alignment) * alignment;
  }
}

template <class ExtentsType, size_t PaddingValue, size_t ExtentToPadIdx>
MDSPAN_INLINE_FUNCTION constexpr size_t get_actual_static_padding_value() {
  constexpr auto rank = ExtentsType::rank();

  if constexpr (rank <= typename ExtentsType::rank_type(1)) {
    return 0;
  } else if constexpr (PaddingValue != dynamic_extent &&
                       ExtentsType::static_extent(ExtentToPadIdx) !=
                           dynamic_extent) {
    static_assert(
        (PaddingValue != 0) ||
            (ExtentsType::static_extent(ExtentToPadIdx) == 0),
        "padding stride can be 0 only if "
        "extents_type::static_extent(extent-to-pad) is 0 or dynamic_extent");
    return find_next_multiple(PaddingValue,
                                ExtentsType::static_extent(ExtentToPadIdx));
  } else {
    return dynamic_extent;
  }
  // Missing return statement warning from NVCC and ICC
#if (defined(__NVCC__) || defined(__INTEL_COMPILER)) && !defined(__NVCOMPILER)
  return 0;
#endif
}

template <size_t PaddingValue, typename Extents, size_t ExtentToPadIdx, size_t Rank, typename Enabled = void>
struct static_array_type_for_padded_extent
{
  static constexpr size_t padding_value = PaddingValue;
  using index_type = typename Extents::index_type;
  using extents_type = Extents;
  using type = ::MDSPAN_IMPL_STANDARD_NAMESPACE::detail::maybe_static_array<
      index_type, size_t, dynamic_extent,
      ::MDSPAN_IMPL_STANDARD_NAMESPACE::MDSPAN_IMPL_PROPOSED_NAMESPACE::detail::get_actual_static_padding_value<extents_type, PaddingValue,
                                                ExtentToPadIdx>()>;
};

template <size_t PaddingValue, typename Extents, size_t ExtentToPadIdx, size_t Rank>
struct static_array_type_for_padded_extent<PaddingValue, Extents,
                                             ExtentToPadIdx, Rank, std::enable_if_t<Rank <= 1>> {
  using index_type = typename Extents::index_type;
  using extents_type = Extents;
  using type =
      ::MDSPAN_IMPL_STANDARD_NAMESPACE::detail::maybe_static_array<
          index_type, size_t, dynamic_extent, 0>;
};

template <size_t PaddingValue, typename Extents, size_t ExtentToPadIdx>
struct padded_extent {
  static constexpr size_t padding_value = PaddingValue;
  using index_type = typename Extents::index_type;
  using extents_type = Extents;
  using static_array_type = typename static_array_type_for_padded_extent<
      padding_value, Extents, ExtentToPadIdx, Extents::rank()>::type;

  MDSPAN_INLINE_FUNCTION
  static constexpr auto static_value() { return static_array_type::static_value(0); }

  MDSPAN_INLINE_FUNCTION
  static constexpr static_array_type
  init_padding(const Extents &exts) {
    if constexpr ((Extents::rank() > 1) && (padding_value == dynamic_extent)) {
      return {exts.extent(ExtentToPadIdx)};
    } else {
      return init_padding(exts, padding_value);
    }
    // Missing return statement warning from NVCC and ICC
#if (defined(__NVCC__) || defined(__INTEL_COMPILER)) && !defined(__NVCOMPILER)
    return {};
#endif
  }

  MDSPAN_INLINE_FUNCTION static constexpr static_array_type
  init_padding([[maybe_unused]] const Extents &exts,
               [[maybe_unused]] index_type pv) {
    if constexpr (Extents::rank() > 1) {
      return {find_next_multiple(pv,
                                   exts.extent(ExtentToPadIdx))};
    } else {
      return {};
    }
    // Missing return statement warning from NVCC and ICC
#if (defined(__NVCC__) || defined(__INTEL_COMPILER)) && !defined(__NVCOMPILER)
    return {};
#endif
  }

  template <typename Mapping, size_t PaddingStrideIdx>
  MDSPAN_INLINE_FUNCTION static constexpr static_array_type
  init_padding([[maybe_unused]] const Mapping &other_mapping,
                      std::integral_constant<size_t, PaddingStrideIdx>) {
    if constexpr (Extents::rank() > 1) {
      return {other_mapping.stride(PaddingStrideIdx)};
    } else {
      return {};
    }
    // Missing return statement warning from NVCC and ICC
#if (defined(__NVCC__) || defined(__INTEL_COMPILER)) && !defined(__NVCOMPILER)
    return {};
#endif
  }
};
} // namespace detail

template <size_t PaddingValue>
template <class Extents>
class layout_left_padded<PaddingValue>::mapping {
public:
  static constexpr size_t padding_value = PaddingValue;

  using extents_type = Extents;
  using index_type = typename extents_type::index_type;
  using size_type = typename extents_type::size_type;
  using rank_type = typename extents_type::rank_type;
  using layout_type = layout_left_padded<padding_value>;

#ifndef MDSPAN_INTERNAL_TEST
private:
#endif // MDSPAN_INTERNAL_TEST

  static constexpr rank_type padded_stride_idx = detail::layout_padded_constants<layout_type, extents_type>::padded_stride_idx;
  static constexpr rank_type extent_to_pad_idx = detail::layout_padded_constants<layout_type, extents_type>::extent_to_pad_idx;

  static_assert((padding_value != 0)
                || (extents_type::static_extent(extent_to_pad_idx) == 0)
                || (extents_type::static_extent(extent_to_pad_idx) == dynamic_extent),
                "out of bounds access for rank 0");

  using padded_stride_type = detail::padded_extent< padding_value, extents_type, extent_to_pad_idx >;

  static constexpr size_t static_padding_stride = padded_stride_type::static_value();

  typename padded_stride_type::static_array_type padded_stride = {};
  extents_type exts = {};

  MDSPAN_INLINE_FUNCTION constexpr index_type
  compute_offset(std::index_sequence<>) const {
    return 0;
  }

  template <size_t Rank, class IndexOffset>
  MDSPAN_INLINE_FUNCTION constexpr index_type
  compute_offset(std::index_sequence<Rank>, IndexOffset index_offset) const {
    return index_offset;
  }

  template <size_t... Ranks, class... IndexOffsets>
  MDSPAN_INLINE_FUNCTION constexpr index_type
  compute_offset(std::index_sequence<Ranks...>,
                 IndexOffsets... index_offsets) const {
    index_type indices[] = {static_cast<index_type>(index_offsets)...};
    // self-recursive fold trick from
    // https://github.com/llvm/llvm-project/blob/96e1914aa2e6d8966acbfbe2f4d184201f1aa318/libcxx/include/mdspan/layout_left.h#L144
    index_type res = 0;
    ((res = indices[extents_type::rank() - 1 - Ranks] +
            ((extents_type::rank() - 1 - Ranks) == extent_to_pad_idx
                 ? padded_stride.value(0)
                 : exts.extent(extents_type::rank() - 1 - Ranks)) *
                res),
     ...);
    return res;
  }

public:
#if !MDSPAN_HAS_CXX_20 || defined(__NVCC__)
  MDSPAN_INLINE_FUNCTION
  constexpr mapping()
      : mapping(extents_type{})
  {}
#else
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr mapping()
    requires(static_padding_stride != dynamic_extent) = default;

  MDSPAN_INLINE_FUNCTION
  constexpr mapping()
    requires(static_padding_stride == dynamic_extent)
      : mapping(extents_type{})
  {}
#endif

  MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mapping(const mapping&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mapping& operator=(const mapping&) noexcept = default;

  /**
   * Initializes the mapping with the given extents.
   *
   * \param ext the given extents
   */
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const extents_type& ext)
    : padded_stride(padded_stride_type::init_padding(ext)), exts(ext)
  {}

  /**
   * Initializes the mapping with the given extents and the specified padding value.
   *
   * This overload participates in overload resolution only if `is_convertible_v<Size, index_type>`
   * is `true` and `is_nothrow_constructible_v<index_type, Size>` is `true`
   *
   * \param ext the given extents
   * \param padding_value the padding value
   */
  MDSPAN_TEMPLATE_REQUIRES(
    class Size,
    /* requires */ (
      std::is_convertible_v<Size, index_type>
      && std::is_nothrow_constructible_v<index_type, Size>
    )
  )
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const extents_type &ext, Size dynamic_padding_value)
      : padded_stride(padded_stride_type::init_padding(ext, dynamic_padding_value)), exts(ext)
  {
    assert((padding_value == dynamic_extent) || (static_cast<index_type>(padding_value) == static_cast<index_type>(dynamic_padding_value)));
  }

  /**
   * Converting constructor from `layout_left::mapping`.
   *
   * This overload participates in overload resolution only if
   * `is_constructible_v<extents_type, OtherExtents>` is true. If
   * `OtherExtents::rank() > 1` then one of `padding_value`, `static_extent(0)`,
   * or `OtherExtents::static_extent(0)` must be `dynamic_extent`; otherwise,
   * `OtherExtents::static_extent(0)` must be equal to the least multiple of
   * `padding_value` greater than or equal to `extents_type::static_extent(0)`
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (std::is_constructible_v<extents_type, OtherExtents>))
  MDSPAN_CONDITIONAL_EXPLICIT(
      (!std::is_convertible_v<OtherExtents, extents_type>))
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const layout_left::mapping<OtherExtents> &other_mapping)
      : padded_stride(padded_stride_type::init_padding(
            other_mapping,
            std::integral_constant<size_t, padded_stride_idx>{})),
        exts(other_mapping.extents()) {
    static_assert(
        (OtherExtents::rank() > 1) ||
        (static_padding_stride != dynamic_extent) ||
        (OtherExtents::static_extent(extent_to_pad_idx) != dynamic_extent) ||
        (static_padding_stride ==
         OtherExtents::static_extent(extent_to_pad_idx)));
  }

  /**
   * Converting constructor from `layout_stride::mapping`.
   *
   * This overload participates in overload resolution only if
   * `is_constructible_v<extents_type, OtherExtents>` is true
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (std::is_constructible_v<extents_type, OtherExtents>))
  MDSPAN_CONDITIONAL_EXPLICIT((extents_type::rank() > 0))
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const layout_stride::mapping<OtherExtents> &other_mapping)
      : padded_stride(padded_stride_type::init_padding(
            other_mapping,
            std::integral_constant<size_t, padded_stride_idx>{})),
        exts(other_mapping.extents()) {}

  /**
   * Converting constructor from `layout_left_padded::mapping`.
   *
   * This overload participates in overload resolution only if
   * `is_constructible_v<extents_type, OtherExtents>` is true. Either
   * `padding_value` or `OtherPaddingStride` must be `std::dynamic_extent`, or
   * `padding_value == OtherPaddingStride`.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class Mapping,
      /* requires */ (detail::is_layout_left_padded_mapping<Mapping>::value
                          &&std::is_constructible_v<
                              extents_type, typename Mapping::extents_type>))
  MDSPAN_CONDITIONAL_EXPLICIT((extents_type::rank() > 1 &&
                               (padding_value == dynamic_extent ||
                                Mapping::padding_value == dynamic_extent)))
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const Mapping &other_mapping)
      : padded_stride(padded_stride_type::init_padding(
            other_mapping,
            std::integral_constant<size_t, padded_stride_idx>{})),
        exts(other_mapping.extents()) {
    static_assert(padding_value == dynamic_extent ||
                  Mapping::padding_value == dynamic_extent ||
                  padding_value == Mapping::padding_value);
  }

  /**
   * Converting constructor from `layout_right_padded::mapping`.
   *
   * This overload participates in overload resolution only if
   * `extents_type::rank()` is 0 or 1 and `is_constructible_v<extents_type,
   * OtherExtents>` is `true`.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class Mapping,
      /* requires */ (detail::is_layout_right_padded_mapping<Mapping>::value
                              &&extents_type::rank() <= 1 &&
                      std::is_constructible_v<extents_type,
                                              typename Mapping::extents_type>))
  MDSPAN_CONDITIONAL_EXPLICIT(
      (!std::is_convertible_v<typename Mapping::extents_type, extents_type>))
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const Mapping &other_mapping) noexcept
      : padded_stride(padded_stride_type::init_padding(
            static_cast<extents_type>(other_mapping.extents()),
            other_mapping.extents().extent(extent_to_pad_idx))),
        exts(other_mapping.extents()) {}

  MDSPAN_INLINE_FUNCTION constexpr const extents_type &
  extents() const noexcept {
    return exts;
  }

  constexpr std::array<index_type, extents_type::rank()>
  strides() const noexcept {
    if constexpr (extents_type::rank() == 0) {
      return {};
    } else if constexpr (extents_type::rank() == 1) {
      return {1};
    } else {
      index_type value = 1;
      std::array<index_type, extents_type::rank()> s{};
      s[extent_to_pad_idx] = value;
      value *= padded_stride.value(0);
      for (rank_type r = extent_to_pad_idx + 1; r < extents_type::rank() - 1;
           ++r) {
        s[r] = value;
        value *= exts.extent(r);
      }
      s[extents_type::rank() - 1] = value;
      return s;
    }
  }

  MDSPAN_INLINE_FUNCTION constexpr index_type
  required_span_size() const noexcept {
    if constexpr (extents_type::rank() == 0) {
      return 1;
    } else if constexpr (extents_type::rank() == 1) {
      return exts.extent(0);
    } else {
      index_type value = padded_stride.value(0);
      for (rank_type r = 1; r < extents_type::rank(); ++r) {
        value *= exts.extent(r);
      }
      return value + exts.extent(0) - padded_stride.value(0);
    }
  }

  /**
   * Return the mapping given the provided indices per rank.
   *
   * This overload participates in overload resolution only if:
   * - `sizeof...(Indices) == extents_type::rank()`,
   * - `(is_convertible_v<Indices, index_type> && ...) is true`, and
   * - (is_nothrow_constructible_v<index_type, Indices> && ...) is true.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class... Indices,
      /* requires */ (sizeof...(Indices) == extents_type::rank() &&
                      (::MDSPAN_IMPL_STANDARD_NAMESPACE::detail::
                           are_valid_indices<index_type, Indices...>())))
  MDSPAN_INLINE_FUNCTION constexpr size_t
  operator()(Indices... idxs) const noexcept {
#if !defined(NDEBUG)
    ::MDSPAN_IMPL_STANDARD_NAMESPACE::detail::check_all_indices(this->extents(),
                                                                idxs...);
#endif // ! NDEBUG
    return compute_offset(std::index_sequence_for<Indices...>{}, idxs...);
  }

  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_unique() noexcept {
    return true;
  }
  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_exhaustive() noexcept {
    return (extents_type::rank() <= rank_type(1)) ||
           (extents_type::static_extent(extent_to_pad_idx) != dynamic_extent &&
            extents_type::static_extent(extent_to_pad_idx) ==
                padded_stride_type::static_value());
  }
  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_strided() noexcept {
    return true;
  }

  MDSPAN_INLINE_FUNCTION static constexpr bool is_unique() noexcept {
    return true;
  }
  MDSPAN_INLINE_FUNCTION constexpr bool is_exhaustive() const noexcept {
    return (extents_type::rank() < 2) ||
           (exts.extent(extent_to_pad_idx) == padded_stride.value(0));
  }
  MDSPAN_INLINE_FUNCTION static constexpr bool is_strided() noexcept {
    return true;
  }

  MDSPAN_INLINE_FUNCTION
  constexpr index_type stride(rank_type r) const noexcept {
    assert(r < extents_type::rank());
    if (r == 0)
      return index_type(1);

    index_type value = padded_stride.value(0);
    for (rank_type k = 1; k < r; k++)
      value *= exts.extent(k);

    return value;
  }

  /**
   * Equality operator between `layout_left_padded`s
   *
   * This overload only participates in overload resolution if
   * `OtherExtents::rank() == extents_type::rank()`.
   *
   * \note There is currently a difference from p2642r2, where this function is
   * specified as taking `layout_left_padded< padding_value >::mapping<
   * Extents>`. However, this makes `padding_value` non-deducible.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class Mapping,
      /* requires */ (detail::is_layout_left_padded_mapping<Mapping>::value &&
                      (Mapping::extents_type::rank() == extents_type::rank())))
  MDSPAN_INLINE_FUNCTION friend constexpr bool
  operator==(const mapping &left, const Mapping &right) noexcept {
    // Workaround for some compilers not short-circuiting properly with
    // compile-time checks i.e. we can't access stride(_padding_stride_idx) of a
    // rank 0 mapping
    bool strides_equal = true;
    if constexpr (extents_type::rank() > rank_type(1)) {
      strides_equal =
          left.stride(padded_stride_idx) == right.stride(padded_stride_idx);
    }
    return (left.extents() == right.extents()) && strides_equal;
  }

#if !MDSPAN_HAS_CXX_20
  /**
   * Inequality operator between `layout_left_padded`s
   *
   * This overload only participates in overload resolution if
   * `OtherExtents::rank() == extents_type::rank()`.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class Mapping,
      /* requires */ (detail::is_layout_left_padded_mapping<Mapping>::value &&
                      (Mapping::extents_type::rank() == extents_type::rank())))
  MDSPAN_INLINE_FUNCTION friend constexpr bool
  operator!=(const mapping &left, const Mapping &right) noexcept {
    return !(left == right);
  }
#endif

   // [mdspan.submdspan.mapping], submdspan mapping specialization
   template<class... SliceSpecifiers>
   MDSPAN_INLINE_FUNCTION
     constexpr auto submdspan_mapping_impl(
       SliceSpecifiers... slices) const;

   template<class... SliceSpecifiers>
   MDSPAN_INLINE_FUNCTION
     friend constexpr auto submdspan_mapping(
       const mapping& src, SliceSpecifiers... slices) {
         return src.submdspan_mapping_impl(slices...);
     }
};

template <size_t PaddingValue>
template <class Extents>
class layout_right_padded<PaddingValue>::mapping {
public:
  static constexpr size_t padding_value = PaddingValue;

  using extents_type = Extents;
  using index_type = typename extents_type::index_type;
  using size_type = typename extents_type::size_type;
  using rank_type = typename extents_type::rank_type;
  using layout_type = layout_right_padded<padding_value>;

#ifndef MDSPAN_INTERNAL_TEST
  private:
#endif // MDSPAN_INTERNAL_TEST

  static constexpr rank_type padded_stride_idx = detail::layout_padded_constants<layout_type, extents_type>::padded_stride_idx;
  static constexpr rank_type extent_to_pad_idx = detail::layout_padded_constants<layout_type, extents_type>::extent_to_pad_idx;

  static_assert((padding_value != 0)
                || (extents_type::static_extent(extent_to_pad_idx) == 0)
                || (extents_type::static_extent(extent_to_pad_idx) == dynamic_extent),
                "if padding stride is 0, static_extent(extent-to-pad-rank) must also be 0 or dynamic_extent");

  using padded_stride_type = detail::padded_extent< padding_value, extents_type, extent_to_pad_idx >;
  static constexpr size_t static_padding_stride = padded_stride_type::static_value();

  typename padded_stride_type::static_array_type padded_stride = {};
  extents_type exts = {};

  MDSPAN_INLINE_FUNCTION constexpr index_type
  compute_offset(std::index_sequence<>) const {
    return 0;
  }

  template <size_t Rank, class IndexOffset>
  MDSPAN_INLINE_FUNCTION constexpr index_type
  compute_offset(std::index_sequence<Rank>, IndexOffset index_offset) const {
    return index_offset;
  }

  template <size_t... Ranks, class... IndexOffsets>
  MDSPAN_INLINE_FUNCTION constexpr index_type
  compute_offset(std::index_sequence<Ranks...>,
                 IndexOffsets... index_offsets) const {
    // self-recursive fold trick from
    // https://github.com/llvm/llvm-project/blob/4d9771741d40cc9cfcccb6b033f43689d36b705a/libcxx/include/mdspan/layout_right.h#L141
    index_type res = 0;
    ((res = static_cast<index_type>(index_offsets) +
            (Ranks == extent_to_pad_idx ? padded_stride.value(0)
                                        : exts.extent(Ranks)) *
                res),
     ...);
    return res;
  }

public:
#if !MDSPAN_HAS_CXX_20 || defined(__NVCC__)
  MDSPAN_INLINE_FUNCTION
      constexpr mapping()
      : mapping(extents_type{})
  {}
#else
  MDSPAN_INLINE_FUNCTION_DEFAULTED
      constexpr mapping()
    requires(static_padding_stride != dynamic_extent) = default;

  MDSPAN_INLINE_FUNCTION
      constexpr mapping()
    requires(static_padding_stride == dynamic_extent)
      : mapping(extents_type{})
  {}
#endif

  MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mapping(const mapping&) noexcept = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mapping& operator=(const mapping&) noexcept = default;

  /**
   * Initializes the mapping with the given extents.
   *
   * \param ext the given extents
   */
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const extents_type &ext)
      : padded_stride(padded_stride_type::init_padding(ext)), exts(ext) {}

  /**
   * Initializes the mapping with the given extents and the specified padding value.
   *
   * This overload participates in overload resolution only if `is_convertible_v<Size, index_type>`
   * is `true` and `is_nothrow_constructible_v<index_type, Size>` is `true`
   *
   * \param ext the given extents
   * \param padding_value the padding value
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class Size,
      /* requires */ (
          std::is_convertible_v<Size, index_type>
              && std::is_nothrow_constructible_v<index_type, Size>
          )
      )
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const extents_type &ext, Size dynamic_padding_value)
      : padded_stride(padded_stride_type::init_padding(ext, static_cast<index_type>(dynamic_padding_value))),
        exts(ext) {
    assert((padding_value == dynamic_extent) ||
           (static_cast<index_type>(padding_value) == static_cast<index_type>(dynamic_padding_value)));
  }

  /**
   * Converting constructor from `layout_right::mapping`.
   *
   * This overload participates in overload resolution only if `is_constructible_v<extents_type, OtherExtents>` is true.
   * If `OtherExtents::rank() > 1` then one of `padding_value`, `static_extent(0)`, or `OtherExtents::static_extent(0)` must be `dynamic_extent`;
   * otherwise, `OtherExtents::static_extent(0)` must be equal to the least multiple of `padding_value` greater than or equal to `extents_type::static_extent(0)`
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (std::is_constructible_v<extents_type, OtherExtents>))
  MDSPAN_CONDITIONAL_EXPLICIT(
      (!std::is_convertible_v<OtherExtents, extents_type>))
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const layout_right::mapping<OtherExtents> &other_mapping)
      : padded_stride(padded_stride_type::init_padding(
            other_mapping,
            std::integral_constant<size_t, padded_stride_idx>{})),
        exts(other_mapping.extents()) {
    static_assert(
        (OtherExtents::rank() > 1) ||
        (padded_stride_type::static_value() != dynamic_extent) ||
        (OtherExtents::static_extent(extent_to_pad_idx) != dynamic_extent) ||
        (padded_stride_type::static_value() ==
         OtherExtents::static_extent(extent_to_pad_idx)));
  }

  /**
   * Converting constructor from `layout_stride::mapping`.
   *
   * This overload participates in overload resolution only if
   * `is_constructible_v<extents_type, OtherExtents>` is true
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (std::is_constructible_v<extents_type, OtherExtents>))
  MDSPAN_CONDITIONAL_EXPLICIT((extents_type::rank() > 0))
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const layout_stride::mapping<OtherExtents> &other_mapping)
      : padded_stride(padded_stride_type::init_padding(
            other_mapping,
            std::integral_constant<size_t, padded_stride_idx>{})),
        exts(other_mapping.extents()) {}

  /**
   * Converting constructor from `layout_right_padded::mapping`.
   *
   * This overload participates in overload resolution only if
   * `is_constructible_v<extents_type, OtherExtents>` is true. Either
   * `padding_value` or `OtherPaddingStride` must be `std::dynamic_extent`, or
   * `padding_value == OtherPaddingStride`.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class Mapping,
      /* requires */ (detail::is_layout_right_padded_mapping<Mapping>::value
                          &&std::is_constructible_v<
                              extents_type, typename Mapping::extents_type>))
  MDSPAN_CONDITIONAL_EXPLICIT((extents_type::rank() > 1 &&
                               (padding_value == dynamic_extent ||
                                Mapping::padding_value == dynamic_extent)))
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const Mapping &other_mapping)
      : padded_stride(padded_stride_type::init_padding(
            other_mapping,
            std::integral_constant<size_t, padded_stride_idx>{})),
        exts(other_mapping.extents()) {
    static_assert(padding_value == dynamic_extent ||
                  Mapping::padding_value == dynamic_extent ||
                  padding_value == Mapping::padding_value);
  }

  /**
   * Converting constructor from `layout_left_padded::mapping`.
   *
   * This overload participates in overload resolution only if
   * `extents_type::rank()` is 0 or 1 and `is_constructible_v<extents_type,
   * OtherExtents>` is `true`.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class Mapping,
      /* requires */ (detail::is_layout_left_padded_mapping<Mapping>::value
                              &&extents_type::rank() <= 1 &&
                      std::is_constructible_v<extents_type,
                                              typename Mapping::extents_type>))
  MDSPAN_CONDITIONAL_EXPLICIT(
      (!std::is_convertible_v<typename Mapping::extents_type, extents_type>))
  MDSPAN_INLINE_FUNCTION
  constexpr mapping(const Mapping &other_mapping) noexcept
      : padded_stride(padded_stride_type::init_padding(
            static_cast<extents_type>(other_mapping.extents()),
            other_mapping.extents().extent(extent_to_pad_idx))),
        exts(other_mapping.extents()) {}

  MDSPAN_INLINE_FUNCTION constexpr const extents_type &
  extents() const noexcept {
    return exts;
  }

  constexpr std::array<index_type, extents_type::rank()>
  strides() const noexcept {
    if constexpr (extents_type::rank() == 0) {
      return {};
    } else if constexpr (extents_type::rank() == 1) {
      return {1};
    } else {
      index_type value = 1;
      std::array<index_type, extents_type::rank()> s{};
      s[extent_to_pad_idx] = value;
      value *= padded_stride.value(0);
      for (rank_type r = extent_to_pad_idx - 1; r > 0; --r) {
        s[r] = value;
        value *= exts.extent(r);
      }
      s[0] = value;
      return s;
    }
  }

  MDSPAN_INLINE_FUNCTION constexpr index_type
  required_span_size() const noexcept {
    if constexpr (extents_type::rank() == 0) {
      return 1;
    } else if constexpr (extents_type::rank() == 1) {
      return exts.extent(0);
    } else {
      index_type value = padded_stride.value(0);
      for (rank_type r = 0; r < extent_to_pad_idx; ++r) {
        value *= exts.extent(r);
      }
      return value + exts.extent(extent_to_pad_idx) - padded_stride.value(0);
    }
  }

  /**
   * Return the mapping given the provided indices per rank.
   *
   * This overload participates in overload resolution only if:
   * - `sizeof...(Indices) == extents_type::rank()`,
   * - `(is_convertible_v<Indices, index_type> && ...) is true`, and
   * - (is_nothrow_constructible_v<index_type, Indices> && ...) is true.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class... Indices,
      /* requires */ (sizeof...(Indices) == extents_type::rank() &&
                      (::MDSPAN_IMPL_STANDARD_NAMESPACE::detail::
                           are_valid_indices<index_type, Indices...>())))
  MDSPAN_INLINE_FUNCTION constexpr size_t
  operator()(Indices... idxs) const noexcept {
    return compute_offset(std::index_sequence_for<Indices...>{}, idxs...);
  }

  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_unique() noexcept {
    return true;
  }
  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_exhaustive() noexcept {
    return (extents_type::rank() <= rank_type(1)) ||
           (extents_type::static_extent(extent_to_pad_idx) != dynamic_extent &&
            extents_type::static_extent(extent_to_pad_idx) ==
                padded_stride_type::static_value());
  }
  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_strided() noexcept {
    return true;
  }

  MDSPAN_INLINE_FUNCTION static constexpr bool is_unique() noexcept {
    return true;
  }
  MDSPAN_INLINE_FUNCTION constexpr bool is_exhaustive() const noexcept {
    return (extents_type::rank() < 2) ||
           (exts.extent(extent_to_pad_idx) == padded_stride.value(0));
  }
  MDSPAN_INLINE_FUNCTION static constexpr bool is_strided() noexcept {
    return true;
  }

  MDSPAN_INLINE_FUNCTION constexpr index_type
  stride(rank_type r) const noexcept {
    assert(r < extents_type::rank());
    if (r == extents_type::rank() - 1)
      return index_type(1);

    index_type value = padded_stride.value(0);
    for (rank_type k = extents_type::rank() - 2; k > r; k--)
      value *= exts.extent(k);

    return value;
  }

  /**
   * Equality operator between `layout_right_padded`s
   *
   * This overload only participates in overload resolution if
   * `OtherExtents::rank() == extents_type::rank()`.
   *
   * \note There is currently a difference from p2642r2, where this function is
   * specified as taking `layout_right_padded< padding_value >::mapping<
   * Extents>`. However, this makes `padding_value` non-deducible.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class Mapping,
      /* requires */ (detail::is_layout_right_padded_mapping<Mapping>::value &&
                      (Mapping::extents_type::rank() == extents_type::rank())))
  MDSPAN_INLINE_FUNCTION friend constexpr bool
  operator==(const mapping &left, const Mapping &right) noexcept {
    // Workaround for some compilers not short-circuiting properly with
    // compile-time checks i.e. we can't access stride(_padding_stride_idx) of a
    // rank 0 mapping
    bool strides_equal = true;
    if constexpr (extents_type::rank() > rank_type(1)) {
      strides_equal =
          left.stride(padded_stride_idx) == right.stride(padded_stride_idx);
    }
    return (left.extents() == right.extents()) && strides_equal;
  }

#if !MDSPAN_HAS_CXX_20
  /**
   * Inequality operator between `layout_right_padded`s
   *
   * This overload only participates in overload resolution if
   * `OtherExtents::rank() == extents_type::rank()`.
   */
  MDSPAN_TEMPLATE_REQUIRES(
      class Mapping,
      /* requires */ (detail::is_layout_right_padded_mapping<Mapping>::value &&
                      (Mapping::extents_type::rank() == extents_type::rank())))
  MDSPAN_INLINE_FUNCTION friend constexpr bool
  operator!=(const mapping &left, const Mapping &right) noexcept {
    return !(left == right);
  }
#endif

   // [mdspan.submdspan.mapping], submdspan mapping specialization
   template<class... SliceSpecifiers>
   MDSPAN_INLINE_FUNCTION
     constexpr auto submdspan_mapping_impl(
       SliceSpecifiers... slices) const;

   template<class... SliceSpecifiers>
   MDSPAN_INLINE_FUNCTION
     friend constexpr auto submdspan_mapping(
       const mapping& src, SliceSpecifiers... slices) {
         return src.submdspan_mapping_impl(slices...);
     }
};
}
}
