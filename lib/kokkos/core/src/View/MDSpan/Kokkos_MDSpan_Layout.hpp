// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_EXPERIMENTAL_MDSPAN_LAYOUT_HPP
#define KOKKOS_EXPERIMENTAL_MDSPAN_LAYOUT_HPP

#include "Kokkos_MDSpan_Extents.hpp"
#include <View/Kokkos_ViewDataAnalysis.hpp>

#ifdef KOKKOS_ENABLE_IMPL_CHECK_POSSIBLY_BREAKING_LAYOUTS
#include <iostream>
#endif

// The difference between a legacy Kokkos array layout and an
// mdspan layout is that the array layouts can have state, but don't have the
// nested mapping. This file provides interoperability helpers.

namespace Kokkos::Impl {
// We do have implementation detail versions of these in our mdspan impl
// However they are not part of the public standard interface
template <class>
struct IsLayoutRightPadded : std::false_type {};

template <size_t Pad>
struct IsLayoutRightPadded<::Kokkos::Experimental::layout_right_padded<Pad>>
    : std::true_type {};

template <class>
struct IsLayoutLeftPadded : std::false_type {};

template <size_t Pad>
struct IsLayoutLeftPadded<::Kokkos::Experimental::layout_left_padded<Pad>>
    : std::true_type {};

template <class ArrayLayout>
struct LayoutFromArrayLayout {
  using type = void;
};

template <>
struct LayoutFromArrayLayout<LayoutLeft> {
  using type = ::Kokkos::Experimental::layout_left_padded<dynamic_extent>;
};

template <>
struct LayoutFromArrayLayout<LayoutRight> {
  using type = ::Kokkos::Experimental::layout_right_padded<dynamic_extent>;
};

template <>
struct LayoutFromArrayLayout<LayoutStride> {
  using type = layout_stride;
};

template <class ArrayLayout, class MDSpanType>
KOKKOS_INLINE_FUNCTION auto array_layout_from_mapping(
    const typename MDSpanType::mapping_type &mapping) {
  using mapping_type = typename MDSpanType::mapping_type;
  using extents_type = typename mapping_type::extents_type;

  constexpr auto rank = extents_type::rank();
  const auto &ext     = mapping.extents();

  static_assert(rank <= ARRAY_LAYOUT_MAX_RANK,
                "Unsupported rank for mdspan (must be <= 8)");

  if constexpr (std::is_same_v<ArrayLayout, LayoutStride>) {
    return Kokkos::LayoutStride{
        rank > 0 ? ext.extent(0) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        rank > 0 ? mapping.stride(0) : 0,
        rank > 1 ? ext.extent(1) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        rank > 1 ? mapping.stride(1) : 0,
        rank > 2 ? ext.extent(2) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        rank > 2 ? mapping.stride(2) : 0,
        rank > 3 ? ext.extent(3) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        rank > 3 ? mapping.stride(3) : 0,
        rank > 4 ? ext.extent(4) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        rank > 4 ? mapping.stride(4) : 0,
        rank > 5 ? ext.extent(5) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        rank > 5 ? mapping.stride(5) : 0,
        rank > 6 ? ext.extent(6) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        rank > 6 ? mapping.stride(6) : 0,
        rank > 7 ? ext.extent(7) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
        rank > 7 ? mapping.stride(7) : 0,
    };
  } else {
    ArrayLayout layout{rank > 0 ? ext.extent(0) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                       rank > 1 ? ext.extent(1) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                       rank > 2 ? ext.extent(2) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                       rank > 3 ? ext.extent(3) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                       rank > 4 ? ext.extent(4) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                       rank > 5 ? ext.extent(5) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                       rank > 6 ? ext.extent(6) : KOKKOS_IMPL_CTOR_DEFAULT_ARG,
                       rank > 7 ? ext.extent(7) : KOKKOS_IMPL_CTOR_DEFAULT_ARG};

    if constexpr (rank > 1 &&
                  std::is_same_v<typename mapping_type::layout_type,
                                 Kokkos::Experimental::layout_left_padded<
                                     dynamic_extent>>) {
      layout.stride = mapping.stride(1);
    }
    if constexpr (std::is_same_v<typename mapping_type::layout_type,
                                 Kokkos::Experimental::layout_right_padded<
                                     dynamic_extent>>) {
// Legacy LayoutRight actually padded the left most dimension, not like
// layout_right_padded the right most dimension. Thus if the stride wasn't
// matching the appropriate extent this conversion doesn't work.
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
      if constexpr (rank == 2) {
        layout.stride = mapping.stride(0);
      }
      if constexpr (rank > 2) {
        if (mapping.stride(rank - 2) != mapping.extents().extent(rank - 1))
          Kokkos::abort(
              "Invalid conversion from layout_right_padded to LayoutRight");
      }
#else
      if constexpr (rank > 1) {
        layout.stride = mapping.stride(rank - 2);
      } else {
        // Just setting the stride to 1 for rank 0/1
        layout.stride = 1;
      }
#endif
    }
    return layout;
  }
}

template <class MappingType, class ArrayLayout, size_t... Idx>
KOKKOS_INLINE_FUNCTION auto mapping_from_array_layout_impl(
    ArrayLayout layout, std::index_sequence<Idx...>) {
  using index_type   = typename MappingType::index_type;
  using extents_type = typename MappingType::extents_type;
  if constexpr (std::is_same_v<typename MappingType::layout_type,
                               layout_left> ||
                std::is_same_v<typename MappingType::layout_type,
                               layout_right>) {
    return MappingType{
        extents_type{dextents<index_type, MappingType::extents_type::rank()>{
            layout.dimension[Idx]...}}};
  } else {
    if (layout.stride == KOKKOS_IMPL_CTOR_DEFAULT_ARG ||
        extents_type::rank() < 2) {
      return MappingType{
          extents_type{dextents<index_type, MappingType::extents_type::rank()>{
              layout.dimension[Idx]...}}};
    } else {
// Handle DEFAULT_ARG, should be layout_dimension 0 or n -1
// assert that this is not default_arg, as a tool for people to
// transition their code and avoid breaking changes
#ifdef KOKKOS_ENABLE_IMPL_CHECK_POSSIBLY_BREAKING_LAYOUTS
      KOKKOS_IF_ON_HOST(
          (if constexpr (std::is_same_v<ArrayLayout, LayoutRight> &&
                         extents_type::rank() > 2) {
            if (layout.stride != KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
              std::cerr
                  << "The layout of values in this Kokkos View may be "
                     "different due "
                     "to a non-defaulted stride. Verify that this is not an "
                     "issue for "
                     "your Views and then disable "
                     "KOKKOS_ENABLE_IMPL_CHECK_POSSIBLY_BREAKING_LAYOUTS.\n";
            }
          }))
#endif

      if (layout.stride == KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
        return MappingType{extents_type{
            dextents<index_type, MappingType::extents_type::rank()>{
                layout.dimension[Idx]...}}};
      } else {
        return MappingType{
            extents_type{
                dextents<index_type, MappingType::extents_type::rank()>{
                    layout.dimension[Idx]...}},
            layout.stride};
      }
    }
  }
}

template <class MappingType, size_t... Idx>
KOKKOS_INLINE_FUNCTION auto mapping_from_array_layout_impl(
    LayoutStride layout, std::index_sequence<Idx...>) {
  static_assert(
      std::is_same_v<typename MappingType::layout_type, layout_stride>);
  using index_type = typename MappingType::index_type;
  index_type strides[MappingType::extents_type::rank()] = {
      layout.stride[Idx]...};
  return MappingType{
      mdspan_non_standard_tag(),
      static_cast<typename MappingType::extents_type>(
          dextents<index_type, MappingType::extents_type::rank()>{
              layout.dimension[Idx]...}),
      strides};
}

// specialization for rank 0 to avoid empty array
template <class MappingType>
KOKKOS_INLINE_FUNCTION auto mapping_from_array_layout_impl(
    LayoutStride, std::index_sequence<>) {
  return MappingType{};
}

template <class MappingType, class ArrayLayout>
KOKKOS_INLINE_FUNCTION auto mapping_from_array_layout(ArrayLayout layout) {
  return mapping_from_array_layout_impl<MappingType>(
      layout, std::make_index_sequence<MappingType::extents_type::rank()>());
}

template <class MDSpanType, class VM>
KOKKOS_INLINE_FUNCTION auto mapping_from_view_mapping(const VM &view_mapping) {
  using mapping_type = typename MDSpanType::mapping_type;
  using extents_type = typename mapping_type::extents_type;

  // std::span is not available in C++17 (our current requirements),
  // so we need to use the std::array constructor for layout mappings.
  // FIXME When C++20 is available, we can use std::span here instead
  std::size_t strides[VM::Rank == 0 ? 1 : VM::Rank];
  view_mapping.stride_fill(&strides[0]);
  if constexpr (std::is_same_v<typename mapping_type::layout_type,
                               Kokkos::layout_stride>) {
    return mapping_type(Kokkos::mdspan_non_standard,
                        extents_from_view_mapping<extents_type>(view_mapping),
                        strides);
  } else if constexpr (VM::Rank > 1 &&
                       std::is_same_v<typename mapping_type::layout_type,
                                      Kokkos::Experimental::layout_left_padded<
                                          Kokkos::dynamic_extent>>) {
    return mapping_type(extents_from_view_mapping<extents_type>(view_mapping),
                        strides[1]);
  } else if constexpr (VM::Rank > 1 &&
                       std::is_same_v<typename mapping_type::layout_type,
                                      Kokkos::Experimental::layout_right_padded<
                                          Kokkos::dynamic_extent>>) {
    return mapping_type(extents_from_view_mapping<extents_type>(view_mapping),
                        strides[VM::Rank - 2]);
  } else {
    return mapping_type(extents_from_view_mapping<extents_type>(view_mapping));
  }
}

template <size_t ScalarSize>
struct Padding {
  static constexpr size_t div =
      ScalarSize == 0 ? 0 : static_cast<size_t>(MEMORY_ALIGNMENT) / ScalarSize;
  static constexpr size_t mod =
      ScalarSize == 0 ? 0 : static_cast<size_t>(MEMORY_ALIGNMENT) % ScalarSize;

  // If memory alignment is a multiple of the trivial scalar size then attempt
  // to align.
  static constexpr size_t align  = ScalarSize != 0 && mod == 0 ? div : 0;
  static constexpr size_t div_ok = (div != 0) ? div : 1;

  KOKKOS_INLINE_FUNCTION
  static constexpr size_t stride(size_t const N) {
    return ((align != 0) &&
            ((static_cast<size_t>(MEMORY_ALIGNMENT_THRESHOLD) * align) < N) &&
            ((N % div_ok) != 0))
               ? N + align - (N % div_ok)
               : N;
  }
};

template <class MappingType, size_t ScalarSize, class ViewCtorProperties,
          class... Sizes>
KOKKOS_INLINE_FUNCTION auto mapping_from_ctor_and_sizes(
    const ViewCtorProperties &, const Sizes... args) {
  using layout_t = typename MappingType::layout_type;
  using ext_t    = typename MappingType::extents_type;
  ext_t ext{args...};
  constexpr bool padded = ViewCtorProperties::allow_padding;
  if constexpr (IsLayoutLeftPadded<layout_t>::value && padded &&
                ext_t::rank() > 1) {
    return MappingType(ext, Padding<ScalarSize>::stride(ext.extent(0)));
  } else if constexpr (IsLayoutRightPadded<layout_t>::value && padded &&
                       ext_t::rank() > 1) {
    return MappingType(
        ext, Padding<ScalarSize>::stride(ext.extent(ext_t::rank() - 1)));
  } else {
    return MappingType(ext);
  }
}

template <class MappingType, size_t ScalarSize, class ViewCtorProperties>
KOKKOS_INLINE_FUNCTION auto mapping_from_ctor_and_8sizes(
    const ViewCtorProperties &arg_prop, [[maybe_unused]] const size_t arg_N0,
    [[maybe_unused]] const size_t arg_N1, [[maybe_unused]] const size_t arg_N2,
    [[maybe_unused]] const size_t arg_N3, [[maybe_unused]] const size_t arg_N4,
    [[maybe_unused]] const size_t arg_N5, [[maybe_unused]] const size_t arg_N6,
    [[maybe_unused]] const size_t arg_N7) {
  if constexpr (MappingType::extents_type::rank() == 0) {
    return mapping_from_ctor_and_sizes<MappingType, ScalarSize>(arg_prop);
  } else if constexpr (MappingType::extents_type::rank() == 1) {
    return mapping_from_ctor_and_sizes<MappingType, ScalarSize>(arg_prop,
                                                                arg_N0);
  } else if constexpr (MappingType::extents_type::rank() == 2) {
    return mapping_from_ctor_and_sizes<MappingType, ScalarSize>(arg_prop,
                                                                arg_N0, arg_N1);
  } else if constexpr (MappingType::extents_type::rank() == 3) {
    return mapping_from_ctor_and_sizes<MappingType, ScalarSize>(
        arg_prop, arg_N0, arg_N1, arg_N2);
  } else if constexpr (MappingType::extents_type::rank() == 4) {
    return mapping_from_ctor_and_sizes<MappingType, ScalarSize>(
        arg_prop, arg_N0, arg_N1, arg_N2, arg_N3);
  } else if constexpr (MappingType::extents_type::rank() == 5) {
    return mapping_from_ctor_and_sizes<MappingType, ScalarSize>(
        arg_prop, arg_N0, arg_N1, arg_N2, arg_N3, arg_N4);
  } else if constexpr (MappingType::extents_type::rank() == 6) {
    return mapping_from_ctor_and_sizes<MappingType, ScalarSize>(
        arg_prop, arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5);
  } else if constexpr (MappingType::extents_type::rank() == 7) {
    return mapping_from_ctor_and_sizes<MappingType, ScalarSize>(
        arg_prop, arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5, arg_N6);
  } else if constexpr (MappingType::extents_type::rank() == 8) {
    return mapping_from_ctor_and_sizes<MappingType, ScalarSize>(
        arg_prop, arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5, arg_N6,
        arg_N7);
  }
}
}  // namespace Kokkos::Impl

#endif  // KOKKOS_EXPERIMENTAL_MDSPAN_LAYOUT_HPP
