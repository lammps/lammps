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

#include <gtest/gtest.h>
#include <type_traits>

#include <Kokkos_Core.hpp>
#include "experimental/__p0009_bits/layout_stride.hpp"

namespace {

template <class T, class ExecutionSpace>
struct TestViewMDSpanConversion {
  using value_type = T;

  template <std::size_t Padding>
  using layout_left_padded = Kokkos::Experimental::layout_left_padded<Padding>;

  template <std::size_t Padding>
  using layout_right_padded =
      Kokkos::Experimental::layout_right_padded<Padding>;

  struct TestAccessor {
    using offset_policy    = TestAccessor;
    using element_type     = value_type;
    using reference        = element_type &;
    using data_handle_type = element_type *;

    constexpr TestAccessor() noexcept = default;
    constexpr reference access(data_handle_type p, std::size_t i) noexcept {
      return p[i];
    }
    constexpr data_handle_type offset(data_handle_type p,
                                      std::size_t i) noexcept {
      return p + i;
    }
  };

  template <class KokkosLayout, class DataType, class MDSpanLayoutMapping,
            class... RefViewProps>
  static void test_conversion_from_mdspan(
      Kokkos::View<DataType, RefViewProps...> ref,
      const MDSpanLayoutMapping &mapping) {
    using unmanaged_view_type =
        Kokkos::View<DataType, KokkosLayout, ExecutionSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
    using natural_mdspan_type = typename Kokkos::Impl::MDSpanViewTraits<
        typename unmanaged_view_type::traits>::mdspan_type;
    using mapping_type       = MDSpanLayoutMapping;
    using mdspan_layout_type = typename MDSpanLayoutMapping::layout_type;
    using extents_type       = typename mapping_type::extents_type;
    using mdspan_type =
        Kokkos::mdspan<value_type, extents_type, mdspan_layout_type>;

    static_assert(std::is_constructible_v<natural_mdspan_type, mdspan_type>);
    static_assert(std::is_convertible_v<mdspan_type, natural_mdspan_type> ==
                  std::is_convertible_v<mdspan_type, unmanaged_view_type>);
    // Manually create an mdspan from ref so we have a valid pointer to play
    // with
    const auto &exts = mapping.extents();
    auto mds         = mdspan_type{ref.data(), mapping};

    auto test_view = unmanaged_view_type(mds);

    ASSERT_EQ(test_view.data(), ref.data());
    ASSERT_EQ(test_view.data(), mds.data_handle());
    ASSERT_EQ(test_view.layout(), ref.layout());
    for (std::size_t r = 0; r < mdspan_type::rank(); ++r) {
      ASSERT_EQ(test_view.extent(r), ref.extent(r));
      ASSERT_EQ(test_view.extent(r), exts.extent(r));
    }
  }

  template <class MDSpanLayoutMapping, class ViewType>
  static void test_conversion_to_mdspan(
      const MDSpanLayoutMapping &ref_layout_mapping, ViewType v) {
    using view_type           = ViewType;
    using natural_mdspan_type = typename Kokkos::Impl::MDSpanViewTraits<
        typename view_type::traits>::mdspan_type;

    static_assert(natural_mdspan_type::rank() == view_type::rank);
    static_assert(std::is_same_v<typename natural_mdspan_type::value_type,
                                 typename view_type::value_type>);
    constexpr bool is_strided_layout =
        std::is_same_v<typename MDSpanLayoutMapping::layout_type,
                       Kokkos::layout_stride>;
    if constexpr (!is_strided_layout) {
      static_assert(natural_mdspan_type::mapping_type::padding_value ==
                    Kokkos::dynamic_extent);
    }
    // test conversion operator to natural mdspan
    {
      natural_mdspan_type cvt = v;
      ASSERT_EQ(cvt.data_handle(), v.data());
      ASSERT_EQ(cvt.mapping(), ref_layout_mapping);

      if constexpr (!is_strided_layout && natural_mdspan_type::rank() > 1) {
        ASSERT_EQ(cvt.mapping().stride(1), ref_layout_mapping.stride(1));
      }
    }
    // test to_mdspan() returning natural mdspan
    {
      auto cvt = v.to_mdspan();
      static_assert(std::is_same_v<natural_mdspan_type, decltype(cvt)>);
      ASSERT_EQ(cvt.data_handle(), v.data());
      ASSERT_EQ(cvt.mapping(), ref_layout_mapping);
    }
    // test conversion operator to different mdspan type
    {
      using element_type   = const typename natural_mdspan_type::element_type;
      using const_acc_type = Kokkos::Impl::SpaceAwareAccessor<
          typename ViewType::memory_space,
          Kokkos::default_accessor<element_type>>;
      using mdspan_type = Kokkos::mdspan<
          element_type,
          Kokkos::dextents<typename natural_mdspan_type::index_type,
                           natural_mdspan_type::rank()>,
          typename natural_mdspan_type::layout_type, const_acc_type>;
      mdspan_type cvt = v;
      ASSERT_EQ(cvt.data_handle(), v.data());
      ASSERT_EQ(cvt.mapping(), ref_layout_mapping);
    }
  }

  template <class MDSpanLayoutMapping, class ViewType, class AccessorType>
  static void test_conversion_to_mdspan_with_accessor(
      const MDSpanLayoutMapping &ref_layout_mapping, ViewType v,
      const AccessorType &a) {
    auto cvt = v.to_mdspan(a);
    static_assert(decltype(cvt)::rank() == ViewType::rank);
    static_assert(std::is_same_v<typename decltype(cvt)::value_type,
                                 typename ViewType::value_type>);
    ASSERT_EQ(cvt.data_handle(), v.data());
    ASSERT_EQ(cvt.mapping(), ref_layout_mapping);
  }

  template <typename ViewType>
  using natural_mdspan_type_for_view = typename Kokkos::Impl::MDSpanViewTraits<
      typename ViewType::traits>::mdspan_type;

  static void run_test() {
    // Verify we can only convert to compatible mdspans
    static_assert(std::is_convertible_v<
                  Kokkos::View<value_type *>,
                  natural_mdspan_type_for_view<Kokkos::View<value_type *>>>);
    static_assert(
        std::is_convertible_v<
            Kokkos::View<value_type *>,
            natural_mdspan_type_for_view<Kokkos::View<const value_type *>>>);

    // Do not cast const away
    static_assert(!std::is_convertible_v<
                  Kokkos::View<const value_type *>,
                  natural_mdspan_type_for_view<Kokkos::View<value_type *>>>);

    // Mismatched dim
    static_assert(!std::is_convertible_v<
                  Kokkos::View<value_type *>,
                  natural_mdspan_type_for_view<Kokkos::View<value_type **>>>);

    // Mismatched layouts
    static_assert(
        !std::is_convertible_v<Kokkos::View<value_type **, Kokkos::LayoutLeft>,
                               natural_mdspan_type_for_view<Kokkos::View<
                                   value_type **, Kokkos::LayoutRight>>>);
    static_assert(
        !std::is_convertible_v<Kokkos::View<value_type **, Kokkos::LayoutRight>,
                               natural_mdspan_type_for_view<Kokkos::View<
                                   value_type **, Kokkos::LayoutLeft>>>);
    // nvcc doesn't do CTAD properly here, making this way more verbose..
    // LayoutLeft
    test_conversion_from_mdspan<Kokkos::LayoutLeft>(
        Kokkos::View<value_type *, Kokkos::LayoutLeft, ExecutionSpace>("ref",
                                                                       7),
        typename layout_left_padded<7>::template mapping<
            Kokkos::dextents<std::size_t, 1>>{
            Kokkos::dextents<std::size_t, 1>(7)});

    test_conversion_from_mdspan<Kokkos::LayoutLeft>(
        Kokkos::View<value_type[7], Kokkos::LayoutLeft, ExecutionSpace>("ref"),
        typename layout_left_padded<7>::template mapping<
            Kokkos::extents<std::size_t, 7>>{
            Kokkos::extents<std::size_t, 7>()});
    test_conversion_from_mdspan<Kokkos::LayoutLeft>(
        Kokkos::View<value_type[7], Kokkos::LayoutLeft, ExecutionSpace>("ref"),
        typename layout_left_padded<7>::template mapping<
            Kokkos::dextents<std::size_t, 1>>{
            Kokkos::dextents<std::size_t, 1>(7)});
    test_conversion_from_mdspan<Kokkos::LayoutLeft>(
        Kokkos::View<value_type *, Kokkos::LayoutLeft, ExecutionSpace>("ref",
                                                                       7),
        typename layout_left_padded<7>::template mapping<
            Kokkos::extents<std::size_t, 7>>{
            Kokkos::extents<std::size_t, 7>()});

    test_conversion_from_mdspan<Kokkos::LayoutLeft>(
        Kokkos::View<value_type **, Kokkos::LayoutLeft, ExecutionSpace>("ref",
                                                                        7, 3),
        typename layout_left_padded<7>::template mapping<
            Kokkos::dextents<std::size_t, 2>>{
            Kokkos::dextents<std::size_t, 2>(7, 3)});
    test_conversion_from_mdspan<Kokkos::LayoutLeft>(
        Kokkos::View<value_type[7][3], Kokkos::LayoutLeft, ExecutionSpace>(
            "ref"),
        typename layout_left_padded<7>::template mapping<
            Kokkos::extents<std::size_t, 7, 3>>{
            Kokkos::extents<std::size_t, 7, 3>()});
    test_conversion_from_mdspan<Kokkos::LayoutLeft>(
        Kokkos::View<value_type[7][3], Kokkos::LayoutLeft, ExecutionSpace>(
            "ref"),
        typename layout_left_padded<7>::template mapping<
            Kokkos::dextents<std::size_t, 2>>{
            Kokkos::dextents<std::size_t, 2>(7, 3)});
    test_conversion_from_mdspan<Kokkos::LayoutLeft>(
        Kokkos::View<value_type **, Kokkos::LayoutLeft, ExecutionSpace>("ref",
                                                                        7, 3),
        typename layout_left_padded<7>::template mapping<
            Kokkos::extents<std::size_t, 7, 3>>{
            Kokkos::extents<std::size_t, 7, 3>()});

    // LayoutRight
    test_conversion_from_mdspan<Kokkos::LayoutRight>(
        Kokkos::View<value_type *, Kokkos::LayoutRight, ExecutionSpace>("ref",
                                                                        7),
        typename layout_right_padded<7>::template mapping<
            Kokkos::dextents<std::size_t, 1>>{
            Kokkos::dextents<std::size_t, 1>(7)});
    test_conversion_from_mdspan<Kokkos::LayoutRight>(
        Kokkos::View<value_type[7], Kokkos::LayoutRight, ExecutionSpace>("ref"),
        typename layout_right_padded<7>::template mapping<
            Kokkos::extents<std::size_t, 7>>{
            Kokkos::extents<std::size_t, 7>()});
    test_conversion_from_mdspan<Kokkos::LayoutRight>(
        Kokkos::View<value_type[7], Kokkos::LayoutRight, ExecutionSpace>("ref"),
        typename layout_right_padded<7>::template mapping<
            Kokkos::dextents<std::size_t, 1>>{
            Kokkos::dextents<std::size_t, 1>(7)});
    test_conversion_from_mdspan<Kokkos::LayoutRight>(
        Kokkos::View<value_type *, Kokkos::LayoutRight, ExecutionSpace>("ref",
                                                                        7),
        typename layout_right_padded<7>::template mapping<
            Kokkos::extents<std::size_t, 7>>{
            Kokkos::extents<std::size_t, 7>()});

    test_conversion_from_mdspan<Kokkos::LayoutRight>(
        Kokkos::View<value_type **, Kokkos::LayoutRight, ExecutionSpace>("ref",
                                                                         3, 7),
        typename layout_right_padded<7>::template mapping<
            Kokkos::dextents<std::size_t, 2>>{
            Kokkos::dextents<std::size_t, 2>(3, 7)});
    test_conversion_from_mdspan<Kokkos::LayoutRight>(
        Kokkos::View<value_type[3][7], Kokkos::LayoutRight, ExecutionSpace>(
            "ref"),
        typename layout_right_padded<7>::template mapping<
            Kokkos::extents<std::size_t, 3, 7>>{
            Kokkos::extents<std::size_t, 3, 7>()});
    test_conversion_from_mdspan<Kokkos::LayoutRight>(
        Kokkos::View<value_type[3][7], Kokkos::LayoutRight, ExecutionSpace>(
            "ref"),
        typename layout_right_padded<7>::template mapping<
            Kokkos::dextents<std::size_t, 2>>{
            Kokkos::dextents<std::size_t, 2>(3, 7)});
    test_conversion_from_mdspan<Kokkos::LayoutRight>(
        Kokkos::View<value_type **, Kokkos::LayoutRight, ExecutionSpace>("ref",
                                                                         3, 7),
        typename layout_right_padded<7>::template mapping<
            Kokkos::extents<std::size_t, 3, 7>>{
            Kokkos::extents<std::size_t, 3, 7>()});

    // LayoutStride
    {
      const size_t strides[] = {2};
      test_conversion_from_mdspan<Kokkos::LayoutStride>(
          Kokkos::View<value_type *, Kokkos::LayoutStride, ExecutionSpace>(
              "ref", Kokkos::LayoutStride{7, 2}),
          Kokkos::layout_stride::mapping<Kokkos::dextents<std::size_t, 1>>{
              Kokkos::mdspan_non_standard, Kokkos::dextents<std::size_t, 1>{7},
              strides});
    }
    {
      const size_t strides[] = {2};
      test_conversion_from_mdspan<Kokkos::LayoutStride>(
          Kokkos::View<value_type[7], Kokkos::LayoutStride, ExecutionSpace>(
              "ref", Kokkos::LayoutStride{7, 2}),
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 7>>{
              Kokkos::mdspan_non_standard, {}, strides});
    }
    {
      const size_t strides[] = {2};
      test_conversion_from_mdspan<Kokkos::LayoutStride>(
          Kokkos::View<value_type[7], Kokkos::LayoutStride, ExecutionSpace>(
              "ref", Kokkos::LayoutStride{7, 2}),
          Kokkos::layout_stride::mapping<Kokkos::dextents<std::size_t, 1>>{
              Kokkos::mdspan_non_standard, Kokkos::dextents<std::size_t, 1>{7},
              strides});
    }
    {
      const size_t strides[] = {2};
      test_conversion_from_mdspan<Kokkos::LayoutStride>(
          Kokkos::View<value_type *, Kokkos::LayoutStride, ExecutionSpace>(
              "ref", Kokkos::LayoutStride{7, 2}),
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 7>>{
              Kokkos::mdspan_non_standard, Kokkos::extents<std::size_t, 7>(),
              strides});
    }

    {
      const size_t strides[] = {2, 4};
      test_conversion_from_mdspan<Kokkos::LayoutStride>(
          Kokkos::View<value_type **, Kokkos::LayoutStride, ExecutionSpace>(
              "ref", Kokkos::LayoutStride{7, 2, 3, 4}),
          Kokkos::layout_stride::mapping<Kokkos::dextents<std::size_t, 2>>{
              Kokkos::mdspan_non_standard,
              Kokkos::dextents<std::size_t, 2>(7, 3), strides});
    }
    {
      const size_t strides[] = {2, 4};
      test_conversion_from_mdspan<Kokkos::LayoutStride>(
          Kokkos::View<value_type[7][3], Kokkos::LayoutStride, ExecutionSpace>(
              "ref", Kokkos::LayoutStride{7, 2, 3, 4}),
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 7, 3>>{
              Kokkos::mdspan_non_standard, Kokkos::extents<std::size_t, 7, 3>(),
              strides});
    }
    {
      const size_t strides[] = {2, 4};
      test_conversion_from_mdspan<Kokkos::LayoutStride>(
          Kokkos::View<value_type[7][3], Kokkos::LayoutStride, ExecutionSpace>(
              "ref", Kokkos::LayoutStride{7, 2, 3, 4}),
          Kokkos::layout_stride::mapping<Kokkos::dextents<std::size_t, 2>>{
              Kokkos::mdspan_non_standard,
              Kokkos::dextents<std::size_t, 2>(7, 3), strides});
    }
    {
      const size_t strides[] = {2, 4};
      test_conversion_from_mdspan<Kokkos::LayoutStride>(
          Kokkos::View<value_type **, Kokkos::LayoutStride, ExecutionSpace>(
              "ref", Kokkos::LayoutStride{7, 2, 3, 4}),
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 7, 3>>{
              Kokkos::mdspan_non_standard, Kokkos::extents<std::size_t, 7, 3>(),
              strides});
    }

    // Conversion to mdspan
    test_conversion_to_mdspan(
        layout_left_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4>>({}, 4),
        Kokkos::View<value_type *, Kokkos::LayoutLeft, ExecutionSpace>("v", 4));
    test_conversion_to_mdspan(
        layout_left_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4, 7>>({}, 4),
        Kokkos::View<value_type **, Kokkos::LayoutLeft, ExecutionSpace>("v", 4,
                                                                        7));

    test_conversion_to_mdspan(
        layout_right_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4>>({}, 4),
        Kokkos::View<value_type *, Kokkos::LayoutRight, ExecutionSpace>("v",
                                                                        4));
    test_conversion_to_mdspan(
        layout_right_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4, 7>>({}, 7),
        Kokkos::View<value_type **, Kokkos::LayoutRight, ExecutionSpace>("v", 4,
                                                                         7));

    {
      const size_t strides[] = {5};
      test_conversion_to_mdspan(
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 4>>(
              Kokkos::mdspan_non_standard, {}, strides),
          Kokkos::View<value_type *, Kokkos::LayoutStride, ExecutionSpace>(
              "v", Kokkos::LayoutStride{4, 5}));
    }
    {
      const size_t strides[] = {5, 9};
      test_conversion_to_mdspan(
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 4, 7>>(
              Kokkos::mdspan_non_standard, {}, strides),
          Kokkos::View<value_type **, Kokkos::LayoutStride, ExecutionSpace>(
              "v", Kokkos::LayoutStride{4, 5, 7, 9}));
    }

    // Aligned types (for padded layouts)
    test_conversion_to_mdspan(
        layout_left_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 127, 7>>({}, 128),
        Kokkos::View<value_type **, Kokkos::LayoutLeft, ExecutionSpace>(
            Kokkos::view_alloc("v", Kokkos::AllowPadding), 127, 7));

    test_conversion_to_mdspan(
        layout_right_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 7, 127>>({}, 128),
        Kokkos::View<value_type **, Kokkos::LayoutRight, ExecutionSpace>(
            Kokkos::view_alloc("v", Kokkos::AllowPadding), 7, 127));

    // Conversion with standard default_accessor

    test_conversion_to_mdspan_with_accessor(
        layout_left_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4>>({}, 4),
        Kokkos::View<value_type *, Kokkos::LayoutLeft, ExecutionSpace>("v", 4),
        Kokkos::default_accessor<value_type>{});
    test_conversion_to_mdspan_with_accessor(
        layout_left_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4, 7>>({}, 4),
        Kokkos::View<value_type **, Kokkos::LayoutLeft, ExecutionSpace>("v", 4,
                                                                        7),
        Kokkos::default_accessor<value_type>{});

    test_conversion_to_mdspan_with_accessor(
        layout_right_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4>>({}, 4),
        Kokkos::View<value_type *, Kokkos::LayoutRight, ExecutionSpace>("v", 4),
        Kokkos::default_accessor<value_type>{});
    test_conversion_to_mdspan_with_accessor(
        layout_right_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4, 7>>({}, 7),
        Kokkos::View<value_type **, Kokkos::LayoutRight, ExecutionSpace>("v", 4,
                                                                         7),
        Kokkos::default_accessor<value_type>{});

    {
      const size_t strides[] = {5};
      test_conversion_to_mdspan_with_accessor(
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 4>>(
              Kokkos::mdspan_non_standard, {}, strides),
          Kokkos::View<value_type *, Kokkos::LayoutStride, ExecutionSpace>(
              "v", Kokkos::LayoutStride{4, 5}),
          Kokkos::default_accessor<value_type>{});
    }
    {
      const size_t strides[] = {5, 9};
      test_conversion_to_mdspan_with_accessor(
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 4, 7>>(
              Kokkos::mdspan_non_standard, {}, strides),
          Kokkos::View<value_type **, Kokkos::LayoutStride, ExecutionSpace>(
              "v", Kokkos::LayoutStride{4, 5, 7, 9}),
          Kokkos::default_accessor<value_type>{});
    }

    // Conversion with a test accessor

    test_conversion_to_mdspan_with_accessor(
        layout_left_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4>>({}, 4),
        Kokkos::View<value_type *, Kokkos::LayoutLeft, ExecutionSpace>("v", 4),
        TestAccessor{});
    test_conversion_to_mdspan_with_accessor(
        layout_left_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4, 7>>({}, 4),
        Kokkos::View<value_type **, Kokkos::LayoutLeft, ExecutionSpace>("v", 4,
                                                                        7),
        TestAccessor{});

    test_conversion_to_mdspan_with_accessor(
        layout_right_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4>>({}, 4),
        Kokkos::View<value_type *, Kokkos::LayoutRight, ExecutionSpace>("v", 4),
        TestAccessor{});
    test_conversion_to_mdspan_with_accessor(
        layout_right_padded<Kokkos::dynamic_extent>::mapping<
            Kokkos::extents<std::size_t, 4, 7>>({}, 7),
        Kokkos::View<value_type **, Kokkos::LayoutRight, ExecutionSpace>("v", 4,
                                                                         7),
        TestAccessor{});

    {
      const size_t strides[] = {5};
      test_conversion_to_mdspan_with_accessor(
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 4>>(
              Kokkos::mdspan_non_standard, {}, strides),
          Kokkos::View<value_type *, Kokkos::LayoutStride, ExecutionSpace>(
              "v", Kokkos::LayoutStride{4, 5}),
          TestAccessor{});
    }
    {
      const size_t strides[] = {5, 9};
      test_conversion_to_mdspan_with_accessor(
          Kokkos::layout_stride::mapping<Kokkos::extents<std::size_t, 4, 7>>(
              Kokkos::mdspan_non_standard, {}, strides),
          Kokkos::View<value_type **, Kokkos::LayoutStride, ExecutionSpace>(
              "v", Kokkos::LayoutStride{4, 5, 7, 9}),
          TestAccessor{});
    }
  }
};

TEST(TEST_CATEGORY, view_mdspan_conversion) {
  TestViewMDSpanConversion<double, TEST_EXECSPACE>::run_test();
  TestViewMDSpanConversion<float, TEST_EXECSPACE>::run_test();
  TestViewMDSpanConversion<int, TEST_EXECSPACE>::run_test();
}

}  // namespace
