// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif
#include <cstddef>

namespace {

template <class NewExts, class OldExts>
constexpr bool check_expected_static_extent_match() {
  static_assert(OldExts::rank() == NewExts::rank());

  // Only check the last extents - old style Views don't support arbitrary
  // mixing of static and dynamic extents; thus the new extents could contain
  // more static values, but never less.
  static_assert(OldExts::rank_dynamic() >= NewExts::rank_dynamic());

  // Silence warnings about pointless comparison
  if constexpr (OldExts::rank() > 0) {
    for (size_t r = OldExts::rank_dynamic(); r < OldExts::rank(); r++)
      if (OldExts::static_extent(r) != NewExts::static_extent(r)) return false;
  }
  return true;
}

template <class Layout, class... Slices>
constexpr bool new_subview_retains_padded_layout_old_uses_strided() {
  if constexpr (std::is_same_v<Layout,
                               Kokkos::Experimental::layout_right_padded<>>) {
    if constexpr (sizeof...(Slices) == 4 &&
                  std::is_same_v<
                      std::tuple<Slices...>,
                      std::tuple<Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::ALL_t,
                                 Kokkos::pair<int, int>>>)
      return true;
    if constexpr (sizeof...(Slices) == 8 &&
                  std::is_same_v<
                      std::tuple<Slices...>,
                      std::tuple<Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::ALL_t,
                                 Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::ALL_t,
                                 Kokkos::ALL_t, Kokkos::pair<int, int>>>)
      return true;
  }
  return false;
}

// Compare behavior of View with mdspan style args with View with classic args
template <class V, class... Slices>
void compare_subview(V v, Slices... slices) {
  using old_view_t =
      Kokkos::View<typename V::data_type, typename V::array_layout,
                   typename V::memory_space, typename V::memory_traits>;
  old_view_t old_v = v;

  auto old_sub = Kokkos::subview(old_v, slices...);
  auto new_sub = Kokkos::subview(v, slices...);

  using old_sub_t = decltype(old_sub);
  using new_sub_t = decltype(new_sub);

  using mdspan_sub_t = decltype(Kokkos::submdspan(
      std::declval<typename V::mdspan_type>(),
      Kokkos::Impl::transform_kokkos_slice_to_mdspan_slice(
          Kokkos::Impl::convert_to_kokkos_pair_if_std_pair(slices))...));
  static_assert(std::is_same_v<mdspan_sub_t, typename new_sub_t::mdspan_type>);

  static_assert(
      std::is_same_v<typename new_sub_t::index_type, typename V::index_type>);
  static_assert(
      check_expected_static_extent_match<typename new_sub_t::extents_type,
                                         typename old_sub_t::extents_type>());

  ASSERT_EQ(old_sub.extents(), new_sub.extents());
  ASSERT_EQ(old_sub.data(), new_sub.data());
  // Comparing the return type is a bit iffy for rank 1 (and technically 0)
  // mdspan drops down to plain layout_left/right instead of returning padded
  // layouts when the original View had padded layouts since for rank-0/1 there
  // can't be padding. Thus the new View will also drop to those guys, while the
  // old View maintains padded layouts.
  if constexpr (new_sub_t::rank() < 2) {
    if constexpr (std::is_same_v<typename old_sub_t::layout_type,
                                 Kokkos::Experimental::layout_right_padded<>>) {
      static_assert(std::is_same_v<Kokkos::layout_right,
                                   typename new_sub_t::layout_type>);
    } else if constexpr (std::is_same_v<
                             typename old_sub_t::layout_type,
                             Kokkos::Experimental::layout_left_padded<>>) {
      static_assert(
          std::is_same_v<Kokkos::layout_left, typename new_sub_t::layout_type>);
    } else {
      static_assert(std::is_same_v<typename old_sub_t::layout_type,
                                   typename new_sub_t::layout_type>);
    }
  } else {  // new_sub_t::rank() >= 2
    // Deal with the case where submdspan and thus subview with new style args
    // optimizes return type better and does not fall back to layout_stride.
    if constexpr (new_subview_retains_padded_layout_old_uses_strided<
                      typename V::layout_type, Slices...>()) {
      static_assert(!std::is_same_v<typename old_sub_t::layout_type,
                                    typename new_sub_t::layout_type>);
      static_assert(std::is_same_v<typename old_sub_t::layout_type,
                                   Kokkos::layout_stride>);
    } else {
      static_assert(std::is_same_v<typename old_sub_t::layout_type,
                                   typename new_sub_t::layout_type>);
    }
  }
}

// TODO: full_extent is not yet supported as an argument for subview with the
// old View style

constexpr std::pair<int, int> pair00    = {0, 0};
constexpr Kokkos::pair<int, int> pair37 = {3, 7};
constexpr Kokkos::ALL_t all{};

template <class V>
  requires(V::rank() == 1)
void test_subview_args(V v) {
  compare_subview(v, all);
  compare_subview(v, pair00);
  if (v.extent(0) != 0) {
    compare_subview(v, 0);
    compare_subview(v, pair37);
  }
}

template <class V>
  requires(V::rank() == 2)
void test_subview_args(V v) {
  compare_subview(v, pair00, 1);
  compare_subview(v, pair00, pair00);
  compare_subview(v, pair00, all);
  compare_subview(v, all, 1);
  compare_subview(v, all, pair00);
  compare_subview(v, all, all);
  if (v.extent(0) != 0) {
    compare_subview(v, 0, 1);
    compare_subview(v, 0, pair37);
    compare_subview(v, 0, all);
    compare_subview(v, pair37, 1);
    compare_subview(v, pair37, pair37);
    compare_subview(v, pair37, all);
    compare_subview(v, all, 1);
    compare_subview(v, all, pair37);
    compare_subview(v, all, all);
  }
}

template <class V>
  requires(V::rank() == 4)
void test_subview_args(V v) {
  compare_subview(v, pair00, 3, all, 1);
  compare_subview(v, pair00, pair37, 4, pair00);
  compare_subview(v, pair00, all, all, all);
  compare_subview(v, all, 1, 3, 4);
  compare_subview(v, all, 2, pair00, 1);
  compare_subview(v, all, pair37, all, all);
  if (v.extent(0) != 0) {
    compare_subview(v, 0, all, pair37, 1);
    compare_subview(v, 0, 3, 3, pair37);
    compare_subview(v, 0, pair37, pair37, all);
    compare_subview(v, pair37, all, all, 1);
    compare_subview(v, pair37, 3, 3, pair37);
    compare_subview(v, pair37, all, all, all);
    compare_subview(v, all, pair37, all, 1);
    compare_subview(v, all, all, all, pair37);
    compare_subview(v, all, 1, 3, all);
  }
}

template <class V>
  requires(V::rank() == 8)
void test_subview_args(V v) {
  compare_subview(v, pair00, 3, 3, all, 1, 1, 1, 1);
  compare_subview(v, pair00, 3, all, pair37, all, pair37, 4, pair00);
  compare_subview(v, pair00, all, all, all, all, all, all, all);
  compare_subview(v, all, 1, 2, 3, 4, 5, 6, 6);
  compare_subview(v, all, 2, pair00, all, 3, 4, pair37, 1);
  compare_subview(v, all, pair37, all, all, all, all, all, all);
  if (v.extent(0) != 0) {
    compare_subview(v, 0, all, pair37, 3, 5, all, pair37, 1);
    compare_subview(v, 0, 3, 3, 4, 5, 1, 3, pair37);
    compare_subview(v, 0, all, all, 3, pair37, pair37, pair37, all);
    compare_subview(v, pair37, all, pair37, all, 3, pair37, all, 1);
    compare_subview(v, pair37, 3, all, 3, all, pair37, 3, pair37);
    compare_subview(v, pair37, all, all, all, all, all, all, all);
    compare_subview(v, all, pair37, pair37, all, 3, 5, all, 1);
// Newer GCC produces warnings deep inside mdspan for these cases
// The warnings are a false positive where access to a staticly
// sized array with an r>=size is only protected by a runtime
// condition.
#if !(defined(KOKKOS_COMPILER_GNU) && __GNUC__ > 11)
    compare_subview(v, all, all, all, all, all, all, all, pair37);
    compare_subview(v, all, all, 1, all, all, pair37, 3, all);
#endif
  }
}

constexpr size_t dyn = Kokkos::dynamic_extent;

template <class E, class A>
void test_subview_expand(const E& exts, const A& acc) {
  using T    = typename A::element_type;
  auto alloc = Kokkos::view_alloc("TestView");
  {
    using layout_t = Kokkos::Experimental::layout_left_padded<dyn>;
    using map_t    = typename layout_t::mapping<E>;
    test_subview_args(Kokkos::View<T, E, layout_t, A>(alloc, map_t(exts), acc));
  }
  {
    using layout_t = Kokkos::Experimental::layout_right_padded<dyn>;
    using map_t    = typename layout_t::mapping<E>;
    test_subview_args(Kokkos::View<T, E, layout_t, A>(alloc, map_t(exts), acc));
  }
}
}  // namespace

TEST(TEST_CATEGORY, view_mdspan_args_subview) {
  using mem_t = typename TEST_EXECSPACE::memory_space;
  // rank 1
  test_subview_expand(
      Kokkos::dextents<int, 1>(10),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::dextents<int64_t, 1>(0),
      Kokkos::Experimental::Accessor<float, mem_t,
                                     Kokkos::MemoryTraits<Kokkos::Atomic>>());
  test_subview_expand(
      Kokkos::extents<unsigned, 10>(),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::extents<unsigned, 0>(),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());

  // rank 2
  test_subview_expand(
      Kokkos::dextents<int, 2>(10, 20),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::dextents<int64_t, 2>(0, 20),
      Kokkos::Experimental::Accessor<float, mem_t,
                                     Kokkos::MemoryTraits<Kokkos::Atomic>>());
  test_subview_expand(
      Kokkos::extents<unsigned, dyn, 10>(10),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::extents<unsigned, 0, 20>(),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::extents<uint64_t, dyn, 20>(10),
      Kokkos::Experimental::Accessor<float, mem_t,
                                     Kokkos::MemoryTraits<Kokkos::Atomic>>());

  // rank 4
  test_subview_expand(
      Kokkos::dextents<int, 4>(10, 20, 30, 40),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::dextents<int64_t, 4>(0, 20, 30, 40),
      Kokkos::Experimental::Accessor<float, mem_t,
                                     Kokkos::MemoryTraits<Kokkos::Atomic>>());
  test_subview_expand(
      Kokkos::extents<unsigned, dyn, dyn, dyn, 40>(10, 20, 30),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::extents<unsigned, 0, 20, 30, 40>(),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::extents<uint64_t, dyn, 20, 30, 40>(10),
      Kokkos::Experimental::Accessor<float, mem_t,
                                     Kokkos::MemoryTraits<Kokkos::Atomic>>());

  // rank 8
  // using 7 for every extent so we don't use too much memory and its the
  // minimum size to work with pair37
  test_subview_expand(
      Kokkos::dextents<int, 8>(7, 7, 7, 7, 7, 7, 7, 7),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::dextents<int64_t, 8>(0, 7, 7, 7, 7, 7, 7, 7),
      Kokkos::Experimental::Accessor<float, mem_t,
                                     Kokkos::MemoryTraits<Kokkos::Atomic>>());
  test_subview_expand(
      Kokkos::extents<unsigned, dyn, dyn, dyn, dyn, dyn, dyn, dyn, 7>(
          7, 7, 7, 7, 7, 7, 7),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::extents<unsigned, 0, 7, 7, 7, 7, 7, 7, 7>(),
      Kokkos::Experimental::Accessor<float, mem_t, Kokkos::MemoryTraits<>>());
  test_subview_expand(
      Kokkos::extents<uint64_t, dyn, dyn, dyn, 7, 7, 7, 7, 7>(7, 7, 7),
      Kokkos::Experimental::Accessor<float, mem_t,
                                     Kokkos::MemoryTraits<Kokkos::Atomic>>());
}
