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
#include <type_traits>

using ExecutionSpace = TEST_EXECSPACE;

namespace {
template <class ExecutionSpace, class Extents>
auto make_spanning_mdrange_policy_from_extents_impl(const Extents &extents,
                                                    std::index_sequence<0>) {
  return Kokkos::RangePolicy<ExecutionSpace>{0, extents.extent(0)};
}

template <class ExecutionSpace, class Extents, std::size_t... Indices>
auto make_spanning_mdrange_policy_from_extents_impl(
    const Extents &extents, std::index_sequence<Indices...>) {
  using index_type    = typename Extents::index_type;
  constexpr auto rank = Extents::rank();
  return Kokkos::MDRangePolicy<ExecutionSpace, Kokkos::Rank<rank>>{
      {(static_cast<index_type>(Indices * 0))...},
      {extents.extent(Indices)...}};
}

template <class ExecutionSpace, class Extents>
auto make_spanning_mdrange_policy_from_extents(const Extents &extents) {
  return make_spanning_mdrange_policy_from_extents_impl<ExecutionSpace>(
      extents, std::make_index_sequence<Extents::rank()>{});
}

template <class T, class ExtentsType>
void test_default_constructor() {
  using extents_type  = ExtentsType;
  using layout_type   = Kokkos::Experimental::layout_right_padded<>;
  using accessor_type = Kokkos::Impl::CheckedReferenceCountedAccessor<
      T, typename ExecutionSpace::memory_space>;
  using view_type =
      Kokkos::Impl::BV::BasicView<T, extents_type, layout_type, accessor_type>;
  view_type view;

  EXPECT_FALSE(view.data_handle().has_record());
  EXPECT_EQ(view.data_handle().get(), nullptr);
  EXPECT_EQ(view.extents(), extents_type{});
  EXPECT_EQ(view.data_handle().use_count(), 0);
  EXPECT_TRUE(view.is_exhaustive());
  EXPECT_EQ(view.data_handle().get_label(), "");
  EXPECT_TRUE(view.empty());
  EXPECT_EQ(view.size(), 0u);
}

TEST(TEST_CATEGORY, basic_view_default_ctor) {
  test_default_constructor<double, Kokkos::extents<std::size_t, 1>>();
}

template <class T, class ExtentsType>
void test_extents_constructor(const ExtentsType &extents) {
  using extents_type  = ExtentsType;
  using layout_type   = Kokkos::Experimental::layout_right_padded<>;
  using accessor_type = Kokkos::Impl::CheckedReferenceCountedAccessor<
      T, typename ExecutionSpace::memory_space>;
  using view_type =
      Kokkos::Impl::BV::BasicView<T, extents_type, layout_type, accessor_type>;

  view_type view("test_view", extents);

  EXPECT_TRUE(view.data_handle().has_record());
  EXPECT_NE(view.data_handle().get(), nullptr);
  EXPECT_EQ(view.extents(), extents);
  EXPECT_EQ(view.data_handle().use_count(), 1);
  EXPECT_TRUE(view.is_exhaustive());
  EXPECT_EQ(view.data_handle().get_label(), "test_view");
  size_t expected_size = 1;
  // Avoid pointless comparison of unsigned warning for rank==0
  for (int r = 0; r < static_cast<int>(view_type::rank()); r++)
    expected_size *= extents.extent(r);
  EXPECT_EQ(view.size(), expected_size);
  EXPECT_EQ(view.empty(), expected_size == 0u);
}

TEST(TEST_CATEGORY, basic_view_extents_ctor) {
  test_extents_constructor<double>(
      Kokkos::extents<std::size_t, 2, Kokkos::dynamic_extent, 4>(8));
  test_extents_constructor<double>(
      Kokkos::extents<std::size_t, 2, Kokkos::dynamic_extent, 4>(0));
  test_extents_constructor<std::size_t>(Kokkos::extents<std::size_t, 2, 4>());
  test_extents_constructor<int>(Kokkos::extents<std::size_t>());
}

template <class T, template <std::size_t> class LayoutType, class ExtentsType>
void test_mapping_constructor(const ExtentsType &extents, std::size_t padding) {
  using extents_type  = ExtentsType;
  using layout_type   = LayoutType<Kokkos::dynamic_extent>;
  using mapping_type  = typename layout_type::template mapping<ExtentsType>;
  using accessor_type = Kokkos::Impl::CheckedReferenceCountedAccessor<
      T, typename ExecutionSpace::memory_space>;
  using view_type =
      Kokkos::Impl::BV::BasicView<T, extents_type, layout_type, accessor_type>;
  static_assert(std::is_same_v<typename view_type::mapping_type, mapping_type>);

  auto mapping = mapping_type(extents, padding);

  view_type view("test_view", mapping);

  EXPECT_TRUE(view.data_handle().has_record());
  EXPECT_NE(view.data_handle().get(), nullptr);
  EXPECT_EQ(view.data_handle().use_count(), 1);
  EXPECT_EQ(view.data_handle().get_label(), "test_view");
  EXPECT_EQ(view.extents(), mapping.extents());
  EXPECT_EQ(view.is_exhaustive(), mapping.is_exhaustive());
  size_t expected_size = 1;
  // Avoid pointless comparison of unsigned warning for rank==0
  for (int r = 0; r < static_cast<int>(view_type::rank()); r++)
    expected_size *= view.extent(r);
  EXPECT_EQ(view.size(), expected_size);
  EXPECT_EQ(view.empty(), expected_size == 0u);
}

TEST(TEST_CATEGORY, basic_view_mapping_ctor_right) {
  test_mapping_constructor<double, Kokkos::Experimental::layout_left_padded>(
      Kokkos::extents<std::size_t, 2, Kokkos::dynamic_extent>(2, 5), 8);
  test_mapping_constructor<std::size_t,
                           Kokkos::Experimental::layout_left_padded>(
      Kokkos::extents<std::size_t>(), 4);
  test_mapping_constructor<double, Kokkos::Experimental::layout_left_padded>(
      Kokkos::extents<std::size_t, 2, 3>(), 9);
  test_mapping_constructor<int, Kokkos::Experimental::layout_right_padded>(
      Kokkos::extents<std::size_t, 2, Kokkos::dynamic_extent>(2, 5), 8);
  test_mapping_constructor<double, Kokkos::Experimental::layout_right_padded>(
      Kokkos::extents<std::size_t>(), 4);
  test_mapping_constructor<unsigned, Kokkos::Experimental::layout_right_padded>(
      Kokkos::extents<std::size_t, 2, 3>(), 9);
}

template <class ViewType>
struct MDRangeTestFunctor {
  ViewType view;
  template <class... Idxs>
  KOKKOS_FUNCTION void operator()(Idxs... idxs) const {
    view(idxs...) = (idxs + ...);
  }
};

template <class T, class LayoutType, class ExtentsType>
void test_access_with_extents(const ExtentsType &extents) {
  using extents_type  = ExtentsType;
  using layout_type   = Kokkos::Experimental::layout_right_padded<>;
  using accessor_type = Kokkos::Impl::CheckedReferenceCountedAccessor<
      T, typename ExecutionSpace::memory_space>;
  using view_type =
      Kokkos::Impl::BV::BasicView<T, extents_type, layout_type, accessor_type>;

  auto view = view_type("test_view", extents);

  EXPECT_TRUE(view.data_handle().has_record());
  EXPECT_NE(view.data_handle().get(), nullptr);

  auto mdrange_policy =
      make_spanning_mdrange_policy_from_extents<ExecutionSpace>(extents);

  Kokkos::parallel_for(mdrange_policy, MDRangeTestFunctor<view_type>{view});
}

template <class T, class LayoutType>
void test_access() {
  test_access_with_extents<T, LayoutType>(Kokkos::extents<std::size_t, 5>());
  test_access_with_extents<T, LayoutType>(
      Kokkos::extents<std::size_t, 5, 10>());
  test_access_with_extents<T, LayoutType>(
      Kokkos::extents<std::size_t, 5, 2, 2, 2, 2, 2>());
}

TEST(TEST_CATEGORY, basic_view_access) {
  test_access<double, Kokkos::Experimental::layout_left_padded<
                          Kokkos::dynamic_extent>>();
  test_access<std::size_t, Kokkos::Experimental::layout_right_padded<
                               Kokkos::dynamic_extent>>();
}

#if 0  // TODO: this test should be active after View is put on top of BasicView
  template <class T, template <std::size_t> class LayoutType, class SrcViewType,
            class ExtentsType>
  void test_construct_from_view(const ExtentsType &extents,
                                       std::size_t padding) {
    using extents_type  = ExtentsType;
    using layout_type   = LayoutType<Kokkos::dynamic_extent>;
    using mapping_type  = typename layout_type::template mapping<ExtentsType>;
    using accessor_type = Kokkos::Impl::CheckedReferenceCountedAccessor<
        T, typename ExecutionSpace::memory_space>;
    using basic_view_type =
        Kokkos::Impl::BV::BasicView<T, extents_type, layout_type, accessor_type>;
    using view_type = SrcViewType;
    static_assert(std::is_constructible_v<basic_view_type, SrcViewType>);
  }
#endif

#if 0  // TODO: this test should be active after View is put on top of BasicView
TEST(TEST_CATEGORY, basic_view_view_ctor) {
    test_construct_from_view<double,
        Kokkos::Experimental::layout_left_padded,
        Kokkos::View<double[3], Kokkos::LayoutLeft, ExecutionSpace>>(
        Kokkos::extents<std::size_t, 3>(), 0);

    test_construct_from_view<size_t,
        Kokkos::Experimental::layout_left_padded,
        Kokkos::View<double[3], Kokkos::LayoutLeft, ExecutionSpace>>(
        Kokkos::extents<std::size_t, Kokkos::dynamic_extent>(3), 0);
}
#endif

template <class T>
void test_atomic_accessor() {
  using extents_type = Kokkos::extents<int, 10, 12, 30>;
  using layout_type  = Kokkos::Experimental::layout_right_padded<>;
  using accessor_type =
      Kokkos::Impl::CheckedReferenceCountedRelaxedAtomicAccessor<
          T, typename ExecutionSpace::memory_space>;
  using view_type =
      Kokkos::Impl::BV::BasicView<T, extents_type, layout_type, accessor_type>;
  using um_accessor_type = Kokkos::Impl::CheckedRelaxedAtomicAccessor<
      T, typename ExecutionSpace::memory_space>;
  using um_view_type = Kokkos::Impl::BV::BasicView<T, extents_type, layout_type,
                                                   um_accessor_type>;

  extents_type extents{};
  auto view = view_type("test_view", extents);
  um_view_type um_view(view);

  EXPECT_TRUE(view.data_handle().has_record());
  EXPECT_NE(view.data_handle().get(), nullptr);

  auto mdrange_policy =
      make_spanning_mdrange_policy_from_extents<ExecutionSpace>(extents);

  Kokkos::parallel_for(mdrange_policy, MDRangeTestFunctor<view_type>{view});
  Kokkos::parallel_for(mdrange_policy,
                       MDRangeTestFunctor<um_view_type>{um_view});
}

TEST(TEST_CATEGORY, basic_view_atomic_accessor) {
  test_atomic_accessor<int>();
  test_atomic_accessor<double>();
// FIXME OPENACC atomics
#ifndef KOKKOS_ENABLE_OPENACC
  test_atomic_accessor<Kokkos::complex<double>>();
#endif
}

}  // namespace
