// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif
#include <Kokkos_TypeInfo.hpp>

#include <gtest/gtest.h>

namespace Test {
namespace Impl {
template <class ViewTypeDst, class ViewTypeSrc>
struct TestAssignability {
  using mapping_type =
      Kokkos::Impl::ViewMapping<typename ViewTypeDst::traits,
                                typename ViewTypeSrc::traits
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
                                ,
                                typename ViewTypeDst::specialize
#else
                                ,
                                void
#endif
                                >;

  template <class MappingType>
  static void try_assign(
      ViewTypeDst& dst, ViewTypeSrc& src,
      std::enable_if_t<MappingType::is_assignable>* = nullptr) {
    dst = src;
  }

  template <class MappingType>
  static void try_assign(
      ViewTypeDst&, ViewTypeSrc&,
      std::enable_if_t<!MappingType::is_assignable>* = nullptr) {
    FAIL() << "TestAssignability::try_assign: Unexpected call path";
  }

  template <class... Dimensions>
  static void test(bool always, bool sometimes, Dimensions... dims) {
    ViewTypeDst dst;
    ViewTypeSrc src("SRC", dims...);

    bool is_always_assignable =
        Kokkos::is_always_assignable<ViewTypeDst, ViewTypeSrc>::value;
    bool is_assignable = Kokkos::is_assignable(dst, src);

    if (sometimes) {
      try_assign<mapping_type>(dst, src);
    }
    ASSERT_EQ(always, is_always_assignable)
        << Kokkos::Impl::TypeInfo<ViewTypeSrc>::name() << " to "
        << Kokkos::Impl::TypeInfo<ViewTypeDst>::name();
    ASSERT_EQ(sometimes, is_assignable)
        << Kokkos::Impl::TypeInfo<ViewTypeSrc>::name() << " to "
        << Kokkos::Impl::TypeInfo<ViewTypeDst>::name();
  }
};

}  // namespace Impl

TEST(TEST_CATEGORY, view_is_assignable) {
  using namespace Kokkos;
  using h_exec = typename DefaultHostExecutionSpace::memory_space;
  using d_exec = typename TEST_EXECSPACE::memory_space;
  using left   = LayoutLeft;
  using right  = LayoutRight;
  using stride = LayoutStride;
  // Static/Dynamic Extents
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<int*, left, d_exec>>::test(true, true, 10);
  Impl::TestAssignability<View<int[10], left, d_exec>,
                          View<int*, left, d_exec>>::test(false, true, 10);
  Impl::TestAssignability<View<int[5], left, d_exec>,
                          View<int*, left, d_exec>>::test(false, false, 10);
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<int[10], left, d_exec>>::test(true, true);
  Impl::TestAssignability<View<int[10], left, d_exec>,
                          View<int[10], left, d_exec>>::test(true, true);
  Impl::TestAssignability<View<int[5], left, d_exec>,
                          View<int[10], left, d_exec>>::test(false, false);
  Impl::TestAssignability<View<int**, left, d_exec>,
                          View<int**, left, d_exec>>::test(true, true, 10, 10);
  Impl::TestAssignability<View<int* [10], left, d_exec>,
                          View<int**, left, d_exec>>::test(false, true, 10, 10);
  Impl::TestAssignability<View<int* [5], left, d_exec>,
                          View<int**, left, d_exec>>::test(false, false, 10,
                                                           10);
  Impl::TestAssignability<View<int**, left, d_exec>,
                          View<int* [10], left, d_exec>>::test(true, true, 10);
  Impl::TestAssignability<View<int* [10], left, d_exec>,
                          View<int* [10], left, d_exec>>::test(true, true, 10);
  Impl::TestAssignability<View<int* [5], left, d_exec>,
                          View<int* [10], left, d_exec>>::test(false, false,
                                                               10);

  // Mismatch value_type
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<double*, left, d_exec>>::test(false, false, 10);

  // Layout assignment
  Impl::TestAssignability<View<int, left, d_exec>,
                          View<int, right, d_exec>>::test(true, true);
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<int*, right, d_exec>>::test(true, true, 10);
  Impl::TestAssignability<View<int[5], left, d_exec>,
                          View<int*, right, d_exec>>::test(false, false, 10);
  Impl::TestAssignability<View<int[10], left, d_exec>,
                          View<int*, right, d_exec>>::test(false, true, 10);
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<int[5], right, d_exec>>::test(true, true);
  Impl::TestAssignability<View<int[5], left, d_exec>,
                          View<int[10], right, d_exec>>::test(false, false);

  // This could be made possible (due to the degenerate nature of the views) but
  // we do not allow this yet
  // TestAssignability<View<int**,left,d_exec>,View<int**,right,d_exec>>::test(false,true,10,1);
  Impl::TestAssignability<View<int**, left, d_exec>,
                          View<int**, right, d_exec>>::test(false, false, 10,
                                                            2);
  Impl::TestAssignability<View<int**, stride, d_exec>,
                          View<int**, right, d_exec>>::test(true, true, 10, 2);
  Impl::TestAssignability<View<int**, stride, d_exec>,
                          View<int**, left, d_exec>>::test(true, true, 10, 2);

  // Space Assignment
  bool expected = Kokkos::Impl::MemorySpaceAccess<d_exec, h_exec>::assignable;
  Impl::TestAssignability<View<int*, left, d_exec>,
                          View<int*, left, h_exec>>::test(expected, expected,
                                                          10);
  expected = Kokkos::Impl::MemorySpaceAccess<h_exec, d_exec>::assignable;
  Impl::TestAssignability<View<int*, left, h_exec>,
                          View<int*, left, d_exec>>::test(expected, expected,
                                                          10);

  // reference type and const-qualified types
  using SomeViewType = View<int*, left, d_exec>;
  static_assert(is_always_assignable_v<SomeViewType, SomeViewType>);
  static_assert(is_always_assignable_v<SomeViewType, SomeViewType&>);
  static_assert(is_always_assignable_v<SomeViewType, SomeViewType const>);
  static_assert(is_always_assignable_v<SomeViewType, SomeViewType const&>);
  static_assert(is_always_assignable_v<SomeViewType&, SomeViewType>);
  static_assert(is_always_assignable_v<SomeViewType&, SomeViewType&>);
  static_assert(is_always_assignable_v<SomeViewType&, SomeViewType const>);
  static_assert(is_always_assignable_v<SomeViewType&, SomeViewType const&>);

#ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
  // Check assignment to const qualified Views is false
  static_assert(!is_always_assignable_v<SomeViewType const, SomeViewType>);
  static_assert(!is_always_assignable_v<SomeViewType const&, SomeViewType>);
  {
    SomeViewType non_const_view("V", 10);
    const SomeViewType const_view;
    ASSERT_FALSE(Kokkos::is_assignable(const_view, non_const_view));
  }
#endif

  // Rank mismatch: is_always_assignable_v must be false (consteval false case)
  using RankMismatchView1D = View<int*, left, d_exec>;
  using RankMismatchView2D = View<int**, left, d_exec>;
  static_assert(
      !is_always_assignable_v<RankMismatchView1D, RankMismatchView2D>);
  static_assert(
      !is_always_assignable_v<RankMismatchView2D, RankMismatchView1D>);

  // Rank mismatch: is_assignable must compile and return false (not cause a
  // compile error due to a non-static impl_runtime_value in the primary
  // template of is_assignable_impl)
  {
    RankMismatchView1D v1d;
    RankMismatchView2D v2d("v2d", 10, 10);
    ASSERT_FALSE(Kokkos::is_assignable(v1d, v2d));
    ASSERT_FALSE(Kokkos::is_assignable(v2d, v1d));
  }

  // CV-qualified element types: adding const to the element type is allowed;
  // removing const is not
  using ConstElementView = View<const int*, left, d_exec>;
  static_assert(is_always_assignable_v<ConstElementView, SomeViewType>);
  static_assert(!is_always_assignable_v<SomeViewType, ConstElementView>);
  {
    ConstElementView cv;
    SomeViewType v("v", 10);
    ASSERT_TRUE(Kokkos::is_assignable(cv, v));
    ASSERT_FALSE(Kokkos::is_assignable(v, cv));
  }
}
}  // namespace Test
