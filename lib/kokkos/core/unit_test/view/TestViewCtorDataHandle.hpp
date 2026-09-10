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

namespace {

template <class DeviceType>
struct TestViewCtorDataHandle {
  static void test_view_data_handle_ctor_ref_counts() {
#ifndef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
    Kokkos::View<int *> a("A", 5);
    ASSERT_EQ(a.use_count(), 1);
    {
      Kokkos::View<int *> b(a.data_handle(), a.extents());
      ASSERT_EQ(a.use_count(), 2);
      ASSERT_EQ(b.use_count(), 2);
    }
    ASSERT_EQ(a.use_count(), 1);
    {
      Kokkos::View<int *> b(a.data_handle(), a.mapping());
      ASSERT_EQ(a.use_count(), 2);
      ASSERT_EQ(b.use_count(), 2);
    }
    ASSERT_EQ(a.use_count(), 1);
    {
      Kokkos::View<int *> b(a.data_handle(), a.mapping(), a.accessor());
      ASSERT_EQ(a.use_count(), 2);
      ASSERT_EQ(b.use_count(), 2);
    }
    ASSERT_EQ(a.use_count(), 1);

    // This is not currently supported but may be added in the future
    //{
    // Kokkos::View<int *> b(a.data_handle(), 5);
    // ASSERT_EQ(a.use_count(), 2);
    // ASSERT_EQ(b.use_count(), 2);
    //}
    static_assert(!std::constructible_from<Kokkos::View<int *>,
                                           decltype(a.data_handle()), int>);
#endif
  }
};

TEST(TEST_CATEGORY, view_ctor_data_handle) {
  TestViewCtorDataHandle<
      TEST_EXECSPACE>::test_view_data_handle_ctor_ref_counts();
}

}  // namespace
