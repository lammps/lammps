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
#ifndef TESTREALLOC_HPP_
#define TESTREALLOC_HPP_

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace TestViewRealloc {

struct Default {};
struct WithoutInitializing {};

template <typename View, typename... Args>
inline void realloc_dispatch(Default, View& v, Args&&... args) {
  Kokkos::realloc(v, std::forward<Args>(args)...);
}

template <typename View, typename... Args>
inline void realloc_dispatch(WithoutInitializing, View& v, Args&&... args) {
  Kokkos::realloc(Kokkos::WithoutInitializing, v, std::forward<Args>(args)...);
}

template <class DeviceType, class Tag = Default>
void impl_testRealloc() {
  const size_t sizes[8] = {2, 3, 4, 5, 6, 7, 8, 9};

  // Check #904 fix (no reallocation if dimensions didn't change).
  {
    using view_type = Kokkos::View<int*, DeviceType>;
    view_type view_1d("view_1d", sizes[0]);
    const int* oldPointer = view_1d.data();
    auto const& oldLabel  = view_1d.label();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_1d, sizes[0]);
    auto const& newLabel = view_1d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_1d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int**, DeviceType>;
    view_type view_2d("view_2d", sizes[0], sizes[1]);
    auto const& oldLabel  = view_2d.label();
    const int* oldPointer = view_2d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_2d, sizes[0], sizes[1]);
    auto const& newLabel = view_2d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_2d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int***, DeviceType>;
    view_type view_3d("view_3d", sizes[0], sizes[1], sizes[2]);
    auto const& oldLabel  = view_3d.label();
    const int* oldPointer = view_3d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_3d, sizes[0], sizes[1], sizes[2]);
    auto const& newLabel = view_3d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_3d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int****, DeviceType>;
    view_type view_4d("view_4d", sizes[0], sizes[1], sizes[2], sizes[3]);
    auto const& oldLabel  = view_4d.label();
    const int* oldPointer = view_4d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_4d, sizes[0], sizes[1], sizes[2], sizes[3]);
    auto const& newLabel = view_4d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_4d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int*****, DeviceType>;
    view_type view_5d("view_5d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4]);
    auto const& oldLabel  = view_5d.label();
    const int* oldPointer = view_5d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_5d, sizes[0], sizes[1], sizes[2], sizes[3],
                     sizes[4]);
    auto const& newLabel = view_5d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_5d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int******, DeviceType>;
    view_type view_6d("view_6d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5]);
    const int* oldPointer = view_6d.data();
    auto const& oldLabel  = view_6d.label();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_6d, sizes[0], sizes[1], sizes[2], sizes[3],
                     sizes[4], sizes[5]);
    auto const& newLabel = view_6d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_6d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int*******, DeviceType>;
    view_type view_7d("view_7d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5], sizes[6]);
    auto const& oldLabel  = view_7d.label();
    const int* oldPointer = view_7d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_7d, sizes[0], sizes[1], sizes[2], sizes[3],
                     sizes[4], sizes[5], sizes[6]);
    auto const& newLabel = view_7d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_7d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int********, DeviceType>;
    view_type view_8d("view_8d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5], sizes[6], sizes[7]);
    auto const& oldLabel  = view_8d.label();
    const int* oldPointer = view_8d.data();
    EXPECT_NE(oldPointer, nullptr);
    realloc_dispatch(Tag{}, view_8d, sizes[0], sizes[1], sizes[2], sizes[3],
                     sizes[4], sizes[5], sizes[6], sizes[7]);
    auto const& newLabel = view_8d.label();
    EXPECT_EQ(oldLabel, newLabel);
    const int* newPointer = view_8d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
}
struct NoDefaultConstructor {
  int value;
  KOKKOS_FUNCTION
  NoDefaultConstructor(int x) : value(x) {}
};

template <class DeviceType>
void testRealloc() {
  {
    impl_testRealloc<DeviceType>();  // with data initialization
  }
  {
    impl_testRealloc<DeviceType,
                     WithoutInitializing>();  // without data initialization
  }
  // Check #6992 fix (no default initialization in realloc without initializing)
  {
    using view_type = Kokkos::View<NoDefaultConstructor*, DeviceType>;
    view_type view_1d_no_default(
        Kokkos::view_alloc(Kokkos::WithoutInitializing, "view_1d_no_default"),
        5);
    realloc_dispatch(WithoutInitializing{}, view_1d_no_default, 3);
  }
}

}  // namespace TestViewRealloc
#endif  // TESTREALLOC_HPP_
