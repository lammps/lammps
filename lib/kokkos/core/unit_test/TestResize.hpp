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
#ifndef TESTRESIZE_HPP_
#define TESTRESIZE_HPP_

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace TestViewResize {

struct Default {};
struct WithoutInitializing {};

template <typename View, typename... Args>
inline void resize_dispatch(Default, View& v, Args&&... args) {
  Kokkos::resize(v, std::forward<Args>(args)...);
}

template <typename View, typename... Args>
inline void resize_dispatch(WithoutInitializing, View& v, Args&&... args) {
  Kokkos::resize(Kokkos::WithoutInitializing, v, std::forward<Args>(args)...);
}

template <class DeviceType, class Tag = Default>
void impl_testResize() {
  const size_t sizes[8] = {2, 3, 4, 5, 6, 7, 8, 9};

  // Check #904 fix (no reallocation if dimensions didn't change).
  {
    using view_type = Kokkos::View<int*, DeviceType>;
    view_type view_1d("view_1d", sizes[0]);
    const int* oldPointer = view_1d.data();
    EXPECT_NE(oldPointer, nullptr);
    resize_dispatch(Tag{}, view_1d, sizes[0]);
    const int* newPointer = view_1d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int**, DeviceType>;
    view_type view_2d("view_2d", sizes[0], sizes[1]);
    const int* oldPointer = view_2d.data();
    EXPECT_NE(oldPointer, nullptr);
    resize_dispatch(Tag{}, view_2d, sizes[0], sizes[1]);
    const int* newPointer = view_2d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int***, DeviceType>;
    view_type view_3d("view_3d", sizes[0], sizes[1], sizes[2]);
    const int* oldPointer = view_3d.data();
    EXPECT_NE(oldPointer, nullptr);
    resize_dispatch(Tag{}, view_3d, sizes[0], sizes[1], sizes[2]);
    const int* newPointer = view_3d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int****, DeviceType>;
    view_type view_4d("view_4d", sizes[0], sizes[1], sizes[2], sizes[3]);
    const int* oldPointer = view_4d.data();
    EXPECT_NE(oldPointer, nullptr);
    resize_dispatch(Tag{}, view_4d, sizes[0], sizes[1], sizes[2], sizes[3]);
    const int* newPointer = view_4d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int*****, DeviceType>;
    view_type view_5d("view_5d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4]);
    const int* oldPointer = view_5d.data();
    EXPECT_NE(oldPointer, nullptr);
    resize_dispatch(Tag{}, view_5d, sizes[0], sizes[1], sizes[2], sizes[3],
                    sizes[4]);
    const int* newPointer = view_5d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int******, DeviceType>;
    view_type view_6d("view_6d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5]);
    const int* oldPointer = view_6d.data();
    EXPECT_NE(oldPointer, nullptr);
    resize_dispatch(Tag{}, view_6d, sizes[0], sizes[1], sizes[2], sizes[3],
                    sizes[4], sizes[5]);
    const int* newPointer = view_6d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int*******, DeviceType>;
    view_type view_7d("view_7d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5], sizes[6]);
    const int* oldPointer = view_7d.data();
    EXPECT_NE(oldPointer, nullptr);
    resize_dispatch(Tag{}, view_7d, sizes[0], sizes[1], sizes[2], sizes[3],
                    sizes[4], sizes[5], sizes[6]);
    const int* newPointer = view_7d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  {
    using view_type = Kokkos::View<int********, DeviceType>;
    view_type view_8d("view_8d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5], sizes[6], sizes[7]);
    const int* oldPointer = view_8d.data();
    EXPECT_NE(oldPointer, nullptr);
    resize_dispatch(Tag{}, view_8d, sizes[0], sizes[1], sizes[2], sizes[3],
                    sizes[4], sizes[5], sizes[6], sizes[7]);
    const int* newPointer = view_8d.data();
    EXPECT_EQ(oldPointer, newPointer);
  }
  // Resize without initialization: check if data preserved
  {
    using view_type = Kokkos::View<int*, DeviceType>;
    view_type view_1d("view_1d", sizes[0]);
    typename view_type::HostMirror h_view_1d_old =
        Kokkos::create_mirror(view_1d);
    Kokkos::deep_copy(view_1d, 111);
    Kokkos::deep_copy(h_view_1d_old, view_1d);
    resize_dispatch(Tag{}, view_1d, 2 * sizes[0]);
    EXPECT_EQ(view_1d.extent(0), 2 * sizes[0]);
    typename view_type::HostMirror h_view_1d =
        Kokkos::create_mirror_view(view_1d);
    Kokkos::deep_copy(h_view_1d, view_1d);
    bool test = true;
    for (size_t i0 = 0; i0 < sizes[0]; ++i0) {
      if (h_view_1d(i0) != h_view_1d_old(i0)) {
        test = false;
        break;
      }
    }
    EXPECT_TRUE(test);
  }
  {
    using view_type = Kokkos::View<int**, DeviceType>;
    view_type view_2d("view_2d", sizes[0], sizes[1]);
    typename view_type::HostMirror h_view_2d_old =
        Kokkos::create_mirror(view_2d);
    Kokkos::deep_copy(view_2d, 222);
    Kokkos::deep_copy(h_view_2d_old, view_2d);
    resize_dispatch(Tag{}, view_2d, 2 * sizes[0], sizes[1]);
    EXPECT_EQ(view_2d.extent(0), 2 * sizes[0]);
    typename view_type::HostMirror h_view_2d =
        Kokkos::create_mirror_view(view_2d);
    Kokkos::deep_copy(h_view_2d, view_2d);
    bool test = true;
    for (size_t i0 = 0; i0 < sizes[0]; ++i0) {
      for (size_t i1 = 0; i1 < sizes[1]; ++i1) {
        if (h_view_2d(i0, i1) != h_view_2d_old(i0, i1)) {
          test = false;
          break;
        }
      }
    }
    EXPECT_TRUE(test);
  }
  {
    using view_type = Kokkos::View<int***, DeviceType>;
    view_type view_3d("view_3d", sizes[0], sizes[1], sizes[2]);
    typename view_type::HostMirror h_view_3d_old =
        Kokkos::create_mirror(view_3d);
    Kokkos::deep_copy(view_3d, 333);
    Kokkos::deep_copy(h_view_3d_old, view_3d);
    resize_dispatch(Tag{}, view_3d, 2 * sizes[0], sizes[1], sizes[2]);
    EXPECT_EQ(view_3d.extent(0), 2 * sizes[0]);
    typename view_type::HostMirror h_view_3d =
        Kokkos::create_mirror_view(view_3d);
    Kokkos::deep_copy(h_view_3d, view_3d);
    bool test = true;
    for (size_t i0 = 0; i0 < sizes[0]; ++i0) {
      for (size_t i1 = 0; i1 < sizes[1]; ++i1) {
        for (size_t i2 = 0; i2 < sizes[2]; ++i2) {
          if (h_view_3d(i0, i1, i2) != h_view_3d_old(i0, i1, i2)) {
            test = false;
            break;
          }
        }
      }
    }
    EXPECT_TRUE(test);
  }
  {
    using view_type = Kokkos::View<int****, DeviceType>;
    view_type view_4d("view_4d", sizes[0], sizes[1], sizes[2], sizes[3]);
    typename view_type::HostMirror h_view_4d_old =
        Kokkos::create_mirror(view_4d);
    Kokkos::deep_copy(view_4d, 444);
    Kokkos::deep_copy(h_view_4d_old, view_4d);
    resize_dispatch(Tag{}, view_4d, 2 * sizes[0], sizes[1], sizes[2], sizes[3]);
    EXPECT_EQ(view_4d.extent(0), 2 * sizes[0]);
    typename view_type::HostMirror h_view_4d =
        Kokkos::create_mirror_view(view_4d);
    Kokkos::deep_copy(h_view_4d, view_4d);
    bool test = true;
    for (size_t i0 = 0; i0 < sizes[0]; ++i0) {
      for (size_t i1 = 0; i1 < sizes[1]; ++i1) {
        for (size_t i2 = 0; i2 < sizes[2]; ++i2) {
          for (size_t i3 = 0; i3 < sizes[3]; ++i3) {
            if (h_view_4d(i0, i1, i2, i3) != h_view_4d_old(i0, i1, i2, i3)) {
              test = false;
              break;
            }
          }
        }
      }
    }
    EXPECT_TRUE(test);
  }
  {
    using view_type = Kokkos::View<int*****, DeviceType>;
    view_type view_5d("view_5d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4]);
    typename view_type::HostMirror h_view_5d_old =
        Kokkos::create_mirror(view_5d);
    Kokkos::deep_copy(view_5d, 555);
    Kokkos::deep_copy(h_view_5d_old, view_5d);
    resize_dispatch(Tag{}, view_5d, 2 * sizes[0], sizes[1], sizes[2], sizes[3],
                    sizes[4]);
    EXPECT_EQ(view_5d.extent(0), 2 * sizes[0]);
    typename view_type::HostMirror h_view_5d =
        Kokkos::create_mirror_view(view_5d);
    Kokkos::deep_copy(h_view_5d, view_5d);
    bool test = true;
    for (size_t i0 = 0; i0 < sizes[0]; ++i0) {
      for (size_t i1 = 0; i1 < sizes[1]; ++i1) {
        for (size_t i2 = 0; i2 < sizes[2]; ++i2) {
          for (size_t i3 = 0; i3 < sizes[3]; ++i3) {
            for (size_t i4 = 0; i4 < sizes[4]; ++i4) {
              if (h_view_5d(i0, i1, i2, i3, i4) !=
                  h_view_5d_old(i0, i1, i2, i3, i4)) {
                test = false;
                break;
              }
            }
          }
        }
      }
    }
    EXPECT_TRUE(test);
  }
  {
    using view_type = Kokkos::View<int******, DeviceType>;
    view_type view_6d("view_6d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5]);
    typename view_type::HostMirror h_view_6d_old =
        Kokkos::create_mirror(view_6d);
    Kokkos::deep_copy(view_6d, 666);
    Kokkos::deep_copy(h_view_6d_old, view_6d);
    resize_dispatch(Tag{}, view_6d, 2 * sizes[0], sizes[1], sizes[2], sizes[3],
                    sizes[4], sizes[5]);
    EXPECT_EQ(view_6d.extent(0), 2 * sizes[0]);
    typename view_type::HostMirror h_view_6d =
        Kokkos::create_mirror_view(view_6d);
    Kokkos::deep_copy(h_view_6d, view_6d);
    bool test = true;
    for (size_t i0 = 0; i0 < sizes[0]; ++i0) {
      for (size_t i1 = 0; i1 < sizes[1]; ++i1) {
        for (size_t i2 = 0; i2 < sizes[2]; ++i2) {
          for (size_t i3 = 0; i3 < sizes[3]; ++i3) {
            for (size_t i4 = 0; i4 < sizes[4]; ++i4) {
              for (size_t i5 = 0; i5 < sizes[5]; ++i5) {
                if (h_view_6d(i0, i1, i2, i3, i4, i5) !=
                    h_view_6d_old(i0, i1, i2, i3, i4, i5)) {
                  test = false;
                  break;
                }
              }
            }
          }
        }
      }
    }
    EXPECT_TRUE(test);
  }
  {
    using view_type = Kokkos::View<int*******, DeviceType>;
    view_type view_7d("view_7d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5], sizes[6]);
    typename view_type::HostMirror h_view_7d_old =
        Kokkos::create_mirror(view_7d);
    Kokkos::deep_copy(view_7d, 777);
    Kokkos::deep_copy(h_view_7d_old, view_7d);
    resize_dispatch(Tag{}, view_7d, 2 * sizes[0], sizes[1], sizes[2], sizes[3],
                    sizes[4], sizes[5], sizes[6]);
    EXPECT_EQ(view_7d.extent(0), 2 * sizes[0]);
    typename view_type::HostMirror h_view_7d =
        Kokkos::create_mirror_view(view_7d);
    Kokkos::deep_copy(h_view_7d, view_7d);
    bool test = true;
    for (size_t i0 = 0; i0 < sizes[0]; ++i0) {
      for (size_t i1 = 0; i1 < sizes[1]; ++i1) {
        for (size_t i2 = 0; i2 < sizes[2]; ++i2) {
          for (size_t i3 = 0; i3 < sizes[3]; ++i3) {
            for (size_t i4 = 0; i4 < sizes[4]; ++i4) {
              for (size_t i5 = 0; i5 < sizes[5]; ++i5) {
                for (size_t i6 = 0; i6 < sizes[6]; ++i6) {
                  if (h_view_7d(i0, i1, i2, i3, i4, i5, i6) !=
                      h_view_7d_old(i0, i1, i2, i3, i4, i5, i6)) {
                    test = false;
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
    EXPECT_TRUE(test);
  }
  {
    using view_type = Kokkos::View<int********, DeviceType>;
    view_type view_8d("view_8d", sizes[0], sizes[1], sizes[2], sizes[3],
                      sizes[4], sizes[5], sizes[6], sizes[7]);
    typename view_type::HostMirror h_view_8d_old =
        Kokkos::create_mirror(view_8d);
    Kokkos::deep_copy(view_8d, 888);
    Kokkos::deep_copy(h_view_8d_old, view_8d);
    resize_dispatch(Tag{}, view_8d, 2 * sizes[0], sizes[1], sizes[2], sizes[3],
                    sizes[4], sizes[5], sizes[6], sizes[7]);
    EXPECT_EQ(view_8d.extent(0), 2 * sizes[0]);
    typename view_type::HostMirror h_view_8d =
        Kokkos::create_mirror_view(view_8d);
    Kokkos::deep_copy(h_view_8d, view_8d);
    bool test = true;
    for (size_t i0 = 0; i0 < sizes[0]; ++i0) {
      for (size_t i1 = 0; i1 < sizes[1]; ++i1) {
        for (size_t i2 = 0; i2 < sizes[2]; ++i2) {
          for (size_t i3 = 0; i3 < sizes[3]; ++i3) {
            for (size_t i4 = 0; i4 < sizes[4]; ++i4) {
              for (size_t i5 = 0; i5 < sizes[5]; ++i5) {
                for (size_t i6 = 0; i6 < sizes[6]; ++i6) {
                  for (size_t i7 = 0; i7 < sizes[7]; ++i7) {
                    if (h_view_8d(i0, i1, i2, i3, i4, i5, i6, i7) !=
                        h_view_8d_old(i0, i1, i2, i3, i4, i5, i6, i7)) {
                      test = false;
                      break;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    EXPECT_TRUE(test);
  }
}

template <class DeviceType>
void testResize() {
  {
    impl_testResize<DeviceType>();  // with data initialization
  }
  {
    impl_testResize<DeviceType,
                    WithoutInitializing>();  // without data initialization
  }
}

}  // namespace TestViewResize
#endif  // TESTRESIZE_HPP_
