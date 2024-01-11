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

#include <sstream>
#include <iostream>
#include <time.h>

#include <Kokkos_Core.hpp>

namespace Test {

TEST(TEST_CATEGORY, view_layoutstride_left_to_layoutleft_assignment) {
  using exec_space = TEST_EXECSPACE;

  srand(123456);  // arbitrary seed for random generator

  {  // Assignment of rank-1 LayoutLeft = LayoutStride
    int ndims   = 1;
    int dims[]  = {10};
    int order[] = {0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*, Kokkos::LayoutStride, exec_space> src("LayoutStride",
                                                                layout);

    Kokkos::View<double*, Kokkos::LayoutStride, exec_space>::HostMirror h_src =
        Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double*, Kokkos::LayoutLeft, exec_space> dst = src;

    Kokkos::View<double*, Kokkos::LayoutLeft, exec_space>::HostMirror h_dst =
        Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-2 LayoutLeft = LayoutStride
    int ndims   = 2;
    int dims[]  = {10, 9};
    int order[] = {0, 1};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double**, Kokkos::LayoutStride, exec_space> src("LayoutStride",
                                                                 layout);

    Kokkos::View<double**, Kokkos::LayoutStride, exec_space>::HostMirror h_src =
        Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double**, Kokkos::LayoutLeft, exec_space> dst = src;

    Kokkos::View<double**, Kokkos::LayoutLeft, exec_space>::HostMirror h_dst =
        Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-3 LayoutLeft = LayoutStride
    int ndims   = 3;
    int dims[]  = {10, 9, 8};
    int order[] = {0, 1, 2};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double***, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double***, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double***, Kokkos::LayoutLeft, exec_space> dst = src;

    Kokkos::View<double***, Kokkos::LayoutLeft, exec_space>::HostMirror h_dst =
        Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-4 LayoutLeft = LayoutStride
    int ndims   = 4;
    int dims[]  = {10, 9, 8, 7};
    int order[] = {0, 1, 2, 3};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double****, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double****, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double****, Kokkos::LayoutLeft, exec_space> dst = src;

    Kokkos::View<double****, Kokkos::LayoutLeft, exec_space>::HostMirror h_dst =
        Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-5 LayoutLeft = LayoutStride
    int ndims   = 5;
    int dims[]  = {10, 9, 8, 7, 6};
    int order[] = {0, 1, 2, 3, 4};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*****, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double*****, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double*****, Kokkos::LayoutLeft, exec_space> dst = src;

    Kokkos::View<double*****, Kokkos::LayoutLeft, exec_space>::HostMirror
        h_dst = Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-6 LayoutLeft = LayoutStride
    int ndims   = 6;
    int dims[]  = {10, 9, 8, 7, 6, 5};
    int order[] = {0, 1, 2, 3, 4, 5};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double******, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double******, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double******, Kokkos::LayoutLeft, exec_space> dst = src;

    Kokkos::View<double******, Kokkos::LayoutLeft, exec_space>::HostMirror
        h_dst = Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-7 LayoutLeft = LayoutStride
    int ndims   = 7;
    int dims[]  = {10, 9, 8, 7, 6, 5, 4};
    int order[] = {0, 1, 2, 3, 4, 5, 6};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*******, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double*******, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double*******, Kokkos::LayoutLeft, exec_space> dst = src;

    Kokkos::View<double*******, Kokkos::LayoutLeft, exec_space>::HostMirror
        h_dst = Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-8 LayoutLeft = LayoutStride
    int ndims   = 8;
    int dims[]  = {10, 9, 8, 7, 6, 5, 4, 3};
    int order[] = {0, 1, 2, 3, 4, 5, 6, 7};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double********, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double********, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double********, Kokkos::LayoutLeft, exec_space> dst = src;

    Kokkos::View<double********, Kokkos::LayoutLeft, exec_space>::HostMirror
        h_dst = Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
}

TEST(TEST_CATEGORY, view_layoutstride_right_to_layoutright_assignment) {
  using exec_space = TEST_EXECSPACE;

  srand(123456);  // arbitrary seed for random generator

  {  // Assignment of rank-1 LayoutRight = LayoutStride
    int ndims   = 1;
    int dims[]  = {10};
    int order[] = {0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*, Kokkos::LayoutStride, exec_space> src("LayoutStride",
                                                                layout);

    Kokkos::View<double*, Kokkos::LayoutStride, exec_space>::HostMirror h_src =
        Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double*, Kokkos::LayoutRight, exec_space> dst = src;

    Kokkos::View<double*, Kokkos::LayoutRight, exec_space>::HostMirror h_dst =
        Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-2 LayoutRight = LayoutStride
    int ndims   = 2;
    int dims[]  = {10, 9};
    int order[] = {1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double**, Kokkos::LayoutStride, exec_space> src("LayoutStride",
                                                                 layout);

    Kokkos::View<double**, Kokkos::LayoutStride, exec_space>::HostMirror h_src =
        Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double**, Kokkos::LayoutRight, exec_space> dst = src;

    Kokkos::View<double**, Kokkos::LayoutRight, exec_space>::HostMirror h_dst =
        Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-3 LayoutRight = LayoutStride
    int ndims   = 3;
    int dims[]  = {10, 9, 8};
    int order[] = {2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double***, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double***, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double***, Kokkos::LayoutRight, exec_space> dst = src;

    Kokkos::View<double***, Kokkos::LayoutRight, exec_space>::HostMirror h_dst =
        Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-4 LayoutRight = LayoutStride
    int ndims   = 4;
    int dims[]  = {10, 9, 8, 7};
    int order[] = {3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double****, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double****, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double****, Kokkos::LayoutRight, exec_space> dst = src;

    Kokkos::View<double****, Kokkos::LayoutRight, exec_space>::HostMirror
        h_dst = Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-5 LayoutRight = LayoutStride
    int ndims   = 5;
    int dims[]  = {10, 9, 8, 7, 6};
    int order[] = {4, 3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*****, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double*****, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double*****, Kokkos::LayoutRight, exec_space> dst = src;

    Kokkos::View<double*****, Kokkos::LayoutRight, exec_space>::HostMirror
        h_dst = Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-6 LayoutRight = LayoutStride
    int ndims   = 6;
    int dims[]  = {10, 9, 8, 7, 6, 5};
    int order[] = {5, 4, 3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double******, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double******, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double******, Kokkos::LayoutRight, exec_space> dst = src;

    Kokkos::View<double******, Kokkos::LayoutRight, exec_space>::HostMirror
        h_dst = Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-7 LayoutRight = LayoutStride
    int ndims   = 7;
    int dims[]  = {10, 9, 8, 7, 6, 5, 4};
    int order[] = {6, 5, 4, 3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*******, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double*******, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double*******, Kokkos::LayoutRight, exec_space> dst = src;

    Kokkos::View<double*******, Kokkos::LayoutRight, exec_space>::HostMirror
        h_dst = Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
  {  // Assignment of rank-8 LayoutRight = LayoutStride
    int ndims   = 8;
    int dims[]  = {10, 9, 8, 7, 6, 5, 4, 3};
    int order[] = {7, 6, 5, 4, 3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double********, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double********, Kokkos::LayoutStride, exec_space>::HostMirror
        h_src = Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double********, Kokkos::LayoutRight, exec_space> dst = src;

    Kokkos::View<double********, Kokkos::LayoutRight, exec_space>::HostMirror
        h_dst = Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
}

TEST(TEST_CATEGORY_DEATH, view_layoutstride_right_to_layoutleft_assignment) {
  using exec_space = TEST_EXECSPACE;

  srand(123456);  // arbitrary seed for random generator

  {  // Assignment of rank-1 LayoutLeft = LayoutStride (LayoutRight compatible)
    int ndims   = 1;
    int dims[]  = {10};
    int order[] = {0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*, Kokkos::LayoutStride, exec_space> src("LayoutStride",
                                                                layout);

    Kokkos::View<double*, Kokkos::LayoutStride, exec_space>::HostMirror h_src =
        Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double*, Kokkos::LayoutLeft, exec_space> dst;

    dst = src;

    Kokkos::View<double*, Kokkos::LayoutLeft, exec_space>::HostMirror h_dst =
        Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
// WORKAROUND OPENMPTARGET : death tests don't seem to work ...
#if defined(KOKKOS_ENABLE_OPENMPTARGET)
  return;
#endif
  {  // Assignment of rank-2 LayoutLeft = LayoutStride (LayoutRight compatible)
    int ndims   = 2;
    int dims[]  = {10, 9};
    int order[] = {1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double**, Kokkos::LayoutStride, exec_space> src("LayoutStride",
                                                                 layout);

    Kokkos::View<double**, Kokkos::LayoutLeft, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-3 LayoutLeft = LayoutStride (LayoutRight compatible)
    int ndims   = 3;
    int dims[]  = {10, 9, 8};
    int order[] = {2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double***, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double***, Kokkos::LayoutLeft, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-4 LayoutLeft = LayoutStride (LayoutRight compatible)
    int ndims   = 4;
    int dims[]  = {10, 9, 8, 7};
    int order[] = {3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double****, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double****, Kokkos::LayoutLeft, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-5 LayoutLeft = LayoutStride (LayoutRight compatible)
    int ndims   = 5;
    int dims[]  = {10, 9, 8, 7, 6};
    int order[] = {4, 3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*****, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double*****, Kokkos::LayoutLeft, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-6 LayoutLeft = LayoutStride (LayoutRight compatible)
    int ndims   = 6;
    int dims[]  = {10, 9, 8, 7, 6, 5};
    int order[] = {5, 4, 3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double******, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double******, Kokkos::LayoutLeft, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-7 LayoutLeft = LayoutStride (LayoutRight compatible)
    int ndims   = 7;
    int dims[]  = {10, 9, 8, 7, 6, 5, 4};
    int order[] = {6, 5, 4, 3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*******, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double*******, Kokkos::LayoutLeft, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-8 LayoutLeft = LayoutStride (LayoutRight compatible)
    int ndims   = 8;
    int dims[]  = {10, 9, 8, 7, 6, 5, 4, 3};
    int order[] = {7, 6, 5, 4, 3, 2, 1, 0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double********, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double********, Kokkos::LayoutLeft, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
}

TEST(TEST_CATEGORY_DEATH, view_layoutstride_left_to_layoutright_assignment) {
  using exec_space = TEST_EXECSPACE;

  srand(123456);  // arbitrary seed for random generator

  {  // Assignment of rank-1 LayoutRight = LayoutStride (LayoutLeft compatible)
    int ndims   = 1;
    int dims[]  = {10};
    int order[] = {0};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*, Kokkos::LayoutStride, exec_space> src("LayoutStride",
                                                                layout);

    Kokkos::View<double*, Kokkos::LayoutStride, exec_space>::HostMirror h_src =
        Kokkos::create_mirror_view(src);

    for (size_t i = 0; i < src.span(); i++)
      h_src.data()[i] = (double)rand() / RAND_MAX * (100);

    Kokkos::deep_copy(src, h_src);

    Kokkos::View<double*, Kokkos::LayoutRight, exec_space> dst;

    dst = src;

    Kokkos::View<double*, Kokkos::LayoutRight, exec_space>::HostMirror h_dst =
        Kokkos::create_mirror_view(dst);

    Kokkos::deep_copy(h_dst, dst);

    bool test = true;
    for (size_t i = 0; i < src.span(); i++) {
      if (h_src.data()[i] != h_dst.data()[i]) {
        test = false;
        break;
      }
    }
    ASSERT_EQ(dst.span(), src.span());
    ASSERT_EQ(test, true);
  }
// WORKAROUND OPENMPTARGET : death tests don't seem to work ...
#if defined(KOKKOS_ENABLE_OPENMPTARGET)
  return;
#endif
  {  // Assignment of rank-2 LayoutRight = LayoutStride (LayoutLeft compatible)
    int ndims   = 2;
    int dims[]  = {10, 9};
    int order[] = {0, 1};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double**, Kokkos::LayoutStride, exec_space> src("LayoutStride",
                                                                 layout);

    Kokkos::View<double**, Kokkos::LayoutRight, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-3 LayoutRight = LayoutStride (LayoutLeft compatible)
    int ndims   = 3;
    int dims[]  = {10, 9, 8};
    int order[] = {0, 1, 2};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double***, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double***, Kokkos::LayoutRight, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-4 LayoutRight = LayoutStride (LayoutLeft compatible)
    int ndims   = 4;
    int dims[]  = {10, 9, 8, 7};
    int order[] = {0, 1, 2, 3};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double****, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double****, Kokkos::LayoutRight, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-5 LayoutRight = LayoutStride (LayoutLeft compatible)
    int ndims   = 5;
    int dims[]  = {10, 9, 8, 7, 6};
    int order[] = {0, 1, 2, 3, 4};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*****, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double*****, Kokkos::LayoutRight, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-6 LayoutRight = LayoutStride (LayoutLeft compatible)
    int ndims   = 6;
    int dims[]  = {10, 9, 8, 7, 6, 5};
    int order[] = {0, 1, 2, 3, 4, 5};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double******, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double******, Kokkos::LayoutRight, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-7 LayoutRight = LayoutStride (LayoutLeft compatible)
    int ndims   = 7;
    int dims[]  = {10, 9, 8, 7, 6, 5, 4};
    int order[] = {0, 1, 2, 3, 4, 5, 6};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double*******, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double*******, Kokkos::LayoutRight, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
  {  // Assignment of rank-8 LayoutRight = LayoutStride (LayoutLeft compatible)
    int ndims   = 8;
    int dims[]  = {10, 9, 8, 7, 6, 5, 4, 3};
    int order[] = {0, 1, 2, 3, 4, 5, 6, 7};
    Kokkos::LayoutStride layout =
        Kokkos::LayoutStride::order_dimensions(ndims, order, dims);
    Kokkos::View<double********, Kokkos::LayoutStride, exec_space> src(
        "LayoutStride", layout);

    Kokkos::View<double********, Kokkos::LayoutRight, exec_space> dst;

    ::testing::FLAGS_gtest_death_test_style = "threadsafe";
    ASSERT_DEATH({ dst = src; },
                 "View assignment must have compatible layouts");
  }
}

}  // namespace Test

#include <TestIrregularLayout.hpp>
