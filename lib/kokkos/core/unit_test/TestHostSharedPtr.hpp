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

#include <impl/Kokkos_HostSharedPtr.hpp>

#include <gtest/gtest.h>

using Kokkos::Impl::HostSharedPtr;

TEST(TEST_CATEGORY, host_shared_ptr_use_count) {
  using T = int;
  {
    HostSharedPtr<T> p1;
    EXPECT_EQ(p1.use_count(), 0);
  }
  {
    HostSharedPtr<T> p1(nullptr);
    EXPECT_EQ(p1.use_count(), 0);
  }
  {
    HostSharedPtr<T> p1(new T());
    EXPECT_EQ(p1.use_count(), 1);
  }
  {
    HostSharedPtr<T> p1(new T(), [](T* p) { delete p; });
    EXPECT_EQ(p1.use_count(), 1);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    EXPECT_EQ(p1.use_count(), 1);
  }
  {
    HostSharedPtr<T> p1(new T());
    HostSharedPtr<T> p2(p1);  // copy construction
    EXPECT_EQ(p1.use_count(), 2);
    EXPECT_EQ(p2.use_count(), 2);
  }
  {
    HostSharedPtr<T> p1(new T());
    HostSharedPtr<T> p2(std::move(p1));  // move construction
    EXPECT_EQ(p2.use_count(), 1);
  }
  {
    HostSharedPtr<T> p1(new T());
    HostSharedPtr<T> p2;
    p2 = p1;  // copy assignment
    EXPECT_EQ(p1.use_count(), 2);
    EXPECT_EQ(p2.use_count(), 2);
  }
  {
    HostSharedPtr<T> p1(new T());
    HostSharedPtr<T> p2;
    p2 = std::move(p1);  // move assignment
    EXPECT_EQ(p2.use_count(), 1);
  }
}

TEST(TEST_CATEGORY, host_shared_ptr_get) {
  using T = int;
  {
    HostSharedPtr<T> p1;
    EXPECT_EQ(p1.get(), nullptr);
  }
  {
    HostSharedPtr<T> p1(nullptr);
    EXPECT_EQ(p1.get(), nullptr);
  }
  {
    T* p_i = new T();
    HostSharedPtr<T> p1(p_i);
    EXPECT_EQ(p1.get(), p_i);
  }
  {
    T* p_i = new T();
    HostSharedPtr<T> p1(p_i, [](T* p) { delete p; });
    EXPECT_EQ(p1.get(), p_i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    EXPECT_EQ(p1.get(), &i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    HostSharedPtr<T> p2(p1);  // copy construction
    EXPECT_EQ(p1.get(), &i);
    EXPECT_EQ(p1.get(), &i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    HostSharedPtr<T> p2(std::move(p1));  // move construction
    EXPECT_EQ(p1.get(), nullptr);
    EXPECT_EQ(p2.get(), &i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    HostSharedPtr<T> p2;
    p2 = p1;  // copy assignment
    EXPECT_EQ(p1.get(), &i);
    EXPECT_EQ(p2.get(), &i);
  }
  {
    T i;
    HostSharedPtr<T> p1(&i, [](T*) {});
    HostSharedPtr<T> p2;
    p2 = std::move(p1);  // move assignment
    EXPECT_EQ(p1.get(), nullptr);
    EXPECT_EQ(p2.get(), &i);
  }
}
