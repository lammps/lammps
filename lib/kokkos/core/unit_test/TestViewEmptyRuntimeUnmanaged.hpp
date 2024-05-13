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

#include <Kokkos_Core.hpp>

namespace {

template <class T>
void test_empty_view_runtime_unmanaged() {
  T d{};
  auto* p = reinterpret_cast<T*>(0xABADBABE);

  (void)Kokkos::View<T>(p);
  (void)Kokkos::View<T>(&d);
  (void)Kokkos::View<T>(nullptr);
  (void)Kokkos::View<T>(NULL);  // NOLINT(modernize-use-nullptr)
  (void)Kokkos::View<T>(0);     // NOLINT(modernize-use-nullptr)

  (void)Kokkos::View<T*>(p, 0);
  (void)Kokkos::View<T*>(&d, 0);
  (void)Kokkos::View<T*>(nullptr, 0);
  (void)Kokkos::View<T*>(NULL, 0);  // NOLINT(modernize-use-nullptr)
  (void)Kokkos::View<T*>(0, 0);     // NOLINT(modernize-use-nullptr)

  (void)Kokkos::View<T**>(p, 0, 0);
  (void)Kokkos::View<T**>(&d, 0, 0);
  (void)Kokkos::View<T**>(nullptr, 0, 0);
  (void)Kokkos::View<T**>(NULL, 0, 0);  // NOLINT(modernize-use-nullptr)
  (void)Kokkos::View<T**>(0, 0, 0);     // NOLINT(modernize-use-nullptr)
}

TEST(TEST_CATEGORY, view_empty_runtime_unmanaged) {
  test_empty_view_runtime_unmanaged<float>();
  test_empty_view_runtime_unmanaged<const double>();
  test_empty_view_runtime_unmanaged<int>();
  test_empty_view_runtime_unmanaged<char>();
  test_empty_view_runtime_unmanaged<const char>();
}

}  // namespace
