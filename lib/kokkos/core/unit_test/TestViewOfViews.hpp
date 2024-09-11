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

// User-defined type with a View data member
template <class V>
class S {
  V v_;

 public:
  template <class... Extents>
  S(std::string label, Extents... extents) : v_(std::move(label), extents...) {}
  S() = default;
};

template <class V>
void test_view_of_views() {
  using VoV = Kokkos::View<V**, Kokkos::HostSpace>;
  {  // assigning a default-constructed view to destruct the inner objects
    VoV vov("vov", 2, 3);
    V a("a");
    V b("b");
    vov(0, 0) = a;
    vov(1, 0) = a;
    vov(0, 1) = b;
#ifndef KOKKOS_ENABLE_IMPL_VIEW_OF_VIEWS_DESTRUCTOR_PRECONDITION_VIOLATION_WORKAROUND
    vov(0, 0) = V();
    vov(1, 0) = V();
    vov(0, 1) = V();
#endif
  }
  {  // using placement new to construct the inner objects and explicitly
     // calling the destructor
    VoV vov(Kokkos::view_alloc("vov", Kokkos::WithoutInitializing), 2, 3);
    V a("a");
    V b("b");
    new (&vov(0, 0)) V(a);
    new (&vov(1, 0)) V(a);
    new (&vov(0, 1)) V(b);
#ifndef KOKKOS_ENABLE_IMPL_VIEW_OF_VIEWS_DESTRUCTOR_PRECONDITION_VIOLATION_WORKAROUND
    vov(0, 0).~V();
    vov(1, 0).~V();
    vov(0, 1).~V();
#else
    // leaks memory
#endif
  }
}

TEST(TEST_CATEGORY, view_of_views) {
  test_view_of_views<Kokkos::View<int, TEST_EXECSPACE>>();
  test_view_of_views<Kokkos::View<int[4], TEST_EXECSPACE>>();
  // User-defined type with View data member
  test_view_of_views<S<Kokkos::View<float, TEST_EXECSPACE>>>();
}

}  // namespace
