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

// User-defined types with a View data member
template <class V>
class S {
  V v_;

 public:
  template <class... Extents>
  S(std::string label, Extents... extents) : v_(std::move(label), extents...) {}
  KOKKOS_DEFAULTED_FUNCTION S() = default;
};

template <class V>
class N {  // not default constructible
  V v_;

 public:
  template <class... Extents>
  N(std::string label, Extents... extents) : v_(std::move(label), extents...) {}
};

template <class V>
class H {  // constructible and destructible only from on the host side
  V v_;

 public:
  template <class... Extents>
  H(std::string label, Extents... extents) : v_(std::move(label), extents...) {}
  H() {}
  ~H() {}
};

template <class V>
void test_view_of_views_default() {
  // assigning a default-constructed view to destruct the inner objects
  using VoV = Kokkos::View<V**, Kokkos::HostSpace>;
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

template <class V>
void test_view_of_views_without_initializing() {
  // using placement new to construct the inner objects and explicitly
  // calling the destructor
  using VoV = Kokkos::View<V**, Kokkos::HostSpace>;
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

template <class V>
void test_view_of_views_sequential_host_init() {
  // inner views value-initialized sequentially on the host, and also
  // sequentially destructed on the host, without the need to cleanup
  using VoV = Kokkos::View<V**, Kokkos::HostSpace>;
  VoV vov(Kokkos::view_alloc("vov", Kokkos::SequentialHostInit), 2, 3);
  V a("a");
  V b("b");
  vov(0, 0) = a;
  vov(1, 0) = a;
  vov(0, 1) = b;
}

TEST(TEST_CATEGORY, view_of_views_default) {
  test_view_of_views_default<Kokkos::View<int, TEST_EXECSPACE>>();
  test_view_of_views_default<Kokkos::View<int[4], TEST_EXECSPACE>>();
  // User-defined type with View data member
  test_view_of_views_default<S<Kokkos::View<float, TEST_EXECSPACE>>>();
}

TEST(TEST_CATEGORY, view_of_views_without_initializing) {
  test_view_of_views_without_initializing<Kokkos::View<int, TEST_EXECSPACE>>();
  test_view_of_views_without_initializing<
      S<Kokkos::View<float, TEST_EXECSPACE>>>();
  test_view_of_views_without_initializing<
      N<Kokkos::View<double, TEST_EXECSPACE>>>();
  test_view_of_views_without_initializing<
      H<Kokkos::View<int, TEST_EXECSPACE>>>();
}

TEST(TEST_CATEGORY, test_view_of_views_sequential_host_init) {
  test_view_of_views_sequential_host_init<Kokkos::View<int, TEST_EXECSPACE>>();
  test_view_of_views_sequential_host_init<
      S<Kokkos::View<float, TEST_EXECSPACE>>>();
  test_view_of_views_sequential_host_init<
      H<Kokkos::View<int, TEST_EXECSPACE>>>();
}

}  // namespace
