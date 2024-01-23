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
#include <sstream>
#include <iostream>

namespace Test {

template <class Space>
struct NestedView {
  Kokkos::View<int *, Space> member;

 public:
  KOKKOS_INLINE_FUNCTION
  NestedView() : member() {}

  KOKKOS_INLINE_FUNCTION
  NestedView &operator=(const Kokkos::View<int *, Space> &lhs) {
    member = lhs;
    if (member.extent(0)) Kokkos::atomic_add(&member(0), 1);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  ~NestedView() {
    if (member.extent(0)) {
      Kokkos::atomic_add(&member(0), -1);
    }
  }
};

template <class Space>
struct NestedViewFunctor {
  Kokkos::View<NestedView<Space> *, Space> nested;
  Kokkos::View<int *, Space> array;

  NestedViewFunctor(const Kokkos::View<NestedView<Space> *, Space> &arg_nested,
                    const Kokkos::View<int *, Space> &arg_array)
      : nested(arg_nested), array(arg_array) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { nested[i] = array; }
};

template <class Space>
void view_nested_view() {
  Kokkos::View<int *, Space> tracking("tracking", 1);

  typename Kokkos::View<int *, Space>::HostMirror host_tracking =
      Kokkos::create_mirror(tracking);

  {
    Kokkos::View<NestedView<Space> *, Space> a("a_nested_view", 2);

    Kokkos::parallel_for(Kokkos::RangePolicy<Space>(0, 2),
                         NestedViewFunctor<Space>(a, tracking));
    Kokkos::deep_copy(host_tracking, tracking);
    ASSERT_EQ(2, host_tracking(0));

    Kokkos::View<NestedView<Space> *, Space> b("b_nested_view", 2);
    Kokkos::parallel_for(Kokkos::RangePolicy<Space>(0, 2),
                         NestedViewFunctor<Space>(b, tracking));
    Kokkos::deep_copy(host_tracking, tracking);
    ASSERT_EQ(4, host_tracking(0));
  }

  Kokkos::deep_copy(host_tracking, tracking);

  ASSERT_EQ(0, host_tracking(0));
}

TEST(TEST_CATEGORY, view_nested_view) { view_nested_view<TEST_EXECSPACE>(); }

}  // namespace Test
