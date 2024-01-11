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
#ifndef TESTVIEWHOOKS_HPP_
#define TESTVIEWHOOKS_HPP_

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

namespace Test {
template <class DeviceType>
struct TestViewHooks {
  struct TestSubscriber;

  static_assert(
      Kokkos::Experimental::is_hooks_policy<
          Kokkos::Experimental::SubscribableViewHooks<TestSubscriber> >::value,
      "Must be a hooks policy");

  using test_view_type =
      Kokkos::View<double **,
                   Kokkos::Experimental::SubscribableViewHooks<TestSubscriber>,
                   DeviceType>;

  struct TestSubscriber {
    static test_view_type *self_ptr;
    static const test_view_type *other_ptr;

    template <typename View>
    static void copy_constructed(View &self, const View &other) {
      self_ptr  = &self;
      other_ptr = &other;
    }

    template <typename View>
    static void move_constructed(View &self, const View &other) {
      self_ptr  = &self;
      other_ptr = &other;
    }

    template <typename View>
    static void copy_assigned(View &self, const View &other) {
      self_ptr  = &self;
      other_ptr = &other;
    }

    template <typename View>
    static void move_assigned(View &self, const View &other) {
      self_ptr  = &self;
      other_ptr = &other;
    }

    static void reset() {
      self_ptr  = nullptr;
      other_ptr = nullptr;
    }
  };

  static void testViewHooksCopyConstruct() {
    TestSubscriber::reset();
    test_view_type testa;

    test_view_type testb(testa);
    EXPECT_EQ(TestSubscriber::self_ptr, &testb);
    EXPECT_EQ(TestSubscriber::other_ptr, &testa);
  }

  static void testViewHooksMoveConstruct() {
    TestSubscriber::reset();
    test_view_type testa;

    test_view_type testb(std::move(testa));
    EXPECT_EQ(TestSubscriber::self_ptr, &testb);

    // This is valid, even if the view is moved-from
    EXPECT_EQ(TestSubscriber::other_ptr, &testa);
  }

  static void testViewHooksCopyAssign() {
    TestSubscriber::reset();
    test_view_type testa;

    test_view_type testb;
    testb = testa;
    EXPECT_EQ(TestSubscriber::self_ptr, &testb);
    EXPECT_EQ(TestSubscriber::other_ptr, &testa);
  }

  static void testViewHooksMoveAssign() {
    TestSubscriber::reset();
    test_view_type testa;

    test_view_type testb;
    testb = std::move(testa);
    EXPECT_EQ(TestSubscriber::self_ptr, &testb);

    // This is valid, even if the view is moved-from
    EXPECT_EQ(TestSubscriber::other_ptr, &testa);
  }
};

template <class DeviceType>
typename TestViewHooks<DeviceType>::test_view_type
    *TestViewHooks<DeviceType>::TestSubscriber::self_ptr = nullptr;

template <class DeviceType>
const typename TestViewHooks<DeviceType>::test_view_type
    *TestViewHooks<DeviceType>::TestSubscriber::other_ptr = nullptr;

TEST(TEST_CATEGORY, view_hooks) {
  using ExecSpace = TEST_EXECSPACE;
  TestViewHooks<ExecSpace>::testViewHooksCopyConstruct();
  TestViewHooks<ExecSpace>::testViewHooksMoveConstruct();
  TestViewHooks<ExecSpace>::testViewHooksCopyAssign();
  TestViewHooks<ExecSpace>::testViewHooksMoveAssign();
}

}  // namespace Test
#endif  // TESTVIEWHOOKS_HPP_
