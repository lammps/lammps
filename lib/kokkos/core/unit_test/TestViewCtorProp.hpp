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

using vcp_empty_t      = Kokkos::Impl::ViewCtorProp<>;
using vcp_label_base_t = Kokkos::Impl::ViewCtorProp<void, std::string>;
using vcp_label_t      = Kokkos::Impl::ViewCtorProp<std::string>;

// Check traits of Kokkos::Impl::ViewCtorProp<>.
TEST(TEST_CATEGORY, vcp_empty_traits) {
  // Check that the empty view constructor properties class is default
  // constructible. This is needed for calls of Kokkos::view_alloc().
  static_assert(std::is_default_constructible_v<vcp_empty_t>);
  static_assert(std::is_same_v<decltype(Kokkos::view_alloc()), vcp_empty_t>);
}

// Check Kokkos::Impl::is_view_label.
TEST(TEST_CATEGORY, is_view_label) {
  static_assert(Kokkos::Impl::is_view_label<std::string>::value);

  static_assert(Kokkos::Impl::is_view_label<const char[3]>::value);
  static_assert(Kokkos::Impl::is_view_label<char[3]>::value);

  // A char* is not a label. Thus, a label is distinguished from a pointer type.
  static_assert(!Kokkos::Impl::is_view_label<char*>::value);
}

// Check traits of base class Kokkos::Impl::ViewCtorProp<void, std::string>.
TEST(TEST_CATEGORY, vcp_label_base_traits) {
  static_assert(std::is_same_v<typename vcp_label_base_t::type, std::string>);

  // Check that the base class is default constructible. The default constructor
  // may be called by the copy constructor of derived classes, such as when
  // copy constructing a view constructor properties object from another view
  // constructor properties object that holds fewer properties.
  static_assert(std::is_default_constructible_v<vcp_label_base_t>);
}

// Check traits of derived class Kokkos::Impl::ViewCtorProp<std::string>.
TEST(TEST_CATEGORY, vcp_label_traits) {
  static_assert(std::is_base_of_v<vcp_label_base_t, vcp_label_t>);

  static_assert(vcp_label_t::has_label);

  // Check that the derived class is not default constructible. It is a design
  // choice to not allow the default constructor to be called.
  static_assert(!std::is_default_constructible_v<vcp_label_t>);
}

// Check that Kokkos::view_alloc perfect forwards a label passed by
// rvalue reference, and check that the constructor
// of Kokkos::Impl::ViewCtorProp<std::string> moves this label.
TEST(TEST_CATEGORY, view_alloc_can_perfect_forward_label) {
  std::string label("our label");

  auto prop = Kokkos::view_alloc(std::move(label));

  // This is not actually guaranteed in the C++ standard
  // https://eel.is/c++draft/basic.string#string.cons-24
  // > left in a valid but unspecified state
  ASSERT_TRUE(label.empty());  // NOLINT(bugprone-use-after-move)
  ASSERT_EQ(Kokkos::Impl::get_property<Kokkos::Impl::LabelTag>(prop),
            "our label");
}

// Check the copy constructor of Kokkos::Impl::ViewCtorProp<std::string>.
TEST(TEST_CATEGORY, vcp_label_copy_constructor) {
  // Copy construction from a view constructor properties object with a label.
  static_assert(std::is_copy_constructible_v<vcp_label_t>);

  vcp_label_t prop = Kokkos::view_alloc("our label");
  vcp_label_t prop_copy(prop);

  ASSERT_EQ(Kokkos::Impl::get_property<Kokkos::Impl::LabelTag>(prop),
            "our label");
  ASSERT_EQ(Kokkos::Impl::get_property<Kokkos::Impl::LabelTag>(prop_copy),
            "our label");
}

TEST(TEST_CATEGORY, vcp_pointer_add_property) {
  double dummy        = 1.;
  auto properties     = Kokkos::view_wrap(&dummy);
  auto new_properties = Kokkos::Impl::with_properties_if_unset(
      properties, std::string("our label"));
  ASSERT_EQ(Kokkos::Impl::get_property<Kokkos::Impl::LabelTag>(new_properties),
            "our label");
}

}  // namespace
