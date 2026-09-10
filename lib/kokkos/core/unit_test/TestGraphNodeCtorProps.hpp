// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_GraphNodeCtorProps.hpp>

#include <gtest/gtest.h>

namespace {

constexpr bool test_node_ctor_prop_traits() {
  static_assert(std::is_constructible_v<Kokkos::Impl::NodeCtorProp<std::string>,
                                        const char(&)[6]>);
  return true;
}
static_assert(test_node_ctor_prop_traits());

TEST(TEST_CATEGORY, node_props_empty) {
  using props_t = decltype(Kokkos::Experimental::node_props());

  static_assert(std::same_as<props_t, Kokkos::Impl::NodeCtorProps<>>);

  static_assert(Kokkos::Impl::is_node_props_v<props_t>);

  static_assert(props_t::has<std::string> == false);
  static_assert(props_t::has<TEST_EXECSPACE> == false);
}

TEST(TEST_CATEGORY, node_props_label) {
  using props_t = decltype(Kokkos::Experimental::node_props("label"));

  static_assert(
      std::same_as<props_t, Kokkos::Impl::NodeCtorProps<std::string>>);

  static_assert(Kokkos::Impl::is_node_props_v<props_t>);

  static_assert(props_t::has<std::string> == true);
  static_assert(props_t::has<TEST_EXECSPACE> == false);

  const auto props = Kokkos::Experimental::node_props("label");

  ASSERT_EQ(Kokkos::Impl::get_property<std::string>(props), "label");
}

TEST(TEST_CATEGORY, node_props_device_handle) {
  using props_t = decltype(Kokkos::Experimental::node_props(
      Kokkos::Experimental::get_device_handle(std::declval<TEST_EXECSPACE>())));

  static_assert(
      std::same_as<props_t, Kokkos::Impl::NodeCtorProps<
                                Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>>);

  static_assert(Kokkos::Impl::is_node_props_v<props_t>);

  static_assert(props_t::has<std::string> == false);
  static_assert(props_t::has<Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>> ==
                true);

  const auto props = Kokkos::Experimental::node_props(
      Kokkos::Experimental::get_device_handle(TEST_EXECSPACE{}));

  ASSERT_EQ(
      Kokkos::Impl::get_property<Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>(
          props),
      Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>{});
}

TEST(TEST_CATEGORY, node_props_device_handle_label) {
  using props_t = decltype(Kokkos::Experimental::node_props(
      Kokkos::Experimental::get_device_handle(std::declval<TEST_EXECSPACE>()),
      "label"));

  static_assert(
      std::same_as<props_t, Kokkos::Impl::NodeCtorProps<
                                Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>,
                                std::string>>);

  static_assert(Kokkos::Impl::is_node_props_v<props_t>);

  static_assert(props_t::has<std::string> == true);
  static_assert(props_t::has<Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>> ==
                true);

  const auto props = Kokkos::Experimental::node_props(
      Kokkos::Experimental::get_device_handle(TEST_EXECSPACE{}), "label");

  ASSERT_EQ(
      Kokkos::Impl::get_property<Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>(
          props),
      Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>{});
  ASSERT_EQ(Kokkos::Impl::get_property<std::string>(props), "label");
}

TEST(TEST_CATEGORY, node_props_label_device_handle) {
  using props_t = decltype(Kokkos::Experimental::node_props(
      "label",
      Kokkos::Experimental::get_device_handle(std::declval<TEST_EXECSPACE>())));

  static_assert(
      std::same_as<props_t, Kokkos::Impl::NodeCtorProps<
                                std::string,
                                Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>>);

  static_assert(Kokkos::Impl::is_node_props_v<props_t>);

  static_assert(props_t::has<std::string> == true);
  static_assert(props_t::has<Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>> ==
                true);

  const auto props = Kokkos::Experimental::node_props(
      "label", Kokkos::Experimental::get_device_handle(TEST_EXECSPACE{}));

  ASSERT_EQ(Kokkos::Impl::get_property<std::string>(props), "label");
  ASSERT_EQ(
      Kokkos::Impl::get_property<Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>(
          props),
      Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>{});
}

TEST(TEST_CATEGORY, node_props_with_properties_if_unset) {
  {
    const auto props_label_device_handle_old = Kokkos::Experimental::node_props(
        "label", Kokkos::Experimental::get_device_handle(TEST_EXECSPACE{}));
    const auto props_label_device_handle_new =
        Kokkos::Impl::with_properties_if_unset(props_label_device_handle_old,
                                               "another label");

    static_assert(
        std::same_as<
            decltype(props_label_device_handle_old),
            const Kokkos::Impl::NodeCtorProps<
                std::string, Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>>);
    static_assert(
        std::same_as<
            decltype(props_label_device_handle_new),
            const Kokkos::Impl::NodeCtorProps<
                std::string, Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>>);

    ASSERT_EQ(
        Kokkos::Impl::get_property<std::string>(props_label_device_handle_old),
        "label");
    ASSERT_EQ(
        Kokkos::Impl::get_property<std::string>(props_label_device_handle_new),
        "label");
  }

  {
    const auto props_empty = Kokkos::Experimental::node_props();
    const auto props_label_device_handle =
        Kokkos::Impl::with_properties_if_unset(
            props_empty, "label", Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>{});

    static_assert(
        std::same_as<
            decltype(props_label_device_handle),
            const Kokkos::Impl::NodeCtorProps<
                std::string, Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>>);

    ASSERT_EQ(
        Kokkos::Impl::get_property<std::string>(props_label_device_handle),
        "label");
    ASSERT_EQ(
        Kokkos::Impl::get_property<Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>(
            props_label_device_handle),
        Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>{});
  }
}

TEST(TEST_CATEGORY, node_props_get_property) {
  auto props = Kokkos::Experimental::node_props(
      "label", Kokkos::Experimental::get_device_handle(TEST_EXECSPACE{}));

  ASSERT_EQ(Kokkos::Impl::get_property<std::string>(props), "label");
  ASSERT_EQ(
      Kokkos::Impl::get_property<Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>>(
          props),
      Kokkos::Impl::DeviceHandle<TEST_EXECSPACE>{});

  const auto label = Kokkos::Impl::extract_property<std::string>(props);

  ASSERT_EQ(label, "label");
  ASSERT_TRUE(Kokkos::Impl::get_property<std::string>(props).empty());
}

}  // end namespace
