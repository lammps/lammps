// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <impl/Kokkos_InitializationSettings.hpp>

namespace {

TEST(defaultdevicetype, initialization_settings) {
  auto const settings = Kokkos::InitializationSettings()
                            .set_num_threads(255)
                            .set_disable_warnings(false)
                            .set_tools_libs("my_custom_tool.so");
  EXPECT_TRUE(settings.has_num_threads());
  EXPECT_EQ(settings.get_num_threads(), 255);
  EXPECT_FALSE(settings.has_device_id());
  EXPECT_TRUE(settings.has_disable_warnings());
  EXPECT_FALSE(settings.get_disable_warnings());
  EXPECT_FALSE(settings.has_tune_internals());
  EXPECT_FALSE(settings.has_tools_help());
  EXPECT_TRUE(settings.has_tools_libs());
  EXPECT_EQ(settings.get_tools_libs(), "my_custom_tool.so");
  EXPECT_FALSE(settings.has_tools_args());
}

constexpr bool test_initialization_settings_getter() {
#define CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(NAME, TYPE)           \
  static_assert(std::is_same_v<                                                \
                decltype(std::declval<Kokkos::InitializationSettings const&>() \
                             .has_##NAME()),                                   \
                bool>);                                                        \
  static_assert(std::is_same_v<                                                \
                decltype(std::declval<Kokkos::InitializationSettings const&>() \
                             .get_##NAME()),                                   \
                TYPE>);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(num_threads, int);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(device_id, int);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(disable_warnings, bool);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(tune_internals, bool);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(tools_help, bool);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(tools_libs, std::string);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(tools_args, std::string);
#undef CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE
  return true;
}

static_assert(test_initialization_settings_getter());

static_assert(std::is_default_constructible_v<Kokkos::InitializationSettings>);

}  // namespace
