/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>

#include <impl/Kokkos_InitializationSettings.hpp>

namespace {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
void take_initialization_settings(Kokkos::InitializationSettings const&) {}

TEST(defaultdevicetype,
     init_arguments_implicit_conversion_to_initialization_settings) {
  Kokkos::InitArguments arguments;
  take_initialization_settings(arguments);  // check that conversion is implicit
  arguments.device_id      = 1;
  arguments.tune_internals = true;
  Kokkos::InitializationSettings settings{arguments};
  EXPECT_FALSE(settings.has_num_threads());
  EXPECT_TRUE(settings.has_device_id());
  EXPECT_EQ(settings.get_device_id(), 1);
  EXPECT_FALSE(settings.has_num_devices());
  EXPECT_FALSE(settings.has_skip_device());
  EXPECT_FALSE(settings.has_disable_warnings());
  EXPECT_TRUE(settings.has_tune_internals());
  EXPECT_TRUE(settings.get_tune_internals());
  EXPECT_FALSE(settings.has_tools_help());
  EXPECT_FALSE(settings.has_tools_libs());
  EXPECT_FALSE(settings.has_tools_args());
}
#endif

TEST(defaultdevicetype, initialization_settings) {
  auto const settings = Kokkos::InitializationSettings()
                            .set_num_threads(255)
                            .set_disable_warnings(false)
                            .set_tools_libs("my_custom_tool.so");
  EXPECT_TRUE(settings.has_num_threads());
  EXPECT_EQ(settings.get_num_threads(), 255);
  EXPECT_FALSE(settings.has_device_id());
  EXPECT_FALSE(settings.has_num_devices());
  EXPECT_FALSE(settings.has_skip_device());
  EXPECT_TRUE(settings.has_disable_warnings());
  EXPECT_FALSE(settings.get_disable_warnings());
  EXPECT_FALSE(settings.has_tune_internals());
  EXPECT_FALSE(settings.has_tools_help());
  EXPECT_TRUE(settings.has_tools_libs());
  EXPECT_EQ(settings.get_tools_libs(), "my_custom_tool.so");
  EXPECT_FALSE(settings.has_tools_args());
}

#define STATIC_ASSERT(...) static_assert(__VA_ARGS__, "")  // FIXME C++17

constexpr bool test_initialization_settings_getter() {
#define CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(NAME, TYPE)           \
  STATIC_ASSERT(std::is_same<                                                  \
                decltype(std::declval<Kokkos::InitializationSettings const&>() \
                             .has_##NAME()),                                   \
                bool>::value);                                                 \
  STATIC_ASSERT(std::is_same<                                                  \
                decltype(std::declval<Kokkos::InitializationSettings const&>() \
                             .get_##NAME()),                                   \
                TYPE>::value);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(num_threads, int);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(device_id, int);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(num_devices, int);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(skip_device, int);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(disable_warnings, bool);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(tune_internals, bool);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(tools_help, bool);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(tools_libs, std::string);
  CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE(tools_args, std::string);
#undef CHECK_INITIALIZATION_SETTINGS_GETTER_RETURN_TYPE
  return true;
}

STATIC_ASSERT(test_initialization_settings_getter());

STATIC_ASSERT(
    std::is_default_constructible<Kokkos::InitializationSettings>::value);

}  // namespace
