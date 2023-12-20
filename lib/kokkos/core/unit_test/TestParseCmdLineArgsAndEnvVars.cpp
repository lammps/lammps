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

#include <impl/Kokkos_ParseCommandLineArgumentsAndEnvironmentVariables.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>
#include <impl/Kokkos_DeviceManagement.hpp>
#include <impl/Kokkos_Command_Line_Parsing.hpp>

#include <cstdlib>
#include <memory>
#include <mutex>
#include <regex>
#include <string>
#include <unordered_map>

namespace {

class EnvVarsHelper {
  // do not let GTest run unit tests that set the environment concurrently
  static std::mutex mutex_;
  std::vector<std::string> vars_;
  // FIXME_CXX17 prefer optional
  // store name of env var that was already set (if any)
  // in which case unit test is skipped
  std::unique_ptr<std::string> skip_;

  void setup(std::unordered_map<std::string, std::string> const& vars) {
    for (auto const& x : vars) {
      auto const& name  = x.first;
      auto const& value = x.second;
      // skip unit test if env var is already set
      if (getenv(name.c_str())) {
        skip_ = std::make_unique<std::string>(name);
        break;
      }
#ifdef _WIN32
      int const error_code = _putenv((name + "=" + value).c_str());
#else
      int const error_code =
          setenv(name.c_str(), value.c_str(), /*overwrite=*/0);
#endif
      if (error_code != 0) {
        std::cerr << "failed to set environment variable '" << name << "="
                  << value << "'\n";
        std::abort();
      }
      vars_.push_back(name);
    }
  }
  void teardown() {
    for (auto const& name : vars_) {
#ifdef _WIN32
      int const error_code = _putenv((name + "=").c_str());
#else
      int const error_code = unsetenv(name.c_str());
#endif
      if (error_code != 0) {
        std::cerr << "failed to unset environment variable '" << name << "'\n";
        std::abort();
      }
    }
  }

 public:
  auto& skip() { return skip_; }
  EnvVarsHelper(std::unordered_map<std::string, std::string> const& vars) {
    mutex_.lock();
    setup(vars);
  }
  EnvVarsHelper& operator=(
      std::unordered_map<std::string, std::string> const& vars) {
    teardown();
    setup(vars);
    return *this;
  }
  ~EnvVarsHelper() {
    teardown();
    mutex_.unlock();
  }
  EnvVarsHelper(EnvVarsHelper&) = delete;
  EnvVarsHelper& operator=(EnvVarsHelper&) = delete;
  friend std::ostream& operator<<(std::ostream& os, EnvVarsHelper const& ev) {
    for (auto const& name : ev.vars_) {
      os << name << '=' << std::getenv(name.c_str()) << '\n';
    }
    return os;
  }
};
std::mutex EnvVarsHelper::mutex_;
#define SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev)       \
  if (ev.skip()) {                                         \
    GTEST_SKIP() << "environment variable '" << *ev.skip() \
                 << "' is already set";                    \
  }                                                        \
  static_assert(true, "no-op to require trailing semicolon")

class CmdLineArgsHelper {
  int argc_;
  std::vector<char*> argv_;
  std::vector<std::unique_ptr<char[]>> args_;

 public:
  CmdLineArgsHelper(std::vector<std::string> const& args) : argc_(args.size()) {
    for (auto const& x : args) {
      args_.emplace_back(new char[x.size() + 1]);
      char* ptr = args_.back().get();
      strcpy(ptr, x.c_str());
      argv_.push_back(ptr);
    }
    argv_.push_back(nullptr);
  }
  int& argc() { return argc_; }
  char** argv() { return argv_.data(); }
};
#define EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, ...) \
  do {                                                    \
    std::vector<std::string> expected_argv = __VA_ARGS__; \
                                                          \
    int expected_argc = expected_argv.size();             \
    EXPECT_EQ(cla.argc(), expected_argc);                 \
    for (int i = 0; i < expected_argc; ++i) {             \
      EXPECT_EQ(cla.argv()[i], expected_argv[i])          \
          << "arguments differ at index " << i;           \
    }                                                     \
    EXPECT_EQ(cla.argv()[cla.argc()], nullptr);           \
  } while (false)

TEST(defaultdevicetype, cmd_line_args_num_threads) {
  CmdLineArgsHelper cla = {{
      "--foo=bar",
      "--kokkos-num-threads=1",
      "--kokkos-num-threads=2",
  }};
  Kokkos::InitializationSettings settings;
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  EXPECT_TRUE(settings.has_num_threads());
  EXPECT_EQ(settings.get_num_threads(), 2);
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {"--foo=bar"});
}

TEST(defaultdevicetype, cmd_line_args_device_id) {
  CmdLineArgsHelper cla = {{
      "--kokkos-device-id=3",
      "--dummy",
      "--kokkos-device-id=4",
  }};
  Kokkos::InitializationSettings settings;
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  EXPECT_TRUE(settings.has_device_id());
  EXPECT_EQ(settings.get_device_id(), 4);
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {"--dummy"});
}

TEST(defaultdevicetype, cmd_line_args_num_devices) {
  CmdLineArgsHelper cla = {{
      "--kokkos-num-devices=5,6",
      "--kokkos-num-devices=7",
      "-v",
  }};
  Kokkos::InitializationSettings settings;
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  EXPECT_TRUE(settings.has_num_devices());
  EXPECT_EQ(settings.get_num_devices(), 7);
  // this is the current behavior, not suggesting this cannot be revisited
  EXPECT_TRUE(settings.has_skip_device()) << "behavior changed see comment";
  EXPECT_EQ(settings.get_skip_device(), 6) << "behavior changed see comment";
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {"-v"});
}

TEST(defaultdevicetype, cmd_line_args_disable_warning) {
  CmdLineArgsHelper cla = {{
      "--kokkos-disable-warnings=1",
      "--kokkos-disable-warnings=false",
  }};
  Kokkos::InitializationSettings settings;
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  EXPECT_TRUE(settings.has_disable_warnings());
  EXPECT_FALSE(settings.get_disable_warnings());
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {});
}

TEST(defaultdevicetype, cmd_line_args_tune_internals) {
  CmdLineArgsHelper cla = {{
      "--kokkos-tune-internals",
      "--kokkos-num-threads=3",
  }};
  Kokkos::InitializationSettings settings;
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  EXPECT_TRUE(settings.has_tune_internals());
  EXPECT_TRUE(settings.get_tune_internals());
  EXPECT_TRUE(settings.has_num_threads());
  EXPECT_EQ(settings.get_num_threads(), 3);
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {});
}

TEST(defaultdevicetype, cmd_line_args_help) {
  CmdLineArgsHelper cla = {{
      "--help",
  }};
  Kokkos::InitializationSettings settings;
  ::testing::internal::CaptureStdout();
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  auto captured = ::testing::internal::GetCapturedStdout();
  // check that error message was only printed once
  EXPECT_EQ(captured.find("--kokkos-help"), captured.rfind("--kokkos-help"));
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {"--help"});
  auto const help_message_length = captured.length();

  cla = {{
      {"--kokkos-help"},
  }};
  ::testing::internal::CaptureStdout();
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  captured = ::testing::internal::GetCapturedStdout();
  EXPECT_EQ(captured.length(), help_message_length);
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {});

  cla = {{
      {"--kokkos-help"},
      {"--help"},
      {"--kokkos-help"},
  }};
  ::testing::internal::CaptureStdout();
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  captured = ::testing::internal::GetCapturedStdout();
  EXPECT_EQ(captured.length(), help_message_length);
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {"--help"});
}

TEST(defaultdevicetype, cmd_line_args_tools_arguments) {
  CmdLineArgsHelper cla = {{
      "--kokkos-tool-libs=ich_tue_nur.so",
  }};
  Kokkos::InitializationSettings settings;
  ::testing::internal::CaptureStderr();
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  auto captured = ::testing::internal::GetCapturedStderr();
  EXPECT_TRUE(captured.find("not recognized") != std::string::npos &&
              captured.find("--kokkos-tool-libs=ich_tue_nur.so") !=
                  std::string::npos &&
              !settings.has_tools_libs())
      << captured;
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(
      cla, {"--kokkos-tool-libs=ich_tue_nur.so"});

  cla      = {{
      "--kokkos-tools-libs=ich_tue_nur.so",
  }};
  settings = {};
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  EXPECT_TRUE(settings.has_tools_libs());
  EXPECT_EQ(settings.get_tools_libs(), "ich_tue_nur.so");
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {});
}

TEST(defaultdevicetype, cmd_line_args_unrecognized_flag) {
  CmdLineArgsHelper cla = {{
      "--kokkos_num_threads=4",  // underscores instead of dashes
  }};
  Kokkos::InitializationSettings settings;
  ::testing::internal::CaptureStderr();
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  auto captured = ::testing::internal::GetCapturedStderr();
  EXPECT_TRUE(captured.find("not recognized") != std::string::npos &&
              captured.find("--kokkos_num_threads=4") != std::string::npos &&
              !settings.has_num_threads())
      << captured;
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {"--kokkos_num_threads=4"});

  cla = {{
      "-kokkos-num-threads=4",  // missing one leading dash
  }};
  ::testing::internal::CaptureStderr();
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  captured = ::testing::internal::GetCapturedStderr();
  EXPECT_TRUE(captured.find("not recognized") != std::string::npos &&
              captured.find("-kokkos-num-threads=4") != std::string::npos &&
              !settings.has_num_threads())
      << captured;
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {"-kokkos-num-threads=4"});

  cla = {{
      "--kokko-num-threads=4",  // no warning when prefix misspelled
  }};
  ::testing::internal::CaptureStderr();
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  captured = ::testing::internal::GetCapturedStderr();
  EXPECT_TRUE(captured.empty() && !settings.has_num_threads()) << captured;
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla, {"--kokko-num-threads=4"});

  Kokkos::Impl::do_not_warn_not_recognized_command_line_argument(
      std::regex{"^--kokkos-extension.*"});
  cla = {{
      "--kokkos-extension-option=value",  // user explicitly asked not to warn
                                          // about that prefix
  }};
  ::testing::internal::CaptureStderr();
  Kokkos::Impl::parse_command_line_arguments(cla.argc(), cla.argv(), settings);
  captured = ::testing::internal::GetCapturedStderr();
  EXPECT_TRUE(captured.empty()) << captured;
  EXPECT_REMAINING_COMMAND_LINE_ARGUMENTS(cla,
                                          {"--kokkos-extension-option=value"});
}

TEST(defaultdevicetype, env_vars_num_threads) {
  EnvVarsHelper ev = {{
      {"KOKKOS_NUM_THREADS", "24"},
      {"KOKKOS_DISABLE_WARNINGS", "1"},
  }};
  SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev);
  Kokkos::InitializationSettings settings;
  Kokkos::Impl::parse_environment_variables(settings);
  EXPECT_TRUE(settings.has_num_threads());
  EXPECT_EQ(settings.get_num_threads(), 24);
  EXPECT_TRUE(settings.has_disable_warnings());
  EXPECT_TRUE(settings.get_disable_warnings());

  ev = {{
      {"KOKKOS_NUM_THREADS", "1ABC"},
  }};
  SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev);
  settings = {};
  Kokkos::Impl::parse_environment_variables(settings);
  EXPECT_TRUE(settings.has_num_threads());
  EXPECT_EQ(settings.get_num_threads(), 1);
}

TEST(defaultdevicetype, env_vars_device_id) {
  EnvVarsHelper ev = {{
      {"KOKKOS_DEVICE_ID", "33"},
  }};
  SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev);
  Kokkos::InitializationSettings settings;
  Kokkos::Impl::parse_environment_variables(settings);
  EXPECT_TRUE(settings.has_device_id());
  EXPECT_EQ(settings.get_device_id(), 33);
}

TEST(defaultdevicetype, env_vars_num_devices) {
  EnvVarsHelper ev = {{
      {"KOKKOS_NUM_DEVICES", "4"},
      {"KOKKOS_SKIP_DEVICE", "1"},
  }};
  SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev);
  Kokkos::InitializationSettings settings;
  Kokkos::Impl::parse_environment_variables(settings);
  EXPECT_TRUE(settings.has_num_devices());
  EXPECT_EQ(settings.get_num_devices(), 4);
  EXPECT_TRUE(settings.has_skip_device());
  EXPECT_EQ(settings.get_skip_device(), 1);
}

TEST(defaultdevicetype, env_vars_disable_warnings) {
  for (auto const& value_true : {"1", "true", "TRUE", "yEs"}) {
    EnvVarsHelper ev = {{
        {"KOKKOS_DISABLE_WARNINGS", value_true},
    }};
    SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev);
    Kokkos::InitializationSettings settings;
    Kokkos::Impl::parse_environment_variables(settings);
    EXPECT_TRUE(settings.has_disable_warnings())
        << "KOKKOS_DISABLE_WARNINGS=" << value_true;
    EXPECT_TRUE(settings.get_disable_warnings())
        << "KOKKOS_DISABLE_WARNINGS=" << value_true;
  }
  for (auto const& value_false : {"0", "fAlse", "No"}) {
    EnvVarsHelper ev = {{
        {"KOKKOS_DISABLE_WARNINGS", value_false},
    }};
    SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev);
    Kokkos::InitializationSettings settings;
    Kokkos::Impl::parse_environment_variables(settings);
    EXPECT_TRUE(settings.has_disable_warnings())
        << "KOKKOS_DISABLE_WARNINGS=" << value_false;
    EXPECT_FALSE(settings.get_disable_warnings())
        << "KOKKOS_DISABLE_WARNINGS=" << value_false;
  }
}

TEST(defaultdevicetype, env_vars_tune_internals) {
  for (auto const& value_true : {"1", "yES", "true", "TRUE", "tRuE"}) {
    EnvVarsHelper ev = {{
        {"KOKKOS_TUNE_INTERNALS", value_true},
    }};
    SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev);
    Kokkos::InitializationSettings settings;
    Kokkos::Impl::parse_environment_variables(settings);
    EXPECT_TRUE(settings.has_tune_internals())
        << "KOKKOS_TUNE_INTERNALS=" << value_true;
    EXPECT_TRUE(settings.get_tune_internals())
        << "KOKKOS_TUNE_INTERNALS=" << value_true;
  }
  for (auto const& value_false : {"0", "false", "no"}) {
    EnvVarsHelper ev = {{
        {"KOKKOS_TUNE_INTERNALS", value_false},
    }};
    SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev);
    Kokkos::InitializationSettings settings;
    Kokkos::Impl::parse_environment_variables(settings);
    EXPECT_TRUE(settings.has_tune_internals())
        << "KOKKOS_TUNE_INTERNALS=" << value_false;
    EXPECT_FALSE(settings.get_tune_internals())
        << "KOKKOS_TUNE_INTERNALS=" << value_false;
  }
}

TEST(defaultdevicetype, visible_devices) {
#define KOKKOS_TEST_VISIBLE_DEVICES(ENV, CNT, DEV)                    \
  do {                                                                \
    EnvVarsHelper ev{ENV};                                            \
    SKIP_IF_ENVIRONMENT_VARIABLE_ALREADY_SET(ev);                     \
    Kokkos::InitializationSettings settings;                          \
    Kokkos::Impl::parse_environment_variables(settings);              \
    auto computed = Kokkos::Impl::get_visible_devices(settings, CNT); \
    std::vector<int> expected = DEV;                                  \
    EXPECT_EQ(expected.size(), computed.size())                       \
        << ev << "device count: " << CNT;                             \
    auto n = std::min<int>(expected.size(), computed.size());         \
    for (int i = 0; i < n; ++i) {                                     \
      EXPECT_EQ(expected[i], computed[i])                             \
          << "devices differ at index " << i << '\n'                  \
          << ev << "device count: " << CNT;                           \
    }                                                                 \
  } while (false)

#define DEV(...) \
  std::vector<int> { __VA_ARGS__ }
#define ENV(...) std::unordered_map<std::string, std::string>{__VA_ARGS__}

  // first test with all environment variables that are involved in determining
  // the visible devices so user set var do not mess up the logic below.
  KOKKOS_TEST_VISIBLE_DEVICES(
      ENV({"KOKKOS_VISIBLE_DEVICES", "2,1"}, {"KOKKOS_NUM_DEVICES", "8"},
          {"KOKKOS_SKIP_DEVICE", "1"}),
      6, DEV(2, 1));
  KOKKOS_TEST_VISIBLE_DEVICES(
      ENV({"KOKKOS_VISIBLE_DEVICES", "2,1"}, {"KOKKOS_NUM_DEVICES", "8"}, ), 6,
      DEV(2, 1));
  KOKKOS_TEST_VISIBLE_DEVICES(ENV({"KOKKOS_NUM_DEVICES", "3"}), 6,
                              DEV(0, 1, 2));
  KOKKOS_TEST_VISIBLE_DEVICES(
      ENV({"KOKKOS_NUM_DEVICES", "4"}, {"KOKKOS_SKIP_DEVICE", "1"}, ), 6,
      DEV(0, 2, 3));
  KOKKOS_TEST_VISIBLE_DEVICES(ENV({"KOKKOS_VISIBLE_DEVICES", "1,3,4"}), 6,
                              DEV(1, 3, 4));
  KOKKOS_TEST_VISIBLE_DEVICES(
      ENV({"KOKKOS_VISIBLE_DEVICES", "2,1"}, {"KOKKOS_SKIP_DEVICE", "1"}, ), 6,
      DEV(2, 1));
  KOKKOS_TEST_VISIBLE_DEVICES(ENV(), 4, DEV(0, 1, 2, 3));

#undef ENV
#undef DEV
#undef KOKKOS_TEST_VISIBLE_DEVICES
}

}  // namespace
