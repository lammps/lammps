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

#ifndef KOKKOS_INITIALIZATION_SETTINGS_HPP
#define KOKKOS_INITIALIZATION_SETTINGS_HPP

#include <Kokkos_Macros.hpp>

#include <optional>
#include <string>

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
struct InitArguments {
  int num_threads;
  int num_numa;
  int device_id;
  int ndevices;
  int skip_device;
  bool disable_warnings;
  bool tune_internals;
  bool tool_help        = false;
  std::string tool_lib  = {};
  std::string tool_args = {};

  KOKKOS_DEPRECATED_WITH_COMMENT("Use InitializationSettings instead!")
  InitArguments(int nt = -1, int nn = -1, int dv = -1, bool dw = false,
                bool ti = false)
      : num_threads{nt},
        num_numa{nn},
        device_id{dv},
        ndevices{-1},
        skip_device{9999},
        disable_warnings{dw},
        tune_internals{ti} {}
};
#endif

class InitializationSettings {
#define KOKKOS_IMPL_DECLARE(TYPE, NAME)                                    \
 private:                                                                  \
  std::optional<TYPE> m_##NAME;                                            \
                                                                           \
 public:                                                                   \
  InitializationSettings& set_##NAME(TYPE NAME) {                          \
    m_##NAME = NAME;                                                       \
    return *this;                                                          \
  }                                                                        \
  bool has_##NAME() const noexcept { return static_cast<bool>(m_##NAME); } \
  TYPE get_##NAME() const noexcept { return *m_##NAME; }                   \
  static_assert(true, "no-op to require trailing semicolon")

 public:
  KOKKOS_IMPL_DECLARE(int, num_threads);
  KOKKOS_IMPL_DECLARE(int, device_id);
  KOKKOS_IMPL_DECLARE(std::string, map_device_id_by);
  KOKKOS_IMPL_DECLARE(int, num_devices);  // deprecated
  KOKKOS_IMPL_DECLARE(int, skip_device);  // deprecated
  KOKKOS_IMPL_DECLARE(bool, disable_warnings);
  KOKKOS_IMPL_DECLARE(bool, print_configuration);
  KOKKOS_IMPL_DECLARE(bool, tune_internals);
  KOKKOS_IMPL_DECLARE(bool, tools_help);
  KOKKOS_IMPL_DECLARE(std::string, tools_libs);
  KOKKOS_IMPL_DECLARE(std::string, tools_args);

#undef KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER_TYPE
#undef KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER
#undef KOKKOS_IMPL_DECLARE

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
 public:
  InitializationSettings() = default;

  InitializationSettings(InitArguments const& old) {
    if (old.num_threads != -1) {
      set_num_threads(old.num_threads);
    }
    if (old.device_id != -1) {
      set_device_id(old.device_id);
    }
    if (old.ndevices != -1) {
      set_num_devices(old.ndevices);
    }
    if (old.skip_device != 9999) {
      set_skip_device(old.skip_device);
    }
    if (old.disable_warnings) {
      set_disable_warnings(true);
    }
    if (old.tune_internals) {
      set_tune_internals(true);
    }
    if (old.tool_help) {
      set_tools_help(true);
    }
    if (!old.tool_lib.empty()) {
      set_tools_libs(old.tool_lib);
    }
    if (!old.tool_args.empty()) {
      set_tools_args(old.tool_args);
    }
  }
#endif
};

}  // namespace Kokkos

#endif
