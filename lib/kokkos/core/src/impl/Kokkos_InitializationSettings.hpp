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

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#define KOKKOS_IMPL_DECLARE_DEPRECATED(TYPE, NAME)                         \
 private:                                                                  \
  std::optional<TYPE> m_##NAME;                                            \
                                                                           \
 public:                                                                   \
  KOKKOS_DEPRECATED InitializationSettings& set_##NAME(TYPE NAME) {        \
    m_##NAME = NAME;                                                       \
    return *this;                                                          \
  }                                                                        \
  KOKKOS_DEPRECATED bool has_##NAME() const noexcept {                     \
    return static_cast<bool>(m_##NAME);                                    \
  }                                                                        \
  KOKKOS_DEPRECATED TYPE get_##NAME() const noexcept { return *m_##NAME; } \
  static_assert(true, "no-op to require trailing semicolon")
#else
#define KOKKOS_IMPL_DECLARE_DEPRECATED(TYPE, NAME) \
  static_assert(true, "no-op to require trailing semicolon")
#endif

 public:
  KOKKOS_IMPL_DECLARE(int, num_threads);
  KOKKOS_IMPL_DECLARE(int, device_id);
  KOKKOS_IMPL_DECLARE(std::string, map_device_id_by);
  KOKKOS_IMPL_DECLARE_DEPRECATED(int, num_devices);
  KOKKOS_IMPL_DECLARE_DEPRECATED(int, skip_device);
  KOKKOS_IMPL_DECLARE(bool, disable_warnings);
  KOKKOS_IMPL_DECLARE(bool, print_configuration);
  KOKKOS_IMPL_DECLARE(bool, tune_internals);
  KOKKOS_IMPL_DECLARE(bool, tools_help);
  KOKKOS_IMPL_DECLARE(std::string, tools_libs);
  KOKKOS_IMPL_DECLARE(std::string, tools_args);

#undef KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER_TYPE
#undef KOKKOS_IMPL_INIT_ARGS_DATA_MEMBER
#undef KOKKOS_IMPL_DECLARE
};

}  // namespace Kokkos

#endif
