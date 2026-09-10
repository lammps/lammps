// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_INITIALIZATION_SETTINGS_HPP
#define KOKKOS_INITIALIZATION_SETTINGS_HPP

#include <Kokkos_Macros.hpp>

#include <optional>
#include <string>

namespace Kokkos {

class InitializationSettings {
#define KOKKOS_IMPL_DECLARE(TYPE, NAME)                                      \
 private:                                                                    \
  std::optional<TYPE> m_##NAME;                                              \
                                                                             \
 public:                                                                     \
  InitializationSettings& set_##NAME(TYPE NAME) {                            \
    m_##NAME = NAME;                                                         \
    return *this;                                                            \
  }                                                                          \
  bool has_##NAME() const noexcept { return static_cast<bool>(m_##NAME); }   \
  TYPE get_##NAME() const noexcept /* NOLINT(bugprone-exception-escape) */ { \
    return *m_##NAME; /* NOLINT(bugprone-unchecked-optional-access) */       \
  }                                                                          \
  static_assert(true, "no-op to require trailing semicolon")

 public:
  KOKKOS_IMPL_DECLARE(int, num_threads);
  KOKKOS_IMPL_DECLARE(int, device_id);
  KOKKOS_IMPL_DECLARE(std::string, map_device_id_by);
  KOKKOS_IMPL_DECLARE(bool, disable_warnings);
  KOKKOS_IMPL_DECLARE(bool, print_configuration);
  KOKKOS_IMPL_DECLARE(bool, tune_internals);
  KOKKOS_IMPL_DECLARE(bool, tools_help);
  KOKKOS_IMPL_DECLARE(std::string, tools_libs);
  KOKKOS_IMPL_DECLARE(std::string, tools_args);

#undef KOKKOS_IMPL_DECLARE
};

}  // namespace Kokkos

#endif
