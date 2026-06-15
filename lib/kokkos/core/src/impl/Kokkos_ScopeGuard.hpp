// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SCOPE_GUARD_HPP
#define KOKKOS_SCOPE_GUARD_HPP

#include <Kokkos_Abort.hpp>
#include <impl/Kokkos_InitializeFinalize.hpp>

#include <string>

namespace Kokkos::Impl {

inline std::string scopeguard_correct_usage() {
  return std::string(
      "Do instead:\n"
      "  std::unique_ptr<Kokkos::ScopeGuard> guard =\n"
      "    !Kokkos::is_initialized() && !Kokkos::is_finalized()?\n"
      "    new ScopeGuard(argc,argv) : nullptr;\n");
}

inline std::string scopeguard_create_while_initialized_warning() {
  return std::string(
             "Kokkos Error: Creating a ScopeGuard while Kokkos is initialized "
             "is illegal.\n")
      .append(scopeguard_correct_usage());
}

inline std::string scopeguard_create_after_finalize_warning() {
  return std::string(
             "Kokkos Error: Creating a ScopeGuard after Kokkos was finalized "
             "is illegal.\n")
      .append(scopeguard_correct_usage());
}

inline std::string scopeguard_destruct_after_finalize_warning() {
  return std::string(
             "Kokkos Error: Destroying a ScopeGuard after Kokkos was finalized "
             "is illegal.\n")
      .append(scopeguard_correct_usage());
}

}  // namespace Kokkos::Impl

namespace Kokkos {

class [[nodiscard]] ScopeGuard {
 public:
  template <class... Args>
  ScopeGuard(Args&&... args) {
    if (is_initialized()) {
      Kokkos::abort(
          Impl::scopeguard_create_while_initialized_warning().c_str());
    }
    if (is_finalized()) {
      Kokkos::abort(Impl::scopeguard_create_after_finalize_warning().c_str());
    }
    initialize(static_cast<Args&&>(args)...);
  }

  ~ScopeGuard() {
    if (is_finalized()) {
      Kokkos::abort(Impl::scopeguard_destruct_after_finalize_warning().c_str());
    }
    finalize();
  }

  ScopeGuard& operator=(const ScopeGuard&) = delete;
  ScopeGuard& operator=(ScopeGuard&&)      = delete;
  ScopeGuard(const ScopeGuard&)            = delete;
  ScopeGuard(ScopeGuard&&)                 = delete;
};

}  // namespace Kokkos

#endif
