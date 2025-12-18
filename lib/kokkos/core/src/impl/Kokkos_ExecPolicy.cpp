// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <sstream>
#include <Kokkos_ExecPolicy.hpp>

namespace Kokkos {
namespace Impl {
PerTeamValue::PerTeamValue(size_t arg) : value(arg) {}

PerThreadValue::PerThreadValue(size_t arg) : value(arg) {}
}  // namespace Impl

Impl::PerTeamValue PerTeam(const size_t& arg) {
  return Impl::PerTeamValue(arg);
}

Impl::PerThreadValue PerThread(const size_t& arg) {
  return Impl::PerThreadValue(arg);
}

void team_policy_check_valid_storage_level_argument(int level) {
  if (!(level == 0 || level == 1)) {
    std::stringstream ss;
    ss << "TeamPolicy::set_scratch_size(/*level*/ " << level
       << ", ...) storage level argument must be 0 or 1 to be valid\n";
    abort(ss.str().c_str());
  }
}

}  // namespace Kokkos
