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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Core.hpp>
#include <sstream>

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
    Impl::throw_runtime_exception(ss.str());
  }
}

}  // namespace Kokkos
