// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEXTSILICON_ABORT_HPP
#define KOKKOS_NEXTSILICON_ABORT_HPP

#include <nextapi/intrinsics.h>
#include <cstdio>

namespace Kokkos {
namespace Impl {

[[noreturn]] inline void nextsilicon_abort(char const* msg) {
  // FIXME_NEXTSILICON: We need to add printing from handed off code also.
  if (!__next_is_in_handed_off_code()) {
    std::fprintf(stderr, "%s", msg);
    std::fflush(stderr);
  }

  __builtin_trap();
}

}  // namespace Impl
}  // namespace Kokkos

#endif
