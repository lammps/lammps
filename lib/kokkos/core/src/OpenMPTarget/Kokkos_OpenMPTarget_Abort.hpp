// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENMPTARGET_ABORT_HPP
#define KOKKOS_OPENMPTARGET_ABORT_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_OPENMPTARGET

namespace Kokkos {
namespace Impl {

KOKKOS_INLINE_FUNCTION void OpenMPTarget_abort(char const *msg) {
  fprintf(stderr, "%s.\n", msg);
  std::abort();
}

}  // namespace Impl
}  // namespace Kokkos

#endif
#endif
