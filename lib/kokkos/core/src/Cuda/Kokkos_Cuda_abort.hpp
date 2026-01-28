// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CUDA_ABORT_HPP
#define KOKKOS_CUDA_ABORT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_CUDA)

#include <cuda.h>

extern "C" {
/*  Cuda runtime function, declared in <crt/device_runtime.h>
 *  Requires capability 2.x or better.
 */
// NOLINTNEXTLINE(bugprone-reserved-identifier)
[[noreturn]] __device__ void __assertfail(const void *message, const void *file,
                                          unsigned int line,
                                          const void *function,
                                          size_t charsize);
}

namespace Kokkos {
namespace Impl {

[[noreturn]] __device__ static void cuda_abort(const char *const message) {
  const char empty[] = "";

  __assertfail((const void *)message, (const void *)empty, (unsigned int)0,
               (const void *)empty, sizeof(char));
}

}  // namespace Impl
}  // namespace Kokkos
#else
void KOKKOS_CORE_SRC_CUDA_ABORT_PREVENT_LINK_ERROR() {}
#endif /* #if defined( KOKKOS_ENABLE_CUDA ) */
#endif /* #ifndef KOKKOS_CUDA_ABORT_HPP */
