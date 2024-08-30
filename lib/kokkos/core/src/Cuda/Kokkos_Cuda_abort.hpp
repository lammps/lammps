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
extern __device__ void __assertfail(const void *message, const void *file,
                                    unsigned int line, const void *function,
                                    size_t charsize);
}

namespace Kokkos {
namespace Impl {

// required to workaround failures in random number generator unit tests with
// pre-volta architectures
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
__device__ inline void cuda_abort(const char *const message) {
#else
[[noreturn]] __device__ inline void cuda_abort(const char *const message) {
#endif
  const char empty[] = "";

  __assertfail((const void *)message, (const void *)empty, (unsigned int)0,
               (const void *)empty, sizeof(char));

  // This loop is never executed. It's intended to suppress warnings that the
  // function returns, even though it does not. This is necessary because
  // __assertfail is not marked as [[noreturn]], even though it does not return.
  //  Disable with KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK to workaround failures
  //  in random number generator unit tests with pre-volta architectures
#if !defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
  while (true)
    ;
#endif
}

}  // namespace Impl
}  // namespace Kokkos
#else
void KOKKOS_CORE_SRC_CUDA_ABORT_PREVENT_LINK_ERROR() {}
#endif /* #if defined( KOKKOS_ENABLE_CUDA ) */
#endif /* #ifndef KOKKOS_CUDA_ABORT_HPP */
