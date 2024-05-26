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

#ifndef KOKKOS_ABORT_HPP
#define KOKKOS_ABORT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Printf.hpp>
#ifdef KOKKOS_ENABLE_CUDA
#include <Cuda/Kokkos_Cuda_abort.hpp>
#endif
#ifdef KOKKOS_ENABLE_HIP
#include <HIP/Kokkos_HIP_Abort.hpp>
#endif
#ifdef KOKKOS_ENABLE_SYCL
#include <SYCL/Kokkos_SYCL_Abort.hpp>
#endif

namespace Kokkos {
namespace Impl {

[[noreturn]] void host_abort(const char *const);

#if defined(KOKKOS_ENABLE_CUDA) && defined(__CUDA_ARCH__)

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
// required to workaround failures in random number generator unit tests with
// pre-volta architectures
#define KOKKOS_IMPL_ABORT_NORETURN
#else
// cuda_abort aborts when building for other platforms than macOS
#define KOKKOS_IMPL_ABORT_NORETURN [[noreturn]]
#endif

#elif defined(KOKKOS_COMPILER_NVHPC)

#define KOKKOS_IMPL_ABORT_NORETURN

#elif defined(KOKKOS_ENABLE_HIP) && defined(__HIP_DEVICE_COMPILE__)
// HIP aborts
#define KOKKOS_IMPL_ABORT_NORETURN [[noreturn]]
#elif defined(KOKKOS_ENABLE_SYCL) && defined(__SYCL_DEVICE_ONLY__)
// FIXME_SYCL SYCL doesn't abort
#define KOKKOS_IMPL_ABORT_NORETURN
#elif !defined(KOKKOS_ENABLE_OPENMPTARGET) && !defined(KOKKOS_ENABLE_OPENACC)
// Host aborts
#define KOKKOS_IMPL_ABORT_NORETURN [[noreturn]]
#else
// Everything else does not abort
#define KOKKOS_IMPL_ABORT_NORETURN
#endif

// FIXME_SYCL
// Accomodate host pass for device functions that are not [[noreturn]]
#if defined(KOKKOS_ENABLE_SYCL) || \
    (defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK))
#define KOKKOS_IMPL_ABORT_NORETURN_DEVICE
#else
#define KOKKOS_IMPL_ABORT_NORETURN_DEVICE KOKKOS_IMPL_ABORT_NORETURN
#endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) ||          \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET) || \
    defined(KOKKOS_ENABLE_OPENACC)
KOKKOS_IMPL_ABORT_NORETURN_DEVICE inline KOKKOS_IMPL_DEVICE_FUNCTION void
device_abort(const char *const msg) {
#if defined(KOKKOS_ENABLE_CUDA)
  ::Kokkos::Impl::cuda_abort(msg);
#elif defined(KOKKOS_ENABLE_HIP)
  ::Kokkos::Impl::hip_abort(msg);
#elif defined(KOKKOS_ENABLE_SYCL)
  ::Kokkos::Impl::sycl_abort(msg);
#elif defined(KOKKOS_ENABLE_OPENMPTARGET) || defined(KOKKOS_ENABLE_OPENACC)
  printf("%s", msg);  // FIXME_OPENMPTARGET FIXME_OPENACC
#else
#error faulty logic
#endif
}
#endif
}  // namespace Impl

KOKKOS_IMPL_ABORT_NORETURN KOKKOS_INLINE_FUNCTION void abort(
    const char *const message) {
  KOKKOS_IF_ON_HOST(::Kokkos::Impl::host_abort(message);)
  KOKKOS_IF_ON_DEVICE(::Kokkos::Impl::device_abort(message);)
}

#undef KOKKOS_IMPL_ABORT_NORETURN

}  // namespace Kokkos

#endif /* #ifndef KOKKOS_ABORT_HPP */
