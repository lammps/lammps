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

#ifndef KOKKOS_CLOCKTIC_HPP
#define KOKKOS_CLOCKTIC_HPP

#include <Kokkos_Macros.hpp>
#include <stdint.h>
#include <chrono>
#ifdef KOKKOS_ENABLE_OPENMPTARGET
#include <omp.h>
#endif

// To use OpenCL(TM) built-in intrinsics inside kernels, we have to
// forward-declare their prototype, also see
// https://github.com/intel/pti-gpu/blob/master/chapters/binary_instrumentation/OpenCLBuiltIn.md
#if defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_ARCH_INTEL_GPU) && \
    defined(__SYCL_DEVICE_ONLY__)
extern SYCL_EXTERNAL unsigned long __attribute__((overloadable))
intel_get_cycle_counter();
#endif

namespace Kokkos {
namespace Impl {

/**\brief  Quick query of clock register tics
 *
 *  Primary use case is to, with low overhead,
 *  obtain a integral value that consistently varies
 *  across concurrent threads of execution within
 *  a parallel algorithm.
 *  This value is often used to "randomly" seed an
 *  attempt to acquire an indexed resource (e.g., bit)
 *  from an array of resources (e.g., bitset) such that
 *  concurrent threads will have high likelihood of
 *  having different index-seed values.
 */

KOKKOS_IMPL_DEVICE_FUNCTION inline uint64_t clock_tic_device() noexcept {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)

  // Return value of 64-bit hi-res clock register.
  return clock64();

#elif defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_ARCH_INTEL_GPU) && \
    defined(__SYCL_DEVICE_ONLY__)

  return intel_get_cycle_counter();

#elif defined(KOKKOS_ENABLE_OPENMPTARGET)

  return omp_get_wtime() * 1.e9;

#else

  return 0;

#endif
}

KOKKOS_IMPL_HOST_FUNCTION inline uint64_t clock_tic_host() noexcept {
#if defined(__i386__) || defined(__x86_64)

  // Return value of 64-bit hi-res clock register.

  unsigned a = 0, d = 0;

  __asm__ volatile("rdtsc" : "=a"(a), "=d"(d));

  return ((uint64_t)a) | (((uint64_t)d) << 32);

#elif defined(__powerpc64__) || defined(__ppc64__)

  unsigned long cycles = 0;

  asm volatile("mftb %0" : "=r"(cycles));

  return (uint64_t)cycles;

#elif defined(__ppc__)
  // see : pages.cs.wisc.edu/~legault/miniproj-736.pdf or
  // cmssdt.cern.ch/lxr/source/FWCore/Utilities/interface/HRRealTime.h

  uint64_t result = 0;
  uint32_t upper, lower, tmp;

  __asm__ volatile(
      "0: \n"
      "\tmftbu %0     \n"
      "\tmftb  %1     \n"
      "\tmftbu %2     \n"
      "\tcmpw  %2, %0 \n"
      "\tbne   0b     \n"
      : "=r"(upper), "=r"(lower), "=r"(tmp));
  result = upper;
  result = result << 32;
  result = result | lower;

  return (result);

#else

  return std::chrono::high_resolution_clock::now().time_since_epoch().count();

#endif
}

KOKKOS_FORCEINLINE_FUNCTION
uint64_t clock_tic() noexcept {
  KOKKOS_IF_ON_DEVICE((return clock_tic_device();))
  KOKKOS_IF_ON_HOST((return clock_tic_host();))
}

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_CLOCKTIC_HPP
