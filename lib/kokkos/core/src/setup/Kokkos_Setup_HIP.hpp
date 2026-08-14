// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SETUP_HIP_HPP_
#define KOKKOS_SETUP_HIP_HPP_

#if defined(KOKKOS_ENABLE_HIP)

#if defined(KOKKOS_ARCH_AMD_GFX942) && defined(KOKKOS_ARCH_AMD_GFX942_APU)
static_assert(false,
              "Kokkos detected both `KOKKOS_ARCH_AMD_GFX942` and "
              "`KOKKOS_ARCH_AMD_GFX942_APU` which is not allowed.");
#endif

#define KOKKOS_IMPL_HIP_CLANG_WORKAROUND

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

#define KOKKOS_LAMBDA [=] __host__ __device__
#define KOKKOS_CLASS_LAMBDA [ =, *this ] __host__ __device__

#define KOKKOS_DEDUCTION_GUIDE __host__ __device__

#define KOKKOS_IMPL_FORCEINLINE_FUNCTION __device__ __host__ __forceinline__
#define KOKKOS_IMPL_INLINE_FUNCTION __device__ __host__ inline
#define KOKKOS_IMPL_FUNCTION __device__ __host__
#define KOKKOS_IMPL_HOST_FUNCTION __host__
#define KOKKOS_IMPL_DEVICE_FUNCTION __device__

// clang-format off
#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION __device__ __host__
#else
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION @"KOKKOS_RELOCATABLE_FUNCTION requires Kokkos_ENABLE_HIP_RELOCATABLE_DEVICE_CODE=ON"
#endif
// clang-format on

#ifdef KOKKOS_ARCH_AMD_GFX942_APU
#define KOKKOS_IMPL_HIP_UNIFIED_MEMORY
#endif

#endif  // #if defined( KOKKOS_ENABLE_HIP )

#endif
