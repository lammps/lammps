// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CUDA_SETUP_HPP_
#define KOKKOS_CUDA_SETUP_HPP_

#if !defined(KOKKOS_ENABLE_CUDA)
#error \
    "KOKKOS_ENABLE_CUDA was not defined, but Kokkos_Setup_Cuda.hpp was included anyway."
#endif

#if defined(KOKKOS_ENABLE_CUDA) && !defined(__CUDACC__)
#error \
    "KOKKOS_ENABLE_CUDA defined but the compiler is not defining the __CUDACC__ macro as expected"
// Some tooling environments will still function better if we do this here.
#define __CUDACC__
#endif /* defined(KOKKOS_ENABLE_CUDA) && !defined(__CUDACC__) */

// Compiling with a CUDA compiler.
//
//  Include <cuda.h> to pick up the CUDA_VERSION macro defined as:
//    CUDA_VERSION = ( MAJOR_VERSION * 1000 ) + ( MINOR_VERSION * 10 )
//
//  When generating device code the __CUDA_ARCH__ macro is defined as:
//    __CUDA_ARCH__ = ( MAJOR_CAPABILITY * 100 ) + ( MINOR_CAPABILITY * 10 )

#include <cuda_runtime.h>
#include <cuda.h>

#if defined(_WIN32)
#define KOKKOS_IMPL_WINDOWS_CUDA
#endif

#if !defined(CUDA_VERSION)
#error "#include <cuda.h> did not define CUDA_VERSION."
#endif

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 300)
// Compiling with CUDA compiler for device code.
#error "Cuda device capability >= 3.0 is required."
#endif

#define KOKKOS_LAMBDA [=] __host__ __device__
#define KOKKOS_CLASS_LAMBDA [ =, *this ] __host__ __device__

#define KOKKOS_DEDUCTION_GUIDE __host__ __device__

#define KOKKOS_IMPL_FORCEINLINE_FUNCTION __device__ __host__ __forceinline__
#define KOKKOS_IMPL_FORCEINLINE __forceinline__
#define KOKKOS_IMPL_INLINE_FUNCTION __device__ __host__ inline
#define KOKKOS_IMPL_FUNCTION __device__ __host__
#define KOKKOS_IMPL_HOST_FUNCTION __host__
#define KOKKOS_IMPL_DEVICE_FUNCTION __device__

// clang-format off
#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION __device__ __host__
#else
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION @"KOKKOS_RELOCATABLE_FUNCTION requires Kokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON"
#endif
// clang-format on

#endif /* KOKKOS_CUDA_SETUP_HPP_ */
