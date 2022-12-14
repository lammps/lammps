/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_MACROS_HPP_
#define DESUL_ATOMICS_MACROS_HPP_

// Macros

#if (!defined(__CUDA_ARCH__) || !defined(__NVCC__)) &&                       \
    (!defined(__HIP_DEVICE_COMPILE) || !defined(__HIP_PLATFORM_HCC__)) &&    \
    !defined(__SYCL_DEVICE_ONLY__) && !defined(DESUL_HAVE_OPENMP_ATOMICS) && \
    !defined(DESUL_HAVE_SERIAL_ATOMICS)
#define DESUL_IMPL_HAVE_GCC_OR_MSVC_ATOMICS
#endif

// ONLY use GNUC atomics if not compiling for the device
// and we didn't explicitly say to use OPENMP or SERIAL atomics
#if defined(__GNUC__) && defined(DESUL_IMPL_HAVE_GCC_OR_MSVC_ATOMICS)
#define DESUL_HAVE_GCC_ATOMICS
#endif

// Equivalent to above: if we are compiling for the device we
// need to use CUDA/HIP/SYCL atomics instead of MSVC atomics
#if defined(_MSC_VER) && defined(DESUL_IMPL_HAVE_GCC_OR_MSVC_ATOMICS)
#define DESUL_HAVE_MSVC_ATOMICS
#endif

#undef DESUL_IMPL_HAVE_GCC_OR_MSVC_ATOMICS

#ifdef __CUDACC__
#define DESUL_HAVE_CUDA_ATOMICS
#endif

#ifdef __HIPCC__
#define DESUL_HAVE_HIP_ATOMICS
#endif

#ifdef __SYCL_DEVICE_ONLY__
#define DESUL_HAVE_SYCL_ATOMICS
#endif

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || \
    defined(__SYCL_DEVICE_ONLY__)
#define DESUL_HAVE_GPU_LIKE_PROGRESS
#endif

#if defined(DESUL_HAVE_CUDA_ATOMICS) || defined(DESUL_HAVE_HIP_ATOMICS)
#define DESUL_FORCEINLINE_FUNCTION inline __host__ __device__
#define DESUL_INLINE_FUNCTION inline __host__ __device__
#define DESUL_FUNCTION __host__ __device__
#else
#define DESUL_FORCEINLINE_FUNCTION inline
#define DESUL_INLINE_FUNCTION inline
#define DESUL_FUNCTION
#endif

#if !defined(DESUL_HAVE_GPU_LIKE_PROGRESS)
#define DESUL_HAVE_FORWARD_PROGRESS
#endif

#endif  // DESUL_ATOMICS_MACROS_HPP_
