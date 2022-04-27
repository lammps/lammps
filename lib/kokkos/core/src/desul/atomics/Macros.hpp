/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_MACROS_HPP_
#define DESUL_ATOMICS_MACROS_HPP_

// Macros

#if defined(__GNUC__) && \
    (!defined(__CUDA_ARCH__) || !defined(__NVCC__)) && \
    (!defined(__HIP_DEVICE_COMPILE) || !defined(__HIP_PLATFORM_HCC__)) && \
    !defined(__SYCL_DEVICE_ONLY__) && \
    !defined(DESUL_HAVE_OPENMP_ATOMICS) && \
    !defined(DESUL_HAVE_SERIAL_ATOMICS)
#define DESUL_HAVE_GCC_ATOMICS
#endif

#ifdef _MSC_VER
#define DESUL_HAVE_MSVC_ATOMICS
#endif

#ifdef __CUDACC__
#define DESUL_HAVE_CUDA_ATOMICS
#endif

#ifdef __HIPCC__
#define DESUL_HAVE_HIP_ATOMICS
#endif

#ifdef __SYCL_DEVICE_ONLY__
#define DESUL_HAVE_SYCL_ATOMICS
#ifdef __clang__
#define DESUL_SYCL_NAMESPACE sycl::ONEAPI
#else
#define DESUL_SYCL_NAMESPACE sycl
#endif
#endif

#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__) || defined(__SYCL_DEVICE_ONLY__)
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
