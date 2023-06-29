/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_MACROS_HPP_
#define DESUL_ATOMICS_MACROS_HPP_

#include <desul/atomics/Config.hpp>

// Intercept incompatible relocatable device code mode which leads to ODR violations
#ifdef DESUL_ATOMICS_ENABLE_CUDA
#if (defined(__clang__) && defined(__CUDA__) && defined(__CLANG_RDC__)) || \
    defined(__CUDACC_RDC__)
#define DESUL_IMPL_CUDA_RDC
#endif

#if (defined(DESUL_ATOMICS_ENABLE_CUDA_SEPARABLE_COMPILATION) &&  \
     !defined(DESUL_IMPL_CUDA_RDC)) ||                            \
    (!defined(DESUL_ATOMICS_ENABLE_CUDA_SEPARABLE_COMPILATION) && \
     defined(DESUL_IMPL_CUDA_RDC))
#error Relocatable device code mode incompatible with desul atomics configuration
#endif

#ifdef DESUL_IMPL_CUDA_RDC
#undef DESUL_IMPL_CUDA_RDC
#endif
#endif

#ifdef DESUL_ATOMICS_ENABLE_HIP
#if (defined(DESUL_ATOMICS_ENABLE_HIP_SEPARABLE_COMPILATION) &&  \
     !defined(__CLANG_RDC__)) ||                                 \
    (!defined(DESUL_ATOMICS_ENABLE_HIP_SEPARABLE_COMPILATION) && \
     defined(__CLANG_RDC__))
#error Relocatable device code mode incompatible with desul atomics configuration
#endif
#endif

// Macros

#if defined(DESUL_ATOMICS_ENABLE_CUDA) && defined(__CUDACC__)
#define DESUL_HAVE_CUDA_ATOMICS
#endif

#if defined(DESUL_ATOMICS_ENABLE_HIP) && defined(__HIPCC__)
#define DESUL_HAVE_HIP_ATOMICS
#endif

#if defined(DESUL_ATOMICS_ENABLE_SYCL) && defined(SYCL_LANGUAGE_VERSION)
#define DESUL_HAVE_SYCL_ATOMICS
#endif

#if defined(DESUL_ATOMICS_ENABLE_OPENMP)
#define DESUL_HAVE_OPENMP_ATOMICS
#endif

// ONLY use GNUC atomics if not explicitly say to use OpenMP atomics
#if !defined(DESUL_HAVE_OPENMP_ATOMICS) && defined(__GNUC__)
#define DESUL_HAVE_GCC_ATOMICS
#endif

// Equivalent to above for MSVC atomics
#if !defined(DESUL_HAVE_OPENMP_ATOMICS) && defined(_MSC_VER)
#define DESUL_HAVE_MSVC_ATOMICS
#endif

#if defined(DESUL_HAVE_CUDA_ATOMICS) || defined(DESUL_HAVE_HIP_ATOMICS)
#define DESUL_FORCEINLINE_FUNCTION inline __host__ __device__
#define DESUL_INLINE_FUNCTION inline __host__ __device__
#define DESUL_FUNCTION __host__ __device__
#define DESUL_IMPL_HOST_FUNCTION __host__
#define DESUL_IMPL_DEVICE_FUNCTION __device__
#else
#define DESUL_FORCEINLINE_FUNCTION inline
#define DESUL_INLINE_FUNCTION inline
#define DESUL_FUNCTION
#define DESUL_IMPL_HOST_FUNCTION
#define DESUL_IMPL_DEVICE_FUNCTION
#endif

#define DESUL_IMPL_STRIP_PARENS(X) DESUL_IMPL_ESC(DESUL_IMPL_ISH X)
#define DESUL_IMPL_ISH(...) DESUL_IMPL_ISH __VA_ARGS__
#define DESUL_IMPL_ESC(...) DESUL_IMPL_ESC_(__VA_ARGS__)
#define DESUL_IMPL_ESC_(...) DESUL_IMPL_VAN_##__VA_ARGS__
#define DESUL_IMPL_VAN_DESUL_IMPL_ISH

#if (defined(DESUL_ATOMICS_ENABLE_CUDA) && defined(__CUDACC__)) && defined(__NVCOMPILER)
#include <nv/target>
#define DESUL_IF_ON_DEVICE(CODE) NV_IF_TARGET(NV_IS_DEVICE, CODE)
#define DESUL_IF_ON_HOST(CODE) NV_IF_TARGET(NV_IS_HOST, CODE)
#endif

// FIXME OpenMP Offload differentiate between device and host, but do we need this?
#if defined(DESUL_HAVE_OPENMP_ATOMICS)
#if 0
// Base function.
static constexpr bool desul_impl_omp_on_host() { return true; }

#pragma omp begin declare variant match(device = {kind(host)})
static constexpr bool desul_impl_omp_on_host() { return true; }
#pragma omp end declare variant

#pragma omp begin declare variant match(device = {kind(nohost)})
static constexpr bool desul_impl_omp_on_host() { return false; }
#pragma omp end declare variant

#define DESUL_IF_ON_DEVICE(CODE)             \
  if constexpr (!desul_impl_omp_on_host()) { \
    DESUL_IMPL_STRIP_PARENS(CODE)            \
  }
#define DESUL_IF_ON_HOST(CODE)              \
  if constexpr (desul_impl_omp_on_host()) { \
    DESUL_IMPL_STRIP_PARENS(CODE)           \
  }
#else
#define DESUL_IF_ON_DEVICE(CODE) \
  {}
#define DESUL_IF_ON_HOST(CODE) \
  { DESUL_IMPL_STRIP_PARENS(CODE) }
#endif
#endif

#if !defined(DESUL_IF_ON_HOST) && !defined(DESUL_IF_ON_DEVICE)
#if (defined(DESUL_ATOMICS_ENABLE_CUDA) && defined(__CUDA_ARCH__)) ||         \
    (defined(DESUL_ATOMICS_ENABLE_HIP) && defined(__HIP_DEVICE_COMPILE__)) || \
    (defined(DESUL_ATOMICS_ENABLE_SYCL) && defined(__SYCL_DEVICE_ONLY__))
#define DESUL_IF_ON_DEVICE(CODE) \
  { DESUL_IMPL_STRIP_PARENS(CODE) }
#define DESUL_IF_ON_HOST(CODE) \
  {}
#else
#define DESUL_IF_ON_DEVICE(CODE) \
  {}
#define DESUL_IF_ON_HOST(CODE) \
  { DESUL_IMPL_STRIP_PARENS(CODE) }
#endif
#endif

#endif  // DESUL_ATOMICS_MACROS_HPP_
