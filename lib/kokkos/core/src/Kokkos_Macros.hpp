// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_MACROS_HPP
#define KOKKOS_MACROS_HPP

//----------------------------------------------------------------------------
/** Pick up configure / build options via #define macros:
 *
 *  KOKKOS_ENABLE_CUDA                Kokkos::Cuda execution and memory spaces
 *  KOKKOS_ENABLE_THREADS             Kokkos::Threads execution space
 *  KOKKOS_ENABLE_HPX                 Kokkos::Experimental::HPX execution space
 *  KOKKOS_ENABLE_OPENMP              Kokkos::OpenMP execution space
 *  KOKKOS_ENABLE_OPENMPTARGET        Kokkos::Experimental::OpenMPTarget
 *                                    execution space
 *  KOKKOS_ENABLE_HIP                 Kokkos::HIP execution space
 *  KOKKOS_ENABLE_SYCL                Kokkos::SYCL execution space
 *  KOKKOS_ENABLE_HWLOC               HWLOC library is available.
 *  KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK  Insert array bounds checks, is expensive!
 *  KOKKOS_ENABLE_CUDA_UVM            Use CUDA UVM for Cuda memory space.
 */

#define KOKKOS_VERSION_LESS(MAJOR, MINOR, PATCH) \
  (KOKKOS_VERSION < ((MAJOR)*10000 + (MINOR)*100 + (PATCH)))

#define KOKKOS_VERSION_LESS_EQUAL(MAJOR, MINOR, PATCH) \
  (KOKKOS_VERSION <= ((MAJOR)*10000 + (MINOR)*100 + (PATCH)))

#define KOKKOS_VERSION_GREATER(MAJOR, MINOR, PATCH) \
  (KOKKOS_VERSION > ((MAJOR)*10000 + (MINOR)*100 + (PATCH)))

#define KOKKOS_VERSION_GREATER_EQUAL(MAJOR, MINOR, PATCH) \
  (KOKKOS_VERSION >= ((MAJOR)*10000 + (MINOR)*100 + (PATCH)))

#define KOKKOS_VERSION_EQUAL(MAJOR, MINOR, PATCH) \
  (KOKKOS_VERSION == ((MAJOR)*10000 + (MINOR)*100 + (PATCH)))

#if !KOKKOS_VERSION_EQUAL(KOKKOS_VERSION_MAJOR, KOKKOS_VERSION_MINOR, \
                          KOKKOS_VERSION_PATCH)
#error implementation bug
#endif

#ifndef KOKKOS_DONT_INCLUDE_CORE_CONFIG_H
#include <KokkosCore_config.h>
#include <impl/Kokkos_DesulAtomicsConfig.hpp>
#include <impl/Kokkos_NvidiaGpuArchitectures.hpp>
#endif

#if __has_include(<version>)
#include <version>
#else
#include <ciso646>
#endif

//----------------------------------------------------------------------------
/** Pick up compiler specific #define macros:
 *
 *  Macros for known compilers evaluate to an integral version value
 *
 *  KOKKOS_COMPILER_NVCC
 *  KOKKOS_COMPILER_GNU
 *  KOKKOS_COMPILER_INTEL_LLVM
 *  KOKKOS_COMPILER_CRAYC
 *  KOKKOS_COMPILER_APPLECC
 *  KOKKOS_COMPILER_CLANG
 *  KOKKOS_COMPILER_NVHPC
 *  KOKKOS_COMPILER_MSVC
 *
 *  A suite of 'KOKKOS_ENABLE_PRAGMA_...' are defined for internal use.
 *
 *  Macros for marking functions to run in an execution space:
 *
 *  KOKKOS_FUNCTION
 *  KOKKOS_INLINE_FUNCTION        request compiler to inline
 *  KOKKOS_FORCEINLINE_FUNCTION   force compiler to inline, use with care!
 */

//----------------------------------------------------------------------------

#if defined(KOKKOS_ENABLE_ATOMICS_BYPASS) &&                              \
    (defined(KOKKOS_ENABLE_THREADS) || defined(KOKKOS_ENABLE_CUDA) ||     \
     defined(KOKKOS_ENABLE_OPENMP) || defined(KOKKOS_ENABLE_HPX) ||       \
     defined(KOKKOS_ENABLE_OPENMPTARGET) || defined(KOKKOS_ENABLE_HIP) || \
     defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENACC))
#error Atomics may only be disabled if neither a host parallel nor a device backend is enabled
#endif

#define KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA

#include <KokkosCore_Config_SetupBackend.hpp>

//----------------------------------------------------------------------------
// Mapping compiler built-ins to KOKKOS_COMPILER_*** macros

#if defined(__NVCC__)
// NVIDIA compiler is being used.
// Code is parsed and separated into host and device code.
// Host code is compiled again with another compiler.
// Device code is compile to 'ptx'.
// NOTE: There is no __CUDACC_VER_PATCH__ officially, its __CUDACC_VER_BUILD__
// which does have more than one digit (potentially undefined number of them).
// This macro definition is in line with our other compiler defs
#define KOKKOS_COMPILER_NVCC \
  __CUDACC_VER_MAJOR__ * 100 + __CUDACC_VER_MINOR__ * 10
#endif  // #if defined( __NVCC__ )

#if !defined(KOKKOS_LAMBDA)
#define KOKKOS_LAMBDA [=]
#endif

#if !defined(KOKKOS_CLASS_LAMBDA)
#define KOKKOS_CLASS_LAMBDA [ =, *this ]
#endif

// #if !defined( __CUDA_ARCH__ ) // Not compiling Cuda code to 'ptx'.

// Intel compiler for host code.

#if defined(__INTEL_LLVM_COMPILER)
#define KOKKOS_COMPILER_INTEL_LLVM __INTEL_LLVM_COMPILER

// Cray compiler for device offload code
#elif defined(__cray__) && defined(__clang__)
#define KOKKOS_COMPILER_CRAY_LLVM \
  __cray_major__ * 100 + __cray_minor__ * 10 + __cray_patchlevel__

#elif defined(_CRAYC)
// CRAY compiler for host code
#define KOKKOS_COMPILER_CRAYC _CRAYC

#elif defined(__APPLE_CC__)
#define KOKKOS_COMPILER_APPLECC __APPLE_CC__

#elif defined(__NVCOMPILER)
#define KOKKOS_COMPILER_NVHPC                                 \
  __NVCOMPILER_MAJOR__ * 10000 + __NVCOMPILER_MINOR__ * 100 + \
      __NVCOMPILER_PATCHLEVEL__

#elif defined(__clang__)
// Check this after the Clang-based proprietary compilers which will also define
// __clang__
#define KOKKOS_COMPILER_CLANG \
  __clang_major__ * 100 + __clang_minor__ * 10 + __clang_patchlevel__

#elif defined(__GNUC__)
// Check this here because many compilers (at least Clang variants and Intel
// classic) define `__GNUC__` for compatibility
#define KOKKOS_COMPILER_GNU \
  __GNUC__ * 100 + __GNUC_MINOR__ * 10 + __GNUC_PATCHLEVEL__

#if (1040 > KOKKOS_COMPILER_GNU)
#error "Compiling with GCC version earlier than 10.4.0 is not supported."
#endif

#elif defined(_MSC_VER)
// Check this after Intel and Clang because those define _MSC_VER for
// compatibility
#define KOKKOS_COMPILER_MSVC _MSC_VER
#endif

#if defined(_OPENMP)
//  Compiling with OpenMP.
//  The value of _OPENMP is an integer value YYYYMM
//  where YYYY and MM are the year and month designation
//  of the supported OpenMP API version.
#endif  // #if defined( _OPENMP )

//----------------------------------------------------------------------------
// Intel compiler macros

#if defined(KOKKOS_COMPILER_INTEL_LLVM)
#if KOKKOS_COMPILER_INTEL_LLVM >= 20230100
#define KOKKOS_ENABLE_PRAGMA_UNROLL 1

#ifndef __SYCL_DEVICE_ONLY__
#define KOKKOS_ENABLE_PRAGMA_LOOPCOUNT 1
#define KOKKOS_ENABLE_PRAGMA_VECTOR 1
#define KOKKOS_ENABLE_PRAGMA_IVDEP 1
#endif
#endif

#if defined(_WIN32)
#define KOKKOS_RESTRICT __restrict
#else
#define KOKKOS_RESTRICT __restrict__
#endif

#ifndef KOKKOS_IMPL_ALIGN_PTR
#if defined(_WIN32)
#define KOKKOS_IMPL_ALIGN_PTR(size)
#else
#define KOKKOS_IMPL_ALIGN_PTR(size) __attribute__((align_value(size)))
#endif
#endif

#if !defined(KOKKOS_ENABLE_ASM) && !defined(_WIN32)
#define KOKKOS_ENABLE_ASM 1
#endif

#if !defined(KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION)
#if !defined(_WIN32)
#define KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION \
  inline __attribute__((always_inline))
#define KOKKOS_IMPL_HOST_FORCEINLINE __attribute__((always_inline))
#else
#define KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION inline
#endif
#endif

#if defined(__MIC__)
// Compiling for Xeon Phi
#endif
#endif

//----------------------------------------------------------------------------
// Cray compiler macros

#if defined(KOKKOS_COMPILER_CRAYC)
#endif

//----------------------------------------------------------------------------
// CLANG compiler macros

#if defined(KOKKOS_COMPILER_CLANG)
// #define KOKKOS_ENABLE_PRAGMA_UNROLL 1
// #define KOKKOS_ENABLE_PRAGMA_IVDEP 1
// #define KOKKOS_ENABLE_PRAGMA_LOOPCOUNT 1
// #define KOKKOS_ENABLE_PRAGMA_VECTOR 1

#if !defined(KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION)
#define KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION \
  inline __attribute__((always_inline))
#define KOKKOS_IMPL_HOST_FORCEINLINE __attribute__((always_inline))
#endif

#if !defined(KOKKOS_IMPL_ALIGN_PTR)
#define KOKKOS_IMPL_ALIGN_PTR(size) __attribute__((aligned(size)))
#endif

#endif

//----------------------------------------------------------------------------
// GNU Compiler macros

#if defined(KOKKOS_COMPILER_GNU)
// #define KOKKOS_ENABLE_PRAGMA_UNROLL 1
// #define KOKKOS_ENABLE_PRAGMA_IVDEP 1
// #define KOKKOS_ENABLE_PRAGMA_LOOPCOUNT 1
// #define KOKKOS_ENABLE_PRAGMA_VECTOR 1

#if !defined(KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION)
#define KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION \
  inline __attribute__((always_inline))
#define KOKKOS_IMPL_HOST_FORCEINLINE __attribute__((always_inline))
#endif

#define KOKKOS_RESTRICT __restrict__

#if !defined(KOKKOS_ENABLE_ASM) && !defined(__PGIC__) &&            \
    (defined(__amd64) || defined(__amd64__) || defined(__x86_64) || \
     defined(__x86_64__) || defined(__PPC64__))
#define KOKKOS_ENABLE_ASM 1
#endif
#endif

//----------------------------------------------------------------------------

#if defined(KOKKOS_COMPILER_NVHPC)
#define KOKKOS_ENABLE_PRAGMA_UNROLL 1
#define KOKKOS_ENABLE_PRAGMA_IVDEP 1
// #define KOKKOS_ENABLE_PRAGMA_LOOPCOUNT 1
#define KOKKOS_ENABLE_PRAGMA_VECTOR 1
#endif

//----------------------------------------------------------------------------

#if defined(KOKKOS_COMPILER_NVCC)
#if defined(__CUDA_ARCH__)
#define KOKKOS_ENABLE_PRAGMA_UNROLL 1
#endif
#endif

//----------------------------------------------------------------------------
// Define function marking macros if compiler specific macros are undefined:

#if !defined(KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION)
#define KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION inline
#endif

#if !defined(KOKKOS_IMPL_HOST_FORCEINLINE)
#define KOKKOS_IMPL_HOST_FORCEINLINE inline
#endif

#if !defined(KOKKOS_IMPL_FORCEINLINE_FUNCTION)
#define KOKKOS_IMPL_FORCEINLINE_FUNCTION KOKKOS_IMPL_HOST_FORCEINLINE_FUNCTION
#endif

#if !defined(KOKKOS_IMPL_FORCEINLINE)
#define KOKKOS_IMPL_FORCEINLINE KOKKOS_IMPL_HOST_FORCEINLINE
#endif

#if !defined(KOKKOS_IMPL_INLINE_FUNCTION)
#define KOKKOS_IMPL_INLINE_FUNCTION inline
#endif

#if !defined(KOKKOS_IMPL_FUNCTION)
#define KOKKOS_IMPL_FUNCTION /**/
#endif

#if !defined(KOKKOS_INLINE_FUNCTION_DELETED)
#define KOKKOS_INLINE_FUNCTION_DELETED
#endif

#if !defined(KOKKOS_DEFAULTED_FUNCTION)
#define KOKKOS_DEFAULTED_FUNCTION
#endif

#if !defined(KOKKOS_DEDUCTION_GUIDE)
#define KOKKOS_DEDUCTION_GUIDE
#endif

#if !defined(KOKKOS_IMPL_HOST_FUNCTION)
#define KOKKOS_IMPL_HOST_FUNCTION
#endif

#if !defined(KOKKOS_IMPL_DEVICE_FUNCTION)
#define KOKKOS_IMPL_DEVICE_FUNCTION
#endif

// FIXME_OPENACC FIXME_OPENMPTARGET
// Move to setup files once there is more content
// clang-format off
#if defined(KOKKOS_ENABLE_OPENACC)
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION @"KOKKOS_RELOCATABLE_FUNCTION is not supported for the OpenACC backend"
#endif
#if defined(KOKKOS_ENABLE_OPENMPTARGET)
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION @"KOKKOS_RELOCATABLE_FUNCTION is not supported for the OpenMPTarget backend"
#endif
// clang-format on

#if !defined(KOKKOS_IMPL_RELOCATABLE_FUNCTION)
#define KOKKOS_IMPL_RELOCATABLE_FUNCTION
#endif

//----------------------------------------------------------------------------
// Define final version of functions. This is so that clang tidy can find these
// macros more easily
#if defined(__clang_analyzer__)
#define KOKKOS_FUNCTION \
  KOKKOS_IMPL_FUNCTION __attribute__((annotate("KOKKOS_FUNCTION")))
#define KOKKOS_INLINE_FUNCTION \
  KOKKOS_IMPL_INLINE_FUNCTION  \
  __attribute__((annotate("KOKKOS_INLINE_FUNCTION")))
#define KOKKOS_FORCEINLINE_FUNCTION \
  KOKKOS_IMPL_FORCEINLINE_FUNCTION  \
  __attribute__((annotate("KOKKOS_FORCEINLINE_FUNCTION")))
#define KOKKOS_RELOCATABLE_FUNCTION \
  KOKKOS_IMPL_RELOCATABLE_FUNCTION  \
  __attribute__((annotate("KOKKOS_RELOCATABLE_FUNCTION")))
#else
#define KOKKOS_FUNCTION KOKKOS_IMPL_FUNCTION
#define KOKKOS_INLINE_FUNCTION KOKKOS_IMPL_INLINE_FUNCTION
#define KOKKOS_FORCEINLINE_FUNCTION KOKKOS_IMPL_FORCEINLINE_FUNCTION
#define KOKKOS_RELOCATABLE_FUNCTION KOKKOS_IMPL_RELOCATABLE_FUNCTION
#endif

//----------------------------------------------------------------------------
// Define empty macro for restrict if necessary:

#if !defined(KOKKOS_RESTRICT)
#define KOKKOS_RESTRICT
#endif

//----------------------------------------------------------------------------
// Define Macro for alignment:

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#ifndef KOKKOS_MEMORY_ALIGNMENT
#define KOKKOS_IMPL_MEMORY_ALIGNMENT 64
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
#define KOKKOS_MEMORY_ALIGNMENT                                         \
  [] {                                                                  \
    int memory_alignment                                                \
        [[deprecated("KOKKOS_MEMORY_ALIGNMENT macro is deprecated")]] = \
            KOKKOS_IMPL_MEMORY_ALIGNMENT;                               \
    return memory_alignment;                                            \
  }();
#else
#define KOKKOS_MEMORY_ALIGNMENT KOKKOS_IMPL_MEMORY_ALIGNMENT
#endif
#else
#define KOKKOS_IMPL_MEMORY_ALIGNMENT KOKKOS_MEMORY_ALIGNMENT
#endif
#else  // KOKKOS_ENABLE_DEPRECATED_CODE_4
#ifdef KOKKOS_MEMORY_ALIGNMENT
static_assert(false,
              "External definition of KOKKOS_MEMORY_ALIGNMENT is not allowed");
#endif
#endif

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#ifndef KOKKOS_MEMORY_ALIGNMENT_THRESHOLD
#define KOKKOS_IMPL_MEMORY_ALIGNMENT_THRESHOLD 1
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
#define KOKKOS_MEMORY_ALIGNMENT_THRESHOLD                            \
  [] {                                                               \
    int memory_alignment [[deprecated(                               \
        "KOKKOS_MEMORY_ALIGNMENT_THRESHOLD macro is deprecated")]] = \
        KOKKOS_IMPL_MEMORY_ALIGNMENT_THRESHOLD;                      \
    return memory_alignment;                                         \
  }();
#else
#define KOKKOS_MEMORY_ALIGNMENT_THRESHOLD KOKKOS_IMPL_MEMORY_ALIGNMENT_THRESHOLD
#endif
#else
#define KOKKOS_IMPL_MEMORY_ALIGNMENT_THRESHOLD KOKKOS_MEMORY_ALIGNMENT_THRESHOLD
#endif
#else  // KOKKOS_ENABLE_DEPRECATED_CODE_4
#ifdef KOKKOS_MEMORY_ALIGNMENT_THRESHOLD
static_assert(
    false,
    "External definition of KOKKOS_MEMORY_ALIGNMENT_THRESHOLD is not allowed");
#endif
#endif

#if !defined(KOKKOS_IMPL_ALIGN_PTR)
#define KOKKOS_IMPL_ALIGN_PTR(size) /* */
#endif

//----------------------------------------------------------------------------
// Determine the default execution space for parallel dispatch.
// There is zero or one default execution space specified.

#if 1 < ((defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA) ? 1 : 0) +         \
         (defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HIP) ? 1 : 0) +          \
         (defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SYCL) ? 1 : 0) +         \
         (defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENACC) ? 1 : 0) +      \
         (defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMPTARGET) ? 1 : 0) + \
         (defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP) ? 1 : 0) +       \
         (defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS) ? 1 : 0) +      \
         (defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HPX) ? 1 : 0) +          \
         (defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL) ? 1 : 0))
#error "More than one KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_* specified."
#endif

// If default is not specified then chose from enabled execution spaces.
// Priority: CUDA, HIP, SYCL, OPENACC, OPENMPTARGET, OPENMP, THREADS, HPX,
// SERIAL
#if defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA)
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HIP)
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SYCL)
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENACC)
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMPTARGET)
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP)
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS)
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HPX)
#elif defined(KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL)
#elif defined(KOKKOS_ENABLE_CUDA)
#define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA
#elif defined(KOKKOS_ENABLE_HIP)
#define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HIP
#elif defined(KOKKOS_ENABLE_SYCL)
#define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SYCL
#elif defined(KOKKOS_ENABLE_OPENACC)
#define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENACC
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
#define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMPTARGET
#elif defined(KOKKOS_ENABLE_OPENMP)
#define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP
#elif defined(KOKKOS_ENABLE_THREADS)
#define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS
#elif defined(KOKKOS_ENABLE_HPX)
#define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_HPX
#else
#define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL
#endif

//----------------------------------------------------------------------------

// Remove surrounding parentheses if present
#define KOKKOS_IMPL_STRIP_PARENS(X) KOKKOS_IMPL_ESC(KOKKOS_IMPL_ISH X)
#define KOKKOS_IMPL_ISH(...) KOKKOS_IMPL_ISH __VA_ARGS__
#define KOKKOS_IMPL_ESC(...) KOKKOS_IMPL_ESC_(__VA_ARGS__)
#define KOKKOS_IMPL_ESC_(...) KOKKOS_IMPL_VAN_##__VA_ARGS__
#define KOKKOS_IMPL_VAN_KOKKOS_IMPL_ISH

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_NVHPC)
#include <nv/target>
#define KOKKOS_IF_ON_DEVICE(CODE) NV_IF_TARGET(NV_IS_DEVICE, CODE)
#define KOKKOS_IF_ON_HOST(CODE) NV_IF_TARGET(NV_IS_HOST, CODE)
#endif

#ifdef KOKKOS_ENABLE_OPENMPTARGET
#ifdef KOKKOS_COMPILER_NVHPC
#define KOKKOS_IF_ON_DEVICE(CODE)   \
  if (__builtin_is_device_code()) { \
    KOKKOS_IMPL_STRIP_PARENS(CODE)  \
  }
#define KOKKOS_IF_ON_HOST(CODE)      \
  if (!__builtin_is_device_code()) { \
    KOKKOS_IMPL_STRIP_PARENS(CODE)   \
  }
#else
// Base function.
static constexpr bool kokkos_omp_on_host() { return true; }

#pragma omp begin declare variant match(device = {kind(host)})
static constexpr bool kokkos_omp_on_host() { return true; }
#pragma omp end declare variant

#pragma omp begin declare variant match(device = {kind(nohost)})
static constexpr bool kokkos_omp_on_host() { return false; }
#pragma omp end declare variant

#define KOKKOS_IF_ON_DEVICE(CODE)        \
  if constexpr (!kokkos_omp_on_host()) { \
    KOKKOS_IMPL_STRIP_PARENS(CODE)       \
  }
#define KOKKOS_IF_ON_HOST(CODE)         \
  if constexpr (kokkos_omp_on_host()) { \
    KOKKOS_IMPL_STRIP_PARENS(CODE)      \
  }
#endif
#endif

#ifdef KOKKOS_ENABLE_OPENACC
#ifdef KOKKOS_COMPILER_NVHPC
#define KOKKOS_IF_ON_DEVICE(CODE)   \
  if (__builtin_is_device_code()) { \
    KOKKOS_IMPL_STRIP_PARENS(CODE)  \
  }
#define KOKKOS_IF_ON_HOST(CODE)      \
  if (!__builtin_is_device_code()) { \
    KOKKOS_IMPL_STRIP_PARENS(CODE)   \
  }
#else
#include <openacc.h>
// FIXME_OPENACC acc_on_device is a non-constexpr function
#define KOKKOS_IF_ON_DEVICE(CODE)                     \
  if constexpr (acc_on_device(acc_device_not_host)) { \
    KOKKOS_IMPL_STRIP_PARENS(CODE)                    \
  }
#define KOKKOS_IF_ON_HOST(CODE)                   \
  if constexpr (acc_on_device(acc_device_host)) { \
    KOKKOS_IMPL_STRIP_PARENS(CODE)                \
  }
#endif
#endif

#if !defined(KOKKOS_IF_ON_HOST) && !defined(KOKKOS_IF_ON_DEVICE)
#if (defined(KOKKOS_ENABLE_CUDA) && defined(__CUDA_ARCH__)) ||         \
    (defined(KOKKOS_ENABLE_HIP) && defined(__HIP_DEVICE_COMPILE__)) || \
    (defined(KOKKOS_ENABLE_SYCL) && defined(__SYCL_DEVICE_ONLY__))
#define KOKKOS_IF_ON_DEVICE(CODE) \
  { KOKKOS_IMPL_STRIP_PARENS(CODE) }
#define KOKKOS_IF_ON_HOST(CODE) \
  {}
#else
#define KOKKOS_IF_ON_DEVICE(CODE) \
  {}
#define KOKKOS_IF_ON_HOST(CODE) \
  { KOKKOS_IMPL_STRIP_PARENS(CODE) }
#endif
#endif

//----------------------------------------------------------------------------

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_DEPRECATED_CODE_4)
#define KOKKOS_ENABLE_CUDA_LDG_INTRINSIC
#endif

#define KOKKOS_INVALID_INDEX (~std::size_t(0))

#define KOKKOS_IMPL_CTOR_DEFAULT_ARG KOKKOS_INVALID_INDEX

// Guard intel compiler version 19 and older
// intel error #2651: attribute does not apply to any entity
// using <deprecated_type> KOKKOS_DEPRECATED = ...
#if defined(KOKKOS_ENABLE_DEPRECATION_WARNINGS) && !defined(__NVCC__)
#define KOKKOS_DEPRECATED [[deprecated]]
#define KOKKOS_DEPRECATED_WITH_COMMENT(comment) [[deprecated(comment)]]
#else
#define KOKKOS_DEPRECATED
#define KOKKOS_DEPRECATED_WITH_COMMENT(comment)
#endif

#define KOKKOS_IMPL_STRINGIFY(x) #x
#define KOKKOS_IMPL_TOSTRING(x) KOKKOS_IMPL_STRINGIFY(x)

#ifdef _MSC_VER
#define KOKKOS_IMPL_DO_PRAGMA(x) __pragma(x)
#define KOKKOS_IMPL_WARNING(desc) \
  KOKKOS_IMPL_DO_PRAGMA(message(  \
      __FILE__ "(" KOKKOS_IMPL_TOSTRING(__LINE__) ") : warning: " #desc))
#else
#define KOKKOS_IMPL_DO_PRAGMA(x) _Pragma(#x)
#define KOKKOS_IMPL_WARNING(desc) KOKKOS_IMPL_DO_PRAGMA(message(#desc))
#endif

// clang-format off
#if defined(__NVCOMPILER)
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH() \
    _Pragma("diag_suppress 1216") \
    _Pragma("diag_suppress deprecated_entity_with_custom_message")
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP() \
    _Pragma("diag_default 1216") \
    _Pragma("diag_suppress deprecated_entity_with_custom_message")
#elif defined(__EDG__)
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH() \
    _Pragma("warning push")                              \
    _Pragma("warning disable 1478")
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP() \
    _Pragma("warning pop")
#elif defined(__GNUC__) || defined(__clang__)
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH() \
    _Pragma("GCC diagnostic push")                       \
    _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP() \
    _Pragma("GCC diagnostic pop")
#elif defined(_MSC_VER)
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH() \
    _Pragma("warning(push)")                             \
    _Pragma("warning(disable: 4996)")
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP() \
    _Pragma("warning(pop)")
#else
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_PUSH()
  #define KOKKOS_IMPL_DISABLE_DEPRECATED_WARNINGS_POP()
#endif

#if defined(__NVCOMPILER)
#define KOKKOS_IMPL_DISABLE_UNREACHABLE_WARNINGS_PUSH() \
  _Pragma("diag_suppress code_is_unreachable")          \
  _Pragma("diag_suppress initialization_not_reachable")
#define KOKKOS_IMPL_DISABLE_UNREACHABLE_WARNINGS_POP() \
  _Pragma("diag_default code_is_unreachable")          \
  _Pragma("diag_default initialization_not_reachable")
#else
#define KOKKOS_IMPL_DISABLE_UNREACHABLE_WARNINGS_PUSH()
#define KOKKOS_IMPL_DISABLE_UNREACHABLE_WARNINGS_POP()
#endif
// clang-format on

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#define KOKKOS_ATTRIBUTE_NODISCARD [[nodiscard]]
#endif

#if (defined(KOKKOS_COMPILER_GNU) || defined(KOKKOS_COMPILER_CLANG) ||         \
     defined(KOKKOS_COMPILER_INTEL_LLVM) || defined(KOKKOS_COMPILER_NVHPC)) && \
    !defined(_WIN32) && !defined(__ANDROID__)
#if __has_include(<execinfo.h>)
#define KOKKOS_IMPL_ENABLE_STACKTRACE
#endif
#define KOKKOS_IMPL_ENABLE_CXXABI
#endif

#if (defined(KOKKOS_IMPL_WINDOWS_CUDA) || defined(KOKKOS_COMPILER_MSVC)) && \
    !defined(KOKKOS_COMPILER_CLANG)
// MSVC (as of 16.5.5 at least) does not do empty base class optimization by
// default when there are multiple bases, even though the standard requires it
// for standard layout types.
#define KOKKOS_IMPL_ENFORCE_EMPTY_BASE_OPTIMIZATION __declspec(empty_bases)
#else
#define KOKKOS_IMPL_ENFORCE_EMPTY_BASE_OPTIMIZATION
#endif

#if defined(KOKKOS_IMPL_BUILD_SHARED_LIBS) && defined(_WIN32)
#ifdef KOKKOS_IMPL_EXPORT_SYMBOLS
#define KOKKOS_IMPL_EXPORT __declspec(dllexport)
#else
#define KOKKOS_IMPL_EXPORT __declspec(dllimport)
#endif
#else
#define KOKKOS_IMPL_EXPORT
#endif

#endif  // #ifndef KOKKOS_MACROS_HPP
