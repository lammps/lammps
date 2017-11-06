/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_MACROS_HPP
#define KOKKOS_MACROS_HPP

//----------------------------------------------------------------------------
/** Pick up configure / build options via #define macros:
 *
 *  KOKKOS_ENABLE_CUDA                Kokkos::Cuda execution and memory spaces
 *  KOKKOS_ENABLE_THREADS             Kokkos::Threads execution space
 *  KOKKOS_ENABLE_QTHREADS            Kokkos::Qthreads execution space
 *  KOKKOS_ENABLE_OPENMP              Kokkos::OpenMP execution space
 *  KOKKOS_ENABLE_OPENMPTARGET        Kokkos::Experimental::OpenMPTarget execution space
 *  KOKKOS_ENABLE_HWLOC               HWLOC library is available.
 *  KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK  Insert array bounds checks, is expensive!
 *  KOKKOS_ENABLE_MPI                 Negotiate MPI/execution space interactions.
 *  KOKKOS_ENABLE_CUDA_UVM            Use CUDA UVM for Cuda memory space.
 */

#ifndef KOKKOS_DONT_INCLUDE_CORE_CONFIG_H
  #include <KokkosCore_config.h>
#endif

#include <impl/Kokkos_OldMacros.hpp>

//----------------------------------------------------------------------------
/** Pick up compiler specific #define macros:
 *
 *  Macros for known compilers evaluate to an integral version value
 *
 *  KOKKOS_COMPILER_NVCC
 *  KOKKOS_COMPILER_GNU
 *  KOKKOS_COMPILER_INTEL
 *  KOKKOS_COMPILER_IBM
 *  KOKKOS_COMPILER_CRAYC
 *  KOKKOS_COMPILER_APPLECC
 *  KOKKOS_COMPILER_CLANG
 *  KOKKOS_COMPILER_PGI
 *
 *  Macros for which compiler extension to use for atomics on intrinsice types
 *
 *  KOKKOS_ENABLE_CUDA_ATOMICS
 *  KOKKOS_ENABLE_GNU_ATOMICS
 *  KOKKOS_ENABLE_INTEL_ATOMICS
 *  KOKKOS_ENABLE_OPENMP_ATOMICS
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

#if defined(KOKKOS_ENABLE_SERIAL) || defined(KOKKOS_ENABLE_THREADS) || \
    defined(KOKKOS_ENABLE_OPENMP) || defined(KOKKOS_ENABLE_QTHREADS) || \
    defined(KOKKOS_ENABLE_ROCM) || defined(KOKKOS_ENABLE_OPENMPTARGET)
  #define KOKKOS_INTERNAL_ENABLE_NON_CUDA_BACKEND
#endif

#define KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA

#if defined( KOKKOS_ENABLE_CUDA ) && defined( __CUDACC__ )
  // Compiling with a CUDA compiler.
  //
  //  Include <cuda.h> to pick up the CUDA_VERSION macro defined as:
  //    CUDA_VERSION = ( MAJOR_VERSION * 1000 ) + ( MINOR_VERSION * 10 )
  //
  //  When generating device code the __CUDA_ARCH__ macro is defined as:
  //    __CUDA_ARCH__ = ( MAJOR_CAPABILITY * 100 ) + ( MINOR_CAPABILITY * 10 )

  #include <cuda_runtime.h>
  #include <cuda.h>

  #if !defined( CUDA_VERSION )
    #error "#include <cuda.h> did not define CUDA_VERSION."
  #endif

  #if ( CUDA_VERSION < 7000 )
    // CUDA supports C++11 in device code starting with version 7.0.
    // This includes auto type and device code internal lambdas.
    #error "Cuda version 7.0 or greater required."
  #endif

  #if defined( __CUDA_ARCH__ ) && ( __CUDA_ARCH__ < 300 )
    // Compiling with CUDA compiler for device code.
    #error "Cuda device capability >= 3.0 is required."
  #endif

  #ifdef KOKKOS_ENABLE_CUDA_LAMBDA
    #if ( CUDA_VERSION < 7050 )
      // CUDA supports C++11 lambdas generated in host code to be given
      // to the device starting with version 7.5. But the release candidate (7.5.6)
      // still identifies as 7.0.
      #error "Cuda version 7.5 or greater required for host-to-device Lambda support."
    #endif

    #if ( CUDA_VERSION < 8000 ) && defined( __NVCC__ )
      #define KOKKOS_LAMBDA [=]__device__
      #if defined( KOKKOS_INTERNAL_ENABLE_NON_CUDA_BACKEND )
        #undef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
      #endif
    #else
      #define KOKKOS_LAMBDA [=]__host__ __device__

      #if defined( KOKKOS_ENABLE_CXX1Z )
        #define KOKKOS_CLASS_LAMBDA        [=,*this] __host__ __device__
      #endif
    #endif

    #if defined( __NVCC__ )
      #define KOKKOS_IMPL_NEED_FUNCTOR_WRAPPER
    #endif
  #else // !defined(KOKKOS_ENABLE_CUDA_LAMBDA)
    #undef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
  #endif // !defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#endif // #if defined( KOKKOS_ENABLE_CUDA ) && defined( __CUDACC__ )

//----------------------------------------------------------------------------
// Language info: C++, CUDA, OPENMP

#if defined( KOKKOS_ENABLE_CUDA )
  // Compiling Cuda code to 'ptx'

  #define KOKKOS_FORCEINLINE_FUNCTION  __device__  __host__  __forceinline__
  #define KOKKOS_INLINE_FUNCTION       __device__  __host__  inline
  #define KOKKOS_FUNCTION              __device__  __host__
  #ifdef KOKKOS_COMPILER_CLANG
  #define KOKKOS_FUNCTION_DEFAULTED KOKKOS_FUNCTION
  #endif
#endif // #if defined( __CUDA_ARCH__ )

#if defined( KOKKOS_ENABLE_ROCM ) && defined( __HCC__ )

  #define KOKKOS_FORCEINLINE_FUNCTION  __attribute__((amp,cpu)) inline
  #define KOKKOS_INLINE_FUNCTION       __attribute__((amp,cpu)) inline
  #define KOKKOS_FUNCTION              __attribute__((amp,cpu))
  #define KOKKOS_LAMBDA                [=] __attribute__((amp,cpu))
  #define KOKKOS_FUNCTION_DEFAULTED    KOKKOS_FUNCTION
#endif

#if defined( _OPENMP )
  //  Compiling with OpenMP.
  //  The value of _OPENMP is an integer value YYYYMM
  //  where YYYY and MM are the year and month designation
  //  of the supported OpenMP API version.
#endif // #if defined( _OPENMP )

//----------------------------------------------------------------------------
// Mapping compiler built-ins to KOKKOS_COMPILER_*** macros

#if defined( __NVCC__ )
  // NVIDIA compiler is being used.
  // Code is parsed and separated into host and device code.
  // Host code is compiled again with another compiler.
  // Device code is compile to 'ptx'.
  #define KOKKOS_COMPILER_NVCC __NVCC__
#endif // #if defined( __NVCC__ )

#if !defined( KOKKOS_LAMBDA )
  #define KOKKOS_LAMBDA [=]
#endif

#if defined( KOKKOS_ENABLE_CXX1Z ) && !defined( KOKKOS_CLASS_LAMBDA )
  #define KOKKOS_CLASS_LAMBDA [=,*this]
#endif

//#if !defined( __CUDA_ARCH__ ) // Not compiling Cuda code to 'ptx'.

// Intel compiler for host code.

#if defined( __INTEL_COMPILER )
  #define KOKKOS_COMPILER_INTEL __INTEL_COMPILER
#elif defined( __ICC )
  // Old define
  #define KOKKOS_COMPILER_INTEL __ICC
#elif defined( __ECC )
  // Very old define
  #define KOKKOS_COMPILER_INTEL __ECC
#endif

// CRAY compiler for host code
#if defined( _CRAYC )
  #define KOKKOS_COMPILER_CRAYC _CRAYC
#endif

#if defined( __IBMCPP__ )
  // IBM C++
  #define KOKKOS_COMPILER_IBM __IBMCPP__
#elif defined( __IBMC__ )
  #define KOKKOS_COMPILER_IBM __IBMC__
#endif

#if defined( __APPLE_CC__ )
  #define KOKKOS_COMPILER_APPLECC __APPLE_CC__
#endif

#if defined( __clang__ ) && !defined( KOKKOS_COMPILER_INTEL )
  #define KOKKOS_COMPILER_CLANG __clang_major__*100+__clang_minor__*10+__clang_patchlevel__
#endif

#if !defined( __clang__ ) && !defined( KOKKOS_COMPILER_INTEL ) &&defined( __GNUC__ )
  #define KOKKOS_COMPILER_GNU __GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__

  #if ( 472 > KOKKOS_COMPILER_GNU )
    #error "Compiling with GCC version earlier than 4.7.2 is not supported."
  #endif
#endif

#if defined( __PGIC__ ) && !defined( __GNUC__ )
  #define KOKKOS_COMPILER_PGI __PGIC__*100+__PGIC_MINOR__*10+__PGIC_PATCHLEVEL__

  #if ( 1540 > KOKKOS_COMPILER_PGI )
    #error "Compiling with PGI version earlier than 15.4 is not supported."
  #endif
#endif

//#endif // #if !defined( __CUDA_ARCH__ )

//----------------------------------------------------------------------------
// Intel compiler macros

#if defined( KOKKOS_COMPILER_INTEL )
  #define KOKKOS_ENABLE_PRAGMA_UNROLL 1
  #define KOKKOS_ENABLE_PRAGMA_LOOPCOUNT 1
  #define KOKKOS_ENABLE_PRAGMA_VECTOR 1
  #define KOKKOS_ENABLE_PRAGMA_SIMD 1

  #if ( __INTEL_COMPILER > 1400 )
    #define KOKKOS_ENABLE_PRAGMA_IVDEP 1
  #endif

  #if ! defined( KOKKOS_MEMORY_ALIGNMENT )
    #define KOKKOS_MEMORY_ALIGNMENT 64
  #endif

  #define KOKKOS_RESTRICT __restrict__

  #ifndef KOKKOS_IMPL_ALIGN_PTR
    #define KOKKOS_IMPL_ALIGN_PTR(size) __attribute__((align_value(size)))
  #endif

  #if ( 1400 > KOKKOS_COMPILER_INTEL )
    #if ( 1300 > KOKKOS_COMPILER_INTEL )
      #error "Compiling with Intel version earlier than 13.0 is not supported. Official minimal version is 14.0."
    #else
      #warning "Compiling with Intel version 13.x probably works but is not officially supported. Official minimal version is 14.0."
    #endif
  #endif

  #if !defined( KOKKOS_ENABLE_ASM ) && !defined( _WIN32 )
    #define KOKKOS_ENABLE_ASM 1
  #endif

  #if !defined( KOKKOS_FORCEINLINE_FUNCTION )
    #if !defined( _WIN32 )
      #define KOKKOS_FORCEINLINE_FUNCTION  inline __attribute__((always_inline))
    #else
      #define KOKKOS_FORCEINLINE_FUNCTION inline
    #endif
  #endif

  #if defined( KOKKOS_ARCH_AVX512MIC )
      #define KOKKOS_ENABLE_RFO_PREFETCH 1
  #endif 

  #if defined( __MIC__ )
    // Compiling for Xeon Phi
  #endif
#endif

//----------------------------------------------------------------------------
// Cray compiler macros

#if defined( KOKKOS_COMPILER_CRAYC )
#endif

//----------------------------------------------------------------------------
// IBM Compiler macros

#if defined( KOKKOS_COMPILER_IBM )
  #define KOKKOS_ENABLE_PRAGMA_UNROLL 1
  //#define KOKKOS_ENABLE_PRAGMA_IVDEP 1
  //#define KOKKOS_ENABLE_PRAGMA_LOOPCOUNT 1
  //#define KOKKOS_ENABLE_PRAGMA_VECTOR 1
  //#define KOKKOS_ENABLE_PRAGMA_SIMD 1

  #if ! defined( KOKKOS_ENABLE_ASM )
    #define KOKKOS_ENABLE_ASM 1
  #endif
#endif

//----------------------------------------------------------------------------
// CLANG compiler macros

#if defined( KOKKOS_COMPILER_CLANG )
  //#define KOKKOS_ENABLE_PRAGMA_UNROLL 1
  //#define KOKKOS_ENABLE_PRAGMA_IVDEP 1
  //#define KOKKOS_ENABLE_PRAGMA_LOOPCOUNT 1
  //#define KOKKOS_ENABLE_PRAGMA_VECTOR 1
  //#define KOKKOS_ENABLE_PRAGMA_SIMD 1

  #if !defined( KOKKOS_FORCEINLINE_FUNCTION )
    #define KOKKOS_FORCEINLINE_FUNCTION  inline __attribute__((always_inline))
  #endif

  #if !defined( KOKKOS_IMPL_ALIGN_PTR )
    #define KOKKOS_IMPL_ALIGN_PTR(size) __attribute__((aligned(size)))
  #endif

#endif

//----------------------------------------------------------------------------
// GNU Compiler macros

#if defined( KOKKOS_COMPILER_GNU )
  //#define KOKKOS_ENABLE_PRAGMA_UNROLL 1
  //#define KOKKOS_ENABLE_PRAGMA_IVDEP 1
  //#define KOKKOS_ENABLE_PRAGMA_LOOPCOUNT 1
  //#define KOKKOS_ENABLE_PRAGMA_VECTOR 1
  //#define KOKKOS_ENABLE_PRAGMA_SIMD 1

  #if defined( KOKKOS_ARCH_AVX512MIC )
      #define KOKKOS_ENABLE_RFO_PREFETCH 1
  #endif

  #if !defined( KOKKOS_FORCEINLINE_FUNCTION )
    #define KOKKOS_FORCEINLINE_FUNCTION inline __attribute__((always_inline))
  #endif

  #if !defined( KOKKOS_ENABLE_ASM ) && !defined( __PGIC__ ) && \
      ( defined( __amd64 ) || defined( __amd64__ ) || \
        defined( __x86_64 ) || defined( __x86_64__ ) || \
	defined(__PPC64__) )
    #define KOKKOS_ENABLE_ASM 1
  #endif
#endif

//----------------------------------------------------------------------------

#if defined( KOKKOS_COMPILER_PGI )
  #define KOKKOS_ENABLE_PRAGMA_UNROLL 1
  #define KOKKOS_ENABLE_PRAGMA_IVDEP 1
  //#define KOKKOS_ENABLE_PRAGMA_LOOPCOUNT 1
  #define KOKKOS_ENABLE_PRAGMA_VECTOR 1
  //#define KOKKOS_ENABLE_PRAGMA_SIMD 1
#endif

//----------------------------------------------------------------------------

#if defined( KOKKOS_COMPILER_NVCC )
  #if defined( __CUDA_ARCH__ )
    #define KOKKOS_ENABLE_PRAGMA_UNROLL 1
  #endif
#endif

//----------------------------------------------------------------------------
// Define function marking macros if compiler specific macros are undefined:

#if !defined( KOKKOS_FORCEINLINE_FUNCTION )
  #define KOKKOS_FORCEINLINE_FUNCTION  inline
#endif

#if !defined( KOKKOS_INLINE_FUNCTION )
  #define KOKKOS_INLINE_FUNCTION  inline
#endif

#if !defined( KOKKOS_FUNCTION )
  #define KOKKOS_FUNCTION /**/
#endif

#if !defined( KOKKOS_FUNCTION_DEFAULTED )
  #define KOKKOS_FUNCTION_DEFAULTED /**/
#endif

//----------------------------------------------------------------------------
// Define empty macro for restrict if necessary:

#if !defined( KOKKOS_RESTRICT )
  #define KOKKOS_RESTRICT
#endif

//----------------------------------------------------------------------------
// Define Macro for alignment:

#if ! defined( KOKKOS_MEMORY_ALIGNMENT )
  #define KOKKOS_MEMORY_ALIGNMENT 16
#endif

#if ! defined( KOKKOS_MEMORY_ALIGNMENT_THRESHOLD )
  #define KOKKOS_MEMORY_ALIGNMENT_THRESHOLD 4
#endif

#if !defined( KOKKOS_IMPL_ALIGN_PTR )
  #define KOKKOS_IMPL_ALIGN_PTR(size) /* */
#endif

//----------------------------------------------------------------------------
// Determine the default execution space for parallel dispatch.
// There is zero or one default execution space specified.

#if 1 < ( ( defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA ) ? 1 : 0 ) + \
          ( defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_ROCM ) ? 1 : 0 ) + \
          ( defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMPTARGET ) ? 1 : 0 ) + \
          ( defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP ) ? 1 : 0 ) + \
          ( defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS ) ? 1 : 0 ) + \
          ( defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_QTHREADS ) ? 1 : 0 ) + \
          ( defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL ) ? 1 : 0 ) )
  #error "More than one KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_* specified."
#endif

// If default is not specified then chose from enabled execution spaces.
// Priority: CUDA, OPENMP, THREADS, QTHREADS, SERIAL
#if   defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA )
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_ROCM )
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMPTARGET )
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP )
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS )
//#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_QTHREADS )
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL )
#elif defined( KOKKOS_ENABLE_CUDA )
  #define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_CUDA
#elif defined( KOKKOS_ENABLE_ROCM )
  #define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_ROCM
#elif defined( KOKKOS_ENABLE_OPENMPTARGET )
  #define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMPTARGET
#elif defined( KOKKOS_ENABLE_OPENMP )
  #define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP
#elif defined( KOKKOS_ENABLE_THREADS )
  #define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS
//#elif defined( KOKKOS_ENABLE_QTHREADS )
//  #define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_QTHREADS
#else
  #define KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_SERIAL
#endif

//----------------------------------------------------------------------------
// Determine for what space the code is being compiled:

#if defined( __CUDACC__ ) && defined( __CUDA_ARCH__ ) && defined( KOKKOS_ENABLE_CUDA )
  #define KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA
#elif   defined( __HCC__ ) && defined( __HCC_ACCELERATOR__ ) && defined( KOKKOS_ENABLE_ROCM )
  #define KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_ROCM_GPU
#else
  #define KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
#endif

//----------------------------------------------------------------------------

#if ( defined( _POSIX_C_SOURCE ) && _POSIX_C_SOURCE >= 200112L ) || \
    ( defined( _XOPEN_SOURCE )   && _XOPEN_SOURCE   >= 600 )
  #if defined( KOKKOS_ENABLE_PERFORMANCE_POSIX_MEMALIGN )
    #define KOKKOS_ENABLE_POSIX_MEMALIGN 1
  #endif
#endif

//----------------------------------------------------------------------------
// If compiling with CUDA then must be using CUDA 8 or better
// and use relocateable device code to enable the task policy.
// nvcc relocatable device code option: --relocatable-device-code=true

#if ( defined( KOKKOS_ENABLE_CUDA ) )
  #if ( 8000 <= CUDA_VERSION ) && defined( KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE )
  #define KOKKOS_ENABLE_TASKDAG
  #endif
#else
  #define KOKKOS_ENABLE_TASKDAG
#endif

#endif // #ifndef KOKKOS_MACROS_HPP

