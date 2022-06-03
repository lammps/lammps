/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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

#ifdef KOKKOS_ENABLE_CUDA_LAMBDA
#define KOKKOS_LAMBDA [=] __host__ __device__

#if defined(KOKKOS_ENABLE_CXX17) || defined(KOKKOS_ENABLE_CXX20)
#define KOKKOS_CLASS_LAMBDA [ =, *this ] __host__ __device__
#endif

#else  // !defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#undef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
#endif  // !defined(KOKKOS_ENABLE_CUDA_LAMBDA)

#if (10000 > CUDA_VERSION)
#define KOKKOS_ENABLE_PRE_CUDA_10_DEPRECATION_API
#endif

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 700)
// PTX atomics with memory order semantics are only available on volta and later
#if !defined(KOKKOS_DISABLE_CUDA_ASM)
#if !defined(KOKKOS_ENABLE_CUDA_ASM)
#define KOKKOS_ENABLE_CUDA_ASM
#if !defined(KOKKOS_DISABLE_CUDA_ASM_ATOMICS) && \
    defined(KOKKOS_ENABLE_GNU_ATOMICS)
#define KOKKOS_ENABLE_CUDA_ASM_ATOMICS
#endif
#endif
#endif
#endif

#define KOKKOS_IMPL_FORCEINLINE_FUNCTION __device__ __host__ __forceinline__
#define KOKKOS_IMPL_FORCEINLINE __forceinline__
#define KOKKOS_IMPL_INLINE_FUNCTION __device__ __host__ inline
#define KOKKOS_IMPL_FUNCTION __device__ __host__
#define KOKKOS_IMPL_HOST_FUNCTION __host__
#define KOKKOS_IMPL_DEVICE_FUNCTION __device__
#if defined(KOKKOS_COMPILER_NVCC)
#define KOKKOS_INLINE_FUNCTION_DELETED inline
#else
#define KOKKOS_INLINE_FUNCTION_DELETED __device__ __host__ inline
#endif
#if (CUDA_VERSION < 10000)
#define KOKKOS_DEFAULTED_FUNCTION __host__ __device__ inline
#else
#define KOKKOS_DEFAULTED_FUNCTION inline
#endif
#define KOKKOS_IMPL_HOST_FUNCTION __host__
#define KOKKOS_IMPL_DEVICE_FUNCTION __device__

#if (CUDA_VERSION >= 10000)
#define KOKKOS_CUDA_ENABLE_GRAPHS
#endif

#endif /* KOKKOS_CUDA_SETUP_HPP_ */
