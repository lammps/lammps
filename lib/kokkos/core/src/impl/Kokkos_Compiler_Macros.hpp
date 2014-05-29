/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_COMPILER_MACROS_HPP
#define KOKKOS_COMPILER_MACROS_HPP

#if defined __ECC || defined __ICC || defined __INTEL_COMPILER
  #define KOKKOS_COMPILER_NAME "Intel C++"
  #if defined __ICC
    #define KOKKOS_COMPILER_VERSION __ICC
  #else
    #if defined __INTEL_COMPILER
      #define KOKKOS_COMPILER_VERSION __INTEL_COMPILER
    #else
      #define KOKKOS_COMPILER_VERSION __ECC
    #endif
  #endif

  #define KOKKOS_COMPILER_INTEL 1

  #define KOKKOS_HAVE_PRAGMA_UNROLL 1
  #define KOKKOS_HAVE_PRAGMA_IVDEP 1
  #define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  #define KOKKOS_HAVE_PRAGMA_VECTOR 1
  #define KOKKOS_HAVE_PRAGMA_SIMD 1
#endif

#if defined __IBMC__ || defined __IBMCPP__
  #define KOKKOS_COMPILER_NAME "IBM C++"
  #if defined __IBMC__
    #define KOKKOS_COMPILER_VERSION __IBMC__
  #else
    #define KOKKOS_COMPILER_VERSION __IBMCPP__
  #endif
  #define KOKKOS_COMPILER_IBM 1

  #define KOKKOS_HAVE_PRAGMA_UNROLL 1
  //#define KOKKOS_HAVE_PRAGMA_IVDEP 1
  //#define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  //#define KOKKOS_HAVE_PRAGMA_VECTOR 1
  //#define KOKKOS_HAVE_PRAGMA_SIMD 1
#endif

#if defined __APPLE_CC__
   /* Apple uses GNU as compiler */
  #define KOKKOS_COMPILER_APPLECC 1
#endif

#if defined __clang__
  #define KOKKOS_COMPILER_NAME "Clang"
  #define KOKKOS_COMPILER_VERSION __clang_major__*100+__clang_minor__*10+__clang_patchlevel__
  #define KOKKOS_COMPILER_CLANG 1

  //#define KOKKOS_HAVE_PRAGMA_UNROLL 1
  //#define KOKKOS_HAVE_PRAGMA_IVDEP 1
  //#define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  //#define KOKKOS_HAVE_PRAGMA_VECTOR 1
  //#define KOKKOS_HAVE_PRAGMA_SIMD 1
#endif

#if defined __GNUC__ && !defined KOKKOS_COMPILER_NAME && !defined __clang__
  #define KOKKOS_COMPILER_NAME "Gnu GCC"
  #define KOKKOS_COMPILER_VERSION __GNUC__*100+__GNUC_MINOR__*10+__GNUC_PATCHLEVEL__
  #define KOKKOS_COMPILER_GCC 1

  //#define KOKKOS_HAVE_PRAGMA_UNROLL 1
  //#define KOKKOS_HAVE_PRAGMA_IVDEP 1
  //#define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  //#define KOKKOS_HAVE_PRAGMA_VECTOR 1
  //#define KOKKOS_HAVE_PRAGMA_SIMD 1
#endif

#if defined __PGIC__ && !defined KOKKOS_COMPILER_NAME
  #define KOKKOS_COMPILER_NAME "PGI C++"
  #define KOKKOS_COMPILER_VERSION __PGIC__*100+__PGIC_MINOR__*10+__PGIC_PATCHLEVEL__
  #define KOKKOS_COMPILER_PGI 1

  #define KOKKOS_HAVE_PRAGMA_UNROLL 1
  #define KOKKOS_HAVE_PRAGMA_IVDEP 1
  //#define KOKKOS_HAVE_PRAGMA_LOOPCOUNT 1
  #define KOKKOS_HAVE_PRAGMA_VECTOR 1
  //#define KOKKOS_HAVE_PRAGMA_SIMD 1
#endif

#if defined __NVCC__
  #define KOKKOS_DEVICE_COMPILER_NAME "NVIDIA NVCC"
  #define KOKKOS_DEVICE_COMPILER_VERSION __NVCC__
  #if (!defined(KOKKOS_HAVE_PRAGMA_UNROLL) && defined(__CUDA_ARCH__))
  #define KOKKOS_HAVE_PRAGMA_UNROLL 1
  #endif
#endif

#if !defined KOKKOS_COMPILER_NAME
  #define KOKKOS_COMPILER_NAME "Unknown compiler"
#endif

#if !defined KOKKOS_COMPILER_VERSION
  #define KOKKOS_COMPILER_VERSION 0
#endif

#if !defined KOKKOS_DEVICE_COMPILER_NAME
  #define KOKKOS_DEVICE_COMPILER_NAME KOKKOS_COMPILER_NAME
#endif

#if !defined KOKKOS_DEVICE_COMPILER_VERSION
  #define KOKKOS_DEVICE_COMPILER_VERSION KOKKOS_COMPILER_VERSION
#endif


#endif
