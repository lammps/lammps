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

#ifndef KOKKOS_MACROS_HPP
#define KOKKOS_MACROS_HPP

//----------------------------------------------------------------------------
/** Pick up configure/build options via #define macros:
 *
 *  KOKKOS_HAVE_CUDA                Kokkos::Cuda execution and memory spaces
 *  KOKKOS_HAVE_PTHREAD             Kokkos::Threads execution space
 *  KOKKOS_HAVE_QTHREAD             Kokkos::Qthread execution space
 *  KOKKOS_HAVE_OPENMP              Kokkos::OpenMP  execution space
 *  KOKKOS_HAVE_HWLOC               HWLOC library is available
 *  KOKKOS_HAVE_EXPRESSION_CHECK    insert array bounds checks, is expensive!
 *  KOKKOS_HAVE_CXX11               enable C++11 features
 *
 *  KOKKOS_HAVE_MPI                 negotiate MPI/execution space interactions
 *
 *  KOKKOS_USE_CUDA_UVM             Use CUDA UVM for Cuda memory space
 */

#ifndef KOKKOS_DONT_INCLUDE_CORE_CONFIG_H
#include <KokkosCore_config.h>
#endif

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
 *  KOKKOS_ATOMICS_USE_CUDA
 *  KOKKOS_ATOMICS_USE_GNU
 *  KOKKOS_ATOMICS_USE_INTEL
 *  KOKKOS_ATOMICS_USE_OPENMP31
 *
 *  A suite of 'KOKKOS_HAVE_PRAGMA_...' are defined for internal use.
 *
 *  Macros for marking functions to run in an execution space:
 *
 *  KOKKOS_FUNCTION
 *  KOKKOS_INLINE_FUNCTION        request compiler to inline
 *  KOKKOS_FORCEINLINE_FUNCTION   force compiler to inline, use with care!
 */

#include <impl/Kokkos_Compiler_Macros.hpp>

/** Define function marking macros if compiler specific macros are undefined: */

#if ! defined( KOKKOS_FORCEINLINE_FUNCTION )
#define KOKKOS_FORCEINLINE_FUNCTION  inline
#endif

#if ! defined( KOKKOS_INLINE_FUNCTION )
#define KOKKOS_INLINE_FUNCTION  inline
#endif

#if ! defined( KOKKOS_FUNCTION )
#define KOKKOS_FUNCTION /**/
#endif

/** These should be part of the Atomics API */

#if ! defined( KOKKOS_NONTEMPORAL_PREFETCH_LOAD )
    #define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) ((void)0)
    #define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) ((void)0)
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Non-macro forward declaration placement in this file to be reconsidered...
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Forward declarations for enabled execution and memory spaces.

namespace Kokkos {

class Serial ;    ///< Execution space for serial on CPU
class HostSpace ; ///< Memory space for Serial, Threads, and OpenMP

class Threads ; ///< Pthreads execution space

#if defined( KOKKOS_HAVE_CUDA )
class CudaSpace ; ///< Cuda memory space
class Cuda ;      ///< Cuda execution space
#endif

#if defined( KOKKOS_HAVE_OPENMP )
class OpenMP ; ///< OpenMP execution space
#endif

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Set the default execution space.

/// Define Kokkos::DefaultExecutionSpace as per configuration option
/// or chosen from the enabled execution spaces in the following order:
/// Kokkos::Cuda, Kokkos::OpenMP, Kokkos::Threads, Kokkos::Serial

namespace Kokkos {

#if 1 < ( ( defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA ) ? 1 : 0 ) + \
          ( defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP ) ? 1 : 0 ) + \
          ( defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS ) ? 1 : 0 ) + \
          ( defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL ) ? 1 : 0 ) )

#error "More than one KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_* specified" ;

#endif

#if   defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA )
  typedef Kokkos::Cuda DefaultExecutionSpace ;
#elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP )
  typedef OpenMP DefaultExecutionSpace ;
#elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS )
  typedef Threads DefaultExecutionSpace ;
#elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
  typedef Serial DefaultExecutionSpace ;
#elif  defined ( KOKKOS_HAVE_CUDA )
  #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA
  typedef Kokkos::Cuda DefaultExecutionSpace ;
#elif defined ( KOKKOS_HAVE_OPENMP )
  #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP
  typedef Kokkos::OpenMP DefaultExecutionSpace ;
#elif defined ( KOKKOS_HAVE_PTHREAD )
  #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS
  typedef Kokkos::Threads DefaultExecutionSpace ;
#else
  #define KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL
  typedef Kokkos::Serial DefaultExecutionSpace ;
#endif

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#if defined( __CUDACC__ ) && defined( __CUDA_ARCH__ )
typedef Kokkos::CudaSpace  ActiveExecutionMemorySpace ;
#define KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA
#else
typedef Kokkos::HostSpace  ActiveExecutionMemorySpace ;
#define KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
#endif

template< class ActiveSpace , class MemorySpace >
struct VerifyExecutionCanAccessMemorySpace {};

template< class Space >
struct VerifyExecutionCanAccessMemorySpace< Space , Space >
{
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void *) {}
};

} // namespace Impl
} // namespace Kokkos

// Currently executing in the CUDA space

#define KOKKOS_RESTRICT_EXECUTION_TO_DATA( DATA_SPACE , DATA_PTR ) \
  Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< \
    Kokkos::Impl::ActiveExecutionMemorySpace , DATA_SPACE >::verify( DATA_PTR )

#define KOKKOS_RESTRICT_EXECUTION_TO_( DATA_SPACE ) \
  Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< \
    Kokkos::Impl::ActiveExecutionMemorySpace , DATA_SPACE >::verify()

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_MACROS_HPP */

