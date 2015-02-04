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

#ifndef KOKKOS_CORE_FWD_HPP
#define KOKKOS_CORE_FWD_HPP

//----------------------------------------------------------------------------
// Kokkos_Macros.hpp does introspection on configuration options
// and compiler environment then sets a collection of #define macros.

#include <Kokkos_Macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Forward declarations for class inter-relationships

namespace Kokkos {

class HostSpace ; ///< Memory space for main process and CPU execution spaces

#if defined( KOKKOS_HAVE_SERIAL )
class Serial ;    ///< Execution space main process on CPU
#endif // defined( KOKKOS_HAVE_SERIAL )

#if defined( KOKKOS_HAVE_PTHREAD )
class Threads ;  ///< Execution space with pthreads back-end
#endif

#if defined( KOKKOS_HAVE_OPENMP )
class OpenMP ; ///< OpenMP execution space
#endif

#if defined( KOKKOS_HAVE_CUDA )
class CudaSpace ;            ///< Memory space on Cuda GPU
class CudaUVMSpace ;         ///< Memory space on Cuda GPU with UVM
class CudaHostPinnedSpace ;  ///< Memory space on Host accessible to Cuda GPU
class Cuda ;                 ///< Execution space for Cuda GPU
#endif

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Set the default execution space.

/// Define Kokkos::DefaultExecutionSpace as per configuration option
/// or chosen from the enabled execution spaces in the following order:
/// Kokkos::Cuda, Kokkos::OpenMP, Kokkos::Threads, Kokkos::Serial

namespace Kokkos {

#if   defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_CUDA )
  typedef Kokkos::Cuda DefaultExecutionSpace ;
#elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP )
  typedef OpenMP DefaultExecutionSpace ;
#elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS )
  typedef Threads DefaultExecutionSpace ;
#elif defined ( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
  typedef Serial DefaultExecutionSpace ;
#else
#  error "At least one of the following execution spaces must be defined in order to use Kokkos: Kokkos::Cuda, Kokkos::OpenMP, Kokkos::Serial, or Kokkos::Threads."
#endif

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Detect the active execution space and define its memory space.
// This is used to verify whether a running kernel can access
// a given memory space.

namespace Kokkos {
namespace Impl {

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA ) && defined (KOKKOS_HAVE_CUDA)
typedef Kokkos::CudaSpace  ActiveExecutionMemorySpace ;
#elif defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
typedef Kokkos::HostSpace  ActiveExecutionMemorySpace ;
#else
typedef void ActiveExecutionMemorySpace ;
#endif

template< class ActiveSpace , class MemorySpace >
struct VerifyExecutionCanAccessMemorySpace {
  enum {value = 0};
};

template< class Space >
struct VerifyExecutionCanAccessMemorySpace< Space , Space >
{
  enum {value = 1};
  KOKKOS_INLINE_FUNCTION static void verify(void) {}
  KOKKOS_INLINE_FUNCTION static void verify(const void *) {}
};

} // namespace Impl
} // namespace Kokkos

#define KOKKOS_RESTRICT_EXECUTION_TO_DATA( DATA_SPACE , DATA_PTR ) \
  Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< \
    Kokkos::Impl::ActiveExecutionMemorySpace , DATA_SPACE >::verify( DATA_PTR )

#define KOKKOS_RESTRICT_EXECUTION_TO_( DATA_SPACE ) \
  Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< \
    Kokkos::Impl::ActiveExecutionMemorySpace , DATA_SPACE >::verify()

#endif /* #ifndef KOKKOS_CORE_FWD_HPP */

