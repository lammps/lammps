/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_CUDASPACE_HPP
#define KOKKOS_CUDASPACE_HPP

#if defined( __CUDACC__ )
#include <cuda_runtime.h>
#endif

#include <iosfwd>
#include <typeinfo>
#include <string>

#include <Kokkos_Macros.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Cuda/Kokkos_Cuda_abort.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Cuda memory management */

class CudaSpace {
public:

  typedef CudaSpace     memory_space ;
  typedef unsigned int  size_type ;

  /** \brief  Allocate a contiguous block of memory on the Cuda device
   *          with size = scalar_size * scalar_count.
   *
   *  The input label is associated with the block of memory.
   *  The block of memory is tracked via reference counting where
   *  allocation gives it a reference count of one.
   *
   *  Allocation may only occur on the master thread of the process.
   */
  static void * allocate( const std::string    & label ,
                          const std::type_info & scalar_type ,
                          const size_t           scalar_size ,
                          const size_t           scalar_count );

  /** \brief  Increment the reference count of the block of memory
   *          in which the input pointer resides.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void increment( const void * );

  /** \brief  Decrement the reference count of the block of memory
   *          in which the input pointer resides.  If the reference
   *          count falls to zero the memory is deallocated.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void decrement( const void * );

  /** \brief  Print all tracked memory to the output stream. */
  static void print_memory_view( std::ostream & );

  /** \brief  Retrieve label associated with the input pointer */
  static std::string query_label( const void * );

  /*--------------------------------*/

  static void access_error();
  static void access_error( const void * const );

  /*--------------------------------*/
};

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
struct DeepCopy<HostSpace,CudaSpace> {
  DeepCopy( void * dst , const void * src , size_t );
};

template<>
struct DeepCopy<CudaSpace,HostSpace> {
  DeepCopy( void * dst , const void * src , size_t );
};

template<>
struct DeepCopy<CudaSpace,CudaSpace> {
  DeepCopy( void * dst , const void * src , size_t );
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/** \brief  Cuda code accessing Cuda data is good. */
template<>
struct VerifyExecutionSpaceCanAccessDataSpace< CudaSpace , CudaSpace >
{
  KOKKOS_INLINE_FUNCTION static void verify( void ) {}
  KOKKOS_INLINE_FUNCTION static void verify( const void * ) {}
};

/** \brief  Cuda code accessing non-Cuda data is bad. */
template<>
struct VerifyExecutionSpaceCanAccessDataSpace< CudaSpace , HostSpace >
{
  KOKKOS_INLINE_FUNCTION static void verify(void)
  { Kokkos::cuda_abort("Cuda code called function restricted to HostSpace"); }

  KOKKOS_INLINE_FUNCTION static void verify( const void * )
  { Kokkos::cuda_abort("Cuda code attempted to access HostSpace memory"); }
};

/** \brief  Produce error message when trying to access Cuda 
 *          memory on the host.
 */
template<>
struct VerifyExecutionSpaceCanAccessDataSpace< HostSpace , CudaSpace >
{
#ifdef KOKKOS_USE_UVM
  inline static void verify( void ) { }
  inline static void verify( const void * p ) { }
#else
  inline static void verify( void ) { CudaSpace::access_error(); }
  inline static void verify( const void * p ) { CudaSpace::access_error(p); }
#endif
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #define KOKKOS_CUDASPACE_HPP */

