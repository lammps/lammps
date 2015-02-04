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

#ifndef KOKKOS_HOSTSPACE_HPP
#define KOKKOS_HOSTSPACE_HPP

#include <iosfwd>
#include <typeinfo>
#include <string>

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/// \class HostSpace
/// \brief Memory management for host memory.
///
/// HostSpace is a memory space that governs host memory.  "Host"
/// memory means the usual CPU-accessible memory.
class HostSpace {
public:

  //! Tag this class as a kokkos memory space
  typedef HostSpace  memory_space ;
  typedef size_t     size_type ;

  /// \typedef execution_space
  /// \brief Default execution space for this memory space.
  ///
  /// Every memory space has a default execution space.  This is
  /// useful for things like initializing a View (which happens in
  /// parallel using the View's default execution space).
#if defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_OPENMP )
  typedef Kokkos::OpenMP   execution_space ;
#elif defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS )
  typedef Kokkos::Threads  execution_space ;
#elif defined( KOKKOS_HAVE_OPENMP )
  typedef Kokkos::OpenMP   execution_space ;
#elif defined( KOKKOS_HAVE_PTHREAD )
  typedef Kokkos::Threads  execution_space ;
#elif defined( KOKKOS_HAVE_SERIAL )
  typedef Kokkos::Serial   execution_space ;
#else
#  error "At least one of the following host execution spaces must be defined: Kokkos::OpenMP, Kokkos::Serial, or Kokkos::Threads.  You might be seeing this message if you disabled the Kokkos::Serial device explicitly using the Kokkos_ENABLE_Serial:BOOL=OFF CMake option, but did not enable any of the other host execution space devices."
#endif

  /** \brief  Allocate a contiguous block of memory.
   *
   *  The input label is associated with the block of memory.
   *  The block of memory is tracked via reference counting where
   *  allocation gives it a reference count of one.
   *
   *  Allocation may only occur on the master thread of the process.
   */
  static void * allocate( const std::string & label , const size_t size );

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

  /** \brief  Get the reference count of the block of memory
   *          in which the input pointer resides.  If the reference
   *          count is zero the memory region is not tracked.
   *
   *          Reference counting only occurs on the master thread.
   */
  static int count( const void * );

  /*--------------------------------*/

  /** \brief  Print all tracked memory to the output stream. */
  static void print_memory_view( std::ostream & );

  /** \brief  Retrieve label associated with the input pointer */
  static std::string query_label( const void * );

  /*--------------------------------*/
  /* Functions unique to the HostSpace */

  static int in_parallel();

  static void register_in_parallel( int (*)() );
};


} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class , class > struct DeepCopy ;

template<>
struct DeepCopy<HostSpace,HostSpace> {
  DeepCopy( void * dst , const void * src , size_t n );
};

} // namespace Impl
} // namespace Kokkos

#endif /* #define KOKKOS_HOSTSPACE_HPP */

