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

#if defined( KOKKOS_HAVE_CUDA )

#include <iosfwd>
#include <typeinfo>
#include <string>

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Cuda/Kokkos_Cuda_abort.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Cuda on-device memory management */

class CudaSpace {
public:

  //! Tag this class as a kokkos memory space
  typedef CudaSpace             memory_space ;
  typedef Kokkos::Cuda          execution_space ;
  typedef unsigned int          size_type ;

  /** \brief  Allocate a contiguous block of memory on the Cuda device.
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

  /** \brief  Print all tracked memory to the output stream. */
  static void print_memory_view( std::ostream & );

  /** \brief  Retrieve label associated with the input pointer */
  static std::string query_label( const void * );

  /*--------------------------------*/
  /** \brief  Cuda specific function to attached texture object to an allocation.
   *          Output the texture object, base pointer, and offset from the input pointer.
   */
#if defined( __CUDACC__ )
  static void texture_object_attach( const void            * const arg_ptr
                                   , const unsigned                arg_type_size
                                   , ::cudaChannelFormatDesc const & arg_desc
                                   , ::cudaTextureObject_t * const arg_tex_obj
                                   , void const           ** const arg_alloc_ptr
                                   , int                   * const arg_offset
                                   );
#endif

  /*--------------------------------*/
  /** \brief  Error reporting for HostSpace attempt to access CudaSpace */
  static void access_error();
  static void access_error( const void * const );
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Cuda memory that is accessible to Host execution space
 *          through Cuda's unified virtual memory (UVM) runtime.
 */
class CudaUVMSpace {
public:

  //! Tag this class as a kokkos memory space
  typedef CudaUVMSpace          memory_space ;
  typedef Cuda                  execution_space ;
  typedef unsigned int          size_type ;

  /** \brief  If UVM capability is available */
  static bool available();

  /** \brief  Allocate a contiguous block of memory on the Cuda device.
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

  /** \brief  Print all tracked memory to the output stream. */
  static void print_memory_view( std::ostream & );

  /** \brief  Retrieve label associated with the input pointer */
  static std::string query_label( const void * );

  /** \brief  Cuda specific function to attached texture object to an allocation.
   *          Output the texture object, base pointer, and offset from the input pointer.
   */
#if defined( __CUDACC__ )
  static void texture_object_attach( const void            * const arg_ptr
                                   , const unsigned                arg_type_size
                                   , ::cudaChannelFormatDesc const & arg_desc
                                   , ::cudaTextureObject_t * const arg_tex_obj
                                   , void const           ** const arg_alloc_ptr
                                   , int                   * const arg_offset
                                   );
#endif
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Host memory that is accessible to Cuda execution space
 *          through Cuda's host-pinned memory allocation.
 */
class CudaHostPinnedSpace {
public:

  //! Tag this class as a kokkos memory space
  typedef CudaHostPinnedSpace         memory_space ;
  typedef unsigned int                size_type ;

  /** \brief  Memory is in HostSpace so use the HostSpace::execution_space */
  typedef HostSpace::execution_space  execution_space ;

  /** \brief  Allocate a contiguous block of memory on the Cuda device.
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

  /** \brief  Get the reference count of the block of memory
   *          in which the input pointer resides.  If the reference
   *          count is zero the memory region is not tracked.
   *
   *          Reference counting only occurs on the master thread.
   */
  static int count( const void * );

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
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template<> struct DeepCopy< CudaSpace , CudaSpace >
{
  DeepCopy( void * dst , const void * src , size_t );
  DeepCopy( const Cuda & , void * dst , const void * src , size_t );
};

template<> struct DeepCopy< CudaSpace , HostSpace >
{
  DeepCopy( void * dst , const void * src , size_t );
  DeepCopy( const Cuda & , void * dst , const void * src , size_t );
};

template<> struct DeepCopy< HostSpace , CudaSpace >
{
  DeepCopy( void * dst , const void * src , size_t );
  DeepCopy( const Cuda & , void * dst , const void * src , size_t );
};

template<> struct DeepCopy< CudaSpace , CudaUVMSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaSpace , CudaHostPinnedSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , HostSpace >( dst , src , n ); }
};


template<> struct DeepCopy< CudaUVMSpace , CudaSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaUVMSpace , CudaUVMSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaUVMSpace , CudaHostPinnedSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , HostSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaUVMSpace , HostSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , HostSpace >( dst , src , n ); }
};


template<> struct DeepCopy< CudaHostPinnedSpace , CudaSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaHostPinnedSpace , CudaUVMSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaHostPinnedSpace , CudaHostPinnedSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , HostSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaHostPinnedSpace , HostSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , HostSpace >( dst , src , n ); }
};


template<> struct DeepCopy< HostSpace , CudaUVMSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< HostSpace , CudaHostPinnedSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , HostSpace >( dst , src , n ); }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** Running in CudaSpace attempting to access HostSpace: error */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::CudaSpace , Kokkos::HostSpace >
{
  enum { value = false };
  KOKKOS_INLINE_FUNCTION static void verify( void )
    { Kokkos::abort("Cuda code attempted to access HostSpace memory"); }

  KOKKOS_INLINE_FUNCTION static void verify( const void * )
    { Kokkos::abort("Cuda code attempted to access HostSpace memory"); }
};

/** Running in CudaSpace accessing CudaUVMSpace: ok */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::CudaSpace , Kokkos::CudaUVMSpace >
{
  enum { value = true };
  KOKKOS_INLINE_FUNCTION static void verify( void ) { }
  KOKKOS_INLINE_FUNCTION static void verify( const void * ) { }
};

/** Running in CudaSpace accessing CudaHostPinnedSpace: ok */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::CudaSpace , Kokkos::CudaHostPinnedSpace >
{
  enum { value = true };
  KOKKOS_INLINE_FUNCTION static void verify( void ) { }
  KOKKOS_INLINE_FUNCTION static void verify( const void * ) { }
};

/** Running in CudaSpace attempting to access an unknown space: error */
template< class OtherSpace >
struct VerifyExecutionCanAccessMemorySpace<
  typename enable_if< ! is_same<Kokkos::CudaSpace,OtherSpace>::value , Kokkos::CudaSpace >::type ,
  OtherSpace >
{
  enum { value = false };
  KOKKOS_INLINE_FUNCTION static void verify( void )
    { Kokkos::abort("Cuda code attempted to access unknown Space memory"); }

  KOKKOS_INLINE_FUNCTION static void verify( const void * )
    { Kokkos::abort("Cuda code attempted to access unknown Space memory"); }
};

//----------------------------------------------------------------------------
/** Running in HostSpace attempting to access CudaSpace */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::HostSpace , Kokkos::CudaSpace >
{
  enum { value = false };
  inline static void verify( void ) { CudaSpace::access_error(); }
  inline static void verify( const void * p ) { CudaSpace::access_error(p); }
};

/** Running in HostSpace accessing CudaUVMSpace is OK */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::HostSpace , Kokkos::CudaUVMSpace >
{
  enum { value = true };
  inline static void verify( void ) { }
  inline static void verify( const void * ) { }
};

/** Running in HostSpace accessing CudaHostPinnedSpace is OK */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::HostSpace , Kokkos::CudaHostPinnedSpace >
{
  enum { value = true };
  KOKKOS_INLINE_FUNCTION static void verify( void ) {}
  KOKKOS_INLINE_FUNCTION static void verify( const void * ) {}
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_HAVE_CUDA ) */
#endif /* #define KOKKOS_CUDASPACE_HPP */

