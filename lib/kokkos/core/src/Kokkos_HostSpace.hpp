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

#ifndef KOKKOS_HOSTSPACE_HPP
#define KOKKOS_HOSTSPACE_HPP

#include <cstring>
#include <string>
#include <iosfwd>
#include <typeinfo>

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Concepts.hpp>
#include <Kokkos_MemoryTraits.hpp>

#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

namespace Impl {

/// \brief Initialize lock array for arbitrary size atomics.
///
/// Arbitrary atomics are implemented using a hash table of locks
/// where the hash value is derived from the address of the
/// object for which an atomic operation is performed.
/// This function initializes the locks to zero (unset).
void init_lock_array_host_space();

/// \brief Aquire a lock for the address
///
/// This function tries to aquire the lock for the hash value derived
/// from the provided ptr. If the lock is successfully aquired the
/// function returns true. Otherwise it returns false.
bool lock_address_host_space(void* ptr);

/// \brief Release lock for the address
///
/// This function releases the lock for the hash value derived
/// from the provided ptr. This function should only be called
/// after previously successfully aquiring a lock with
/// lock_address.
void unlock_address_host_space( void* ptr );

} // namespace Impl

} // namespace Kokkos

namespace Kokkos {

/// \class HostSpace
/// \brief Memory management for host memory.
///
/// HostSpace is a memory space that governs host memory.  "Host"
/// memory means the usual CPU-accessible memory.
class HostSpace {
public:
  //! Tag this class as a kokkos memory space
  typedef HostSpace  memory_space;
  typedef size_t     size_type;

  /// \typedef execution_space
  /// \brief Default execution space for this memory space.
  ///
  /// Every memory space has a default execution space.  This is
  /// useful for things like initializing a View (which happens in
  /// parallel using the View's default execution space).
#if defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_OPENMP )
  typedef Kokkos::OpenMP    execution_space;
#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_THREADS )
  typedef Kokkos::Threads   execution_space;
//#elif defined( KOKKOS_ENABLE_DEFAULT_DEVICE_TYPE_QTHREADS )
//  typedef Kokkos::Qthreads  execution_space;
#elif defined( KOKKOS_ENABLE_OPENMP )
  typedef Kokkos::OpenMP    execution_space;
#elif defined( KOKKOS_ENABLE_THREADS )
  typedef Kokkos::Threads   execution_space;
//#elif defined( KOKKOS_ENABLE_QTHREADS )
//  typedef Kokkos::Qthreads  execution_space;
#elif defined( KOKKOS_ENABLE_SERIAL )
  typedef Kokkos::Serial    execution_space;
#else
#  error "At least one of the following host execution spaces must be defined: Kokkos::OpenMP, Kokkos::Threads, Kokkos::Qthreads, or Kokkos::Serial.  You might be seeing this message if you disabled the Kokkos::Serial device explicitly using the Kokkos_ENABLE_Serial:BOOL=OFF CMake option, but did not enable any of the other host execution space devices."
#endif

  //! This memory space preferred device_type
  typedef Kokkos::Device< execution_space, memory_space > device_type;

  /**\brief  Default memory space instance */
  HostSpace();
  HostSpace( HostSpace && rhs ) = default;
  HostSpace( const HostSpace & rhs ) = default;
  HostSpace & operator = ( HostSpace && ) = default;
  HostSpace & operator = ( const HostSpace & ) = default;
  ~HostSpace() = default;

  /**\brief  Non-default memory space instance to choose allocation mechansim, if available */

  enum AllocationMechanism { STD_MALLOC, POSIX_MEMALIGN, POSIX_MMAP, INTEL_MM_ALLOC };

  explicit
  HostSpace( const AllocationMechanism & );

  /**\brief  Allocate untracked memory in the space */
  void * allocate( const size_t arg_alloc_size ) const;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate( void * const arg_alloc_ptr
                 , const size_t arg_alloc_size ) const;

  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name() { return m_name; }

private:
  AllocationMechanism  m_alloc_mech;
  static constexpr const char* m_name = "Host";
  friend class Kokkos::Impl::SharedAllocationRecord< Kokkos::HostSpace, void >;
};

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

static_assert( Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::HostSpace >::assignable, "" );

template< typename S >
struct HostMirror {
private:
  // If input execution space can access HostSpace then keep it.
  // Example: Kokkos::OpenMP can access, Kokkos::Cuda cannot
  enum { keep_exe = Kokkos::Impl::MemorySpaceAccess
                      < typename S::execution_space::memory_space, Kokkos::HostSpace >::accessible };

  // If HostSpace can access memory space then keep it.
  // Example:  Cannot access Kokkos::CudaSpace, can access Kokkos::CudaUVMSpace
  enum { keep_mem = Kokkos::Impl::MemorySpaceAccess
                      < Kokkos::HostSpace, typename S::memory_space >::accessible };

public:

  typedef typename std::conditional
    < keep_exe && keep_mem /* Can keep whole space */
    , S
    , typename std::conditional
        < keep_mem /* Can keep memory space, use default Host execution space */
        , Kokkos::Device< Kokkos::HostSpace::execution_space
                        , typename S::memory_space >
        , Kokkos::HostSpace
        >::type
    >::type  Space;
};

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template<>
class SharedAllocationRecord< Kokkos::HostSpace, void >
  : public SharedAllocationRecord< void, void >
{
private:
  friend Kokkos::HostSpace;

  typedef SharedAllocationRecord< void, void >  RecordBase;

  SharedAllocationRecord( const SharedAllocationRecord & ) = delete;
  SharedAllocationRecord & operator = ( const SharedAllocationRecord & ) = delete;

  static void deallocate( RecordBase * );

  /**\brief  Root record for tracked allocations from this HostSpace instance */
  static RecordBase s_root_record;

  const Kokkos::HostSpace m_space;

protected:
  ~SharedAllocationRecord();
  SharedAllocationRecord() = default;

  SharedAllocationRecord( const Kokkos::HostSpace        & arg_space
                        , const std::string              & arg_label
                        , const size_t                     arg_alloc_size
                        , const RecordBase::function_type  arg_dealloc = & deallocate
                        );

public:

  inline
  std::string get_label() const
  {
    return std::string( RecordBase::head()->m_label );
  }

  KOKKOS_INLINE_FUNCTION static
  SharedAllocationRecord * allocate( const Kokkos::HostSpace &  arg_space
                                   , const std::string       &  arg_label
                                   , const size_t               arg_alloc_size
                                   )
  {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    return new SharedAllocationRecord( arg_space, arg_label, arg_alloc_size );
#else
    return (SharedAllocationRecord *) 0;
#endif
  }


  /**\brief  Allocate tracked memory in the space */
  static
  void * allocate_tracked( const Kokkos::HostSpace & arg_space
                         , const std::string & arg_label
                         , const size_t arg_alloc_size );

  /**\brief  Reallocate tracked memory in the space */
  static
  void * reallocate_tracked( void * const arg_alloc_ptr
                           , const size_t arg_alloc_size );

  /**\brief  Deallocate tracked memory in the space */
  static
  void deallocate_tracked( void * const arg_alloc_ptr );

  static SharedAllocationRecord * get_record( void * arg_alloc_ptr );

  static void print_records( std::ostream &, const Kokkos::HostSpace &, bool detail = false );
};

} // namespace Impl

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template< class DstSpace, class SrcSpace, class ExecutionSpace = typename DstSpace::execution_space > struct DeepCopy;

template< class ExecutionSpace >
struct DeepCopy< HostSpace, HostSpace, ExecutionSpace > {
  DeepCopy( void * dst, const void * src, size_t n ) {
    memcpy( dst, src, n );
  }

  DeepCopy( const ExecutionSpace& exec, void * dst, const void * src, size_t n ) {
    exec.fence();
    memcpy( dst, src, n );
  }
};

} // namespace Impl

} // namespace Kokkos

#endif // #define KOKKOS_HOSTSPACE_HPP

