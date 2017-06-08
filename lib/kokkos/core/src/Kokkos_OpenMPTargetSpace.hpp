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

#ifndef KOKKOS_OPENMPTARGETSPACE_HPP
#define KOKKOS_OPENMPTARGETSPACE_HPP

#include <cstring>
#include <string>
#include <iosfwd>
#include <typeinfo>

#include <Kokkos_Core_fwd.hpp>

#ifdef KOKKOS_ENABLE_OPENMPTARGET

#include <Kokkos_HostSpace.hpp>
#include <omp.h>
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

/// \brief Initialize lock array for arbitrary size atomics.
///
/// Arbitrary atomics are implemented using a hash table of locks
/// where the hash value is derived from the address of the
/// object for which an atomic operation is performed.
/// This function initializes the locks to zero (unset).
//void init_lock_array_host_space();

/// \brief Aquire a lock for the address
///
/// This function tries to aquire the lock for the hash value derived
/// from the provided ptr. If the lock is successfully aquired the
/// function returns true. Otherwise it returns false.
//bool lock_address_host_space(void* ptr);

/// \brief Release lock for the address
///
/// This function releases the lock for the hash value derived
/// from the provided ptr. This function should only be called
/// after previously successfully aquiring a lock with
/// lock_address.
//void unlock_address_host_space(void* ptr);

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {
namespace Experimental {

/// \class OpenMPTargetSpace
/// \brief Memory management for host memory.
///
/// OpenMPTargetSpace is a memory space that governs host memory.  "Host"
/// memory means the usual CPU-accessible memory.
class OpenMPTargetSpace {
public:

  //! Tag this class as a kokkos memory space
  typedef OpenMPTargetSpace  memory_space ;
  typedef size_t     size_type ;

  /// \typedef execution_space
  /// \brief Default execution space for this memory space.
  ///
  /// Every memory space has a default execution space.  This is
  /// useful for things like initializing a View (which happens in
  /// parallel using the View's default execution space).
  typedef Kokkos::Experimental::OpenMPTarget   execution_space ;

  //! This memory space preferred device_type
  typedef Kokkos::Device<execution_space,memory_space> device_type;

  /*--------------------------------*/

  /**\brief  Default memory space instance */
  OpenMPTargetSpace();
  OpenMPTargetSpace( OpenMPTargetSpace && rhs ) = default ;
  OpenMPTargetSpace( const OpenMPTargetSpace & rhs ) = default ;
  OpenMPTargetSpace & operator = ( OpenMPTargetSpace && ) = default ;
  OpenMPTargetSpace & operator = ( const OpenMPTargetSpace & ) = default ;
  ~OpenMPTargetSpace() = default ;

  /**\brief  Allocate untracked memory in the space */
  void * allocate( const size_t arg_alloc_size ) const ;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate( void * const arg_alloc_ptr 
                 , const size_t arg_alloc_size ) const ;

private:

  friend class Kokkos::Impl::SharedAllocationRecord< Kokkos::Experimental::OpenMPTargetSpace , void > ;
};
} // namespace Experimental
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
class SharedAllocationRecord< Kokkos::Experimental::OpenMPTargetSpace , void >
  : public SharedAllocationRecord< void , void >
{
private:

  friend Kokkos::Experimental::OpenMPTargetSpace ;

  typedef SharedAllocationRecord< void , void >  RecordBase ;

  SharedAllocationRecord( const SharedAllocationRecord & ) = delete ;
  SharedAllocationRecord & operator = ( const SharedAllocationRecord & ) = delete ;

  static void deallocate( RecordBase * );

  /**\brief  Root record for tracked allocations from this OpenMPTargetSpace instance */
  static RecordBase s_root_record ;

  const Kokkos::Experimental::OpenMPTargetSpace m_space ;

protected:

  ~SharedAllocationRecord();
  SharedAllocationRecord() = default ;

  SharedAllocationRecord( const Kokkos::Experimental::OpenMPTargetSpace        & arg_space
                        , const std::string              & arg_label
                        , const size_t                     arg_alloc_size
                        , const RecordBase::function_type  arg_dealloc = & deallocate
                        );

public:

  std::string get_label() const;

  KOKKOS_INLINE_FUNCTION static
  SharedAllocationRecord * allocate( const Kokkos::Experimental::OpenMPTargetSpace &  arg_space
                                   , const std::string       &  arg_label
                                   , const size_t               arg_alloc_size
                                   );

  /**\brief  Allocate tracked memory in the space */
  static
  void * allocate_tracked( const Kokkos::Experimental::OpenMPTargetSpace & arg_space
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

  static void print_records( std::ostream & , const Kokkos::Experimental::OpenMPTargetSpace & , bool detail = false );
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//TODO: implement all possible deep_copies
template<class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::OpenMPTargetSpace,Kokkos::Experimental::OpenMPTargetSpace,ExecutionSpace> {
  DeepCopy( void * dst , const void * src , size_t n ) {
    omp_target_memcpy( dst , const_cast<void*> (src) , n, 0, 0, omp_get_default_device(), omp_get_default_device());
  }
  DeepCopy( const ExecutionSpace& exec, void * dst , const void * src , size_t n ) {
    exec.fence();
    omp_target_memcpy( dst , const_cast<void*> (src) , n, 0, 0, omp_get_default_device(), omp_get_default_device());
  }
};


template<class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::OpenMPTargetSpace,HostSpace,ExecutionSpace> {
  DeepCopy( void * dst , const void * src , size_t n ) {
    omp_target_memcpy( dst , const_cast<void*> (src) , n, 0, 0, omp_get_default_device(), omp_get_initial_device());
  }
  DeepCopy( const ExecutionSpace& exec, void * dst , const void * src , size_t n ) {
    exec.fence();
    omp_target_memcpy( dst , const_cast<void*> (src) , n, 0, 0, omp_get_default_device(), omp_get_initial_device());
  }
};

template<class ExecutionSpace>
struct DeepCopy<HostSpace,Kokkos::Experimental::OpenMPTargetSpace,ExecutionSpace> {
  DeepCopy( void * dst , const void * src , size_t n ) {
    omp_target_memcpy( dst , const_cast<void*> (src) , n, 0, 0, omp_get_initial_device(), omp_get_default_device());
  }
  DeepCopy( const ExecutionSpace& exec, void * dst , const void * src , size_t n ) {
    exec.fence();
    omp_target_memcpy( dst , const_cast<void*> (src) , n, 0, 0, omp_get_initial_device(), omp_get_default_device());
  }
};


template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::HostSpace , Kokkos::Experimental::OpenMPTargetSpace >
{
  enum { value = false };
  inline static void verify( void ) { }
  inline static void verify( const void * ) { }
};

} // namespace Impl
} // namespace Kokkos

#endif
#endif /* #define KOKKOS_OPENMPTARGETSPACE_HPP */

