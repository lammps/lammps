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


#include <Kokkos_Macros.hpp>

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstring>

#include <iostream>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <Kokkos_HBWSpace.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Kokkos_Atomic.hpp>
#ifdef KOKKOS_ENABLE_HBWSPACE
#include <memkind.h>
#endif

#if defined(KOKKOS_ENABLE_PROFILING)
#include <impl/Kokkos_Profiling_Interface.hpp>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#ifdef KOKKOS_ENABLE_HBWSPACE
#define MEMKIND_TYPE MEMKIND_HBW //hbw_get_kind(HBW_PAGESIZE_4KB)

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {

/* Default allocation mechanism */
HBWSpace::HBWSpace()
  : m_alloc_mech(
     HBWSpace::STD_MALLOC
    )
{
printf("Init\n");
setenv("MEMKIND_HBW_NODES", "1", 0);
}

/* Default allocation mechanism */
HBWSpace::HBWSpace( const HBWSpace::AllocationMechanism & arg_alloc_mech )
  : m_alloc_mech( HBWSpace::STD_MALLOC )
{
printf("Init2\n");
setenv("MEMKIND_HBW_NODES", "1", 0);
  if ( arg_alloc_mech == STD_MALLOC ) {
    m_alloc_mech = HBWSpace::STD_MALLOC ;
  }
}

void * HBWSpace::allocate( const size_t arg_alloc_size ) const
{
  static_assert( sizeof(void*) == sizeof(uintptr_t)
               , "Error sizeof(void*) != sizeof(uintptr_t)" );

  static_assert( Kokkos::Impl::power_of_two< Kokkos::Impl::MEMORY_ALIGNMENT >::value
               , "Memory alignment must be power of two" );

  constexpr uintptr_t alignment = Kokkos::Impl::MEMORY_ALIGNMENT ;
  constexpr uintptr_t alignment_mask = alignment - 1 ;

  void * ptr = 0 ;

  if ( arg_alloc_size ) {

    if ( m_alloc_mech == STD_MALLOC ) {
      // Over-allocate to and round up to guarantee proper alignment.
      size_t size_padded = arg_alloc_size + sizeof(void*) + alignment ;

      void * alloc_ptr = memkind_malloc(MEMKIND_TYPE, size_padded );

      if (alloc_ptr) {
        uintptr_t address = reinterpret_cast<uintptr_t>(alloc_ptr);

        // offset enough to record the alloc_ptr
        address += sizeof(void *);
        uintptr_t rem = address % alignment;
        uintptr_t offset = rem ? (alignment - rem) : 0u;
        address += offset;
        ptr = reinterpret_cast<void *>(address);
        // record the alloc'd pointer
        address -= sizeof(void *);
        *reinterpret_cast<void **>(address) = alloc_ptr;
      }
    }
  }

  if ( ( ptr == 0 ) || ( reinterpret_cast<uintptr_t>(ptr) == ~uintptr_t(0) )
       || ( reinterpret_cast<uintptr_t>(ptr) & alignment_mask ) ) {
    std::ostringstream msg ;
    msg << "Kokkos::Experimental::HBWSpace::allocate[ " ;
    switch( m_alloc_mech ) {
    case STD_MALLOC: msg << "STD_MALLOC" ; break ;
    }
    msg << " ]( " << arg_alloc_size << " ) FAILED" ;
    if ( ptr == NULL ) { msg << " NULL" ; }
    else { msg << " NOT ALIGNED " << ptr ; }

    std::cerr << msg.str() << std::endl ;
    std::cerr.flush();

    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }

  return ptr;
}


void HBWSpace::deallocate( void * const arg_alloc_ptr , const size_t arg_alloc_size ) const
{
  if ( arg_alloc_ptr ) {

    if ( m_alloc_mech == STD_MALLOC ) {
      void * alloc_ptr = *(reinterpret_cast<void **>(arg_alloc_ptr) -1);
      memkind_free(MEMKIND_TYPE, alloc_ptr );
    }

  }
}

} // namespace Experimental
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

SharedAllocationRecord< void , void >
SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::s_root_record ;

void
SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::
deallocate( SharedAllocationRecord< void , void > * arg_rec )
{
  delete static_cast<SharedAllocationRecord*>(arg_rec);
}

SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::
~SharedAllocationRecord()
{
  #if defined(KOKKOS_ENABLE_PROFILING)
  if(Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::deallocateData(
      Kokkos::Profiling::SpaceHandle(Kokkos::Experimental::HBWSpace::name()),RecordBase::m_alloc_ptr->m_label,
      data(),size());
  }
  #endif

  m_space.deallocate( SharedAllocationRecord< void , void >::m_alloc_ptr
                    , SharedAllocationRecord< void , void >::m_alloc_size
                    );
}

SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::
SharedAllocationRecord( const Kokkos::Experimental::HBWSpace & arg_space
                      , const std::string       & arg_label
                      , const size_t              arg_alloc_size
                      , const SharedAllocationRecord< void , void >::function_type arg_dealloc
                      )
  // Pass through allocated [ SharedAllocationHeader , user_memory ]
  // Pass through deallocation function
  : SharedAllocationRecord< void , void >
      ( & SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::s_root_record
      , reinterpret_cast<SharedAllocationHeader*>( arg_space.allocate( sizeof(SharedAllocationHeader) + arg_alloc_size ) )
      , sizeof(SharedAllocationHeader) + arg_alloc_size
      , arg_dealloc
      )
  , m_space( arg_space )
{
  #if defined(KOKKOS_ENABLE_PROFILING)
  if(Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::allocateData(Kokkos::Profiling::SpaceHandle(arg_space.name()),arg_label,data(),arg_alloc_size);
  }
  #endif

  // Fill in the Header information
  RecordBase::m_alloc_ptr->m_record = static_cast< SharedAllocationRecord< void , void > * >( this );

  strncpy( RecordBase::m_alloc_ptr->m_label
          , arg_label.c_str()
          , SharedAllocationHeader::maximum_label_length
          );
}

//----------------------------------------------------------------------------

void * SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::
allocate_tracked( const Kokkos::Experimental::HBWSpace & arg_space
                , const std::string & arg_alloc_label
                , const size_t arg_alloc_size )
{
  if ( ! arg_alloc_size ) return (void *) 0 ;

  SharedAllocationRecord * const r =
    allocate( arg_space , arg_alloc_label , arg_alloc_size );

  RecordBase::increment( r );

  return r->data();
}

void SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::
deallocate_tracked( void * const arg_alloc_ptr )
{
  if ( arg_alloc_ptr != 0 ) {
    SharedAllocationRecord * const r = get_record( arg_alloc_ptr );

    RecordBase::decrement( r );
  }
}

void * SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::
reallocate_tracked( void * const arg_alloc_ptr
                  , const size_t arg_alloc_size )
{
  SharedAllocationRecord * const r_old = get_record( arg_alloc_ptr );
  SharedAllocationRecord * const r_new = allocate( r_old->m_space , r_old->get_label() , arg_alloc_size );

  Kokkos::Impl::DeepCopy<Kokkos::Experimental::HBWSpace,Kokkos::Experimental::HBWSpace>( r_new->data() , r_old->data()
                                             , std::min( r_old->size() , r_new->size() ) );

  RecordBase::increment( r_new );
  RecordBase::decrement( r_old );

  return r_new->data();
}

SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void > *
SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::get_record( void * alloc_ptr )
{
  typedef SharedAllocationHeader  Header ;
  typedef SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >  RecordHost ;

  SharedAllocationHeader const * const head   = alloc_ptr ? Header::get_header( alloc_ptr ) : (SharedAllocationHeader *)0 ;
  RecordHost                   * const record = head ? static_cast< RecordHost * >( head->m_record ) : (RecordHost *) 0 ;

  if ( ! alloc_ptr || record->m_alloc_ptr != head ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::Impl::SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::get_record ERROR" ) );
  }

  return record ;
}

// Iterate records to print orphaned memory ...
void SharedAllocationRecord< Kokkos::Experimental::HBWSpace , void >::
print_records( std::ostream & s , const Kokkos::Experimental::HBWSpace & space , bool detail )
{
  SharedAllocationRecord< void , void >::print_host_accessible_records( s , "HBWSpace" , & s_root_record , detail );
}

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {
namespace {
  const unsigned HBW_SPACE_ATOMIC_MASK = 0xFFFF;
  const unsigned HBW_SPACE_ATOMIC_XOR_MASK = 0x5A39;
  static int HBW_SPACE_ATOMIC_LOCKS[HBW_SPACE_ATOMIC_MASK+1];
}

namespace Impl {
void init_lock_array_hbw_space() {
  static int is_initialized = 0;
  if(! is_initialized)
    for(int i = 0; i < static_cast<int> (HBW_SPACE_ATOMIC_MASK+1); i++)
      HBW_SPACE_ATOMIC_LOCKS[i] = 0;
}

bool lock_address_hbw_space(void* ptr) {
  return 0 == atomic_compare_exchange( &HBW_SPACE_ATOMIC_LOCKS[
      (( size_t(ptr) >> 2 ) & HBW_SPACE_ATOMIC_MASK) ^ HBW_SPACE_ATOMIC_XOR_MASK] ,
                                  0 , 1);
}

void unlock_address_hbw_space(void* ptr) {
   atomic_exchange( &HBW_SPACE_ATOMIC_LOCKS[
      (( size_t(ptr) >> 2 ) & HBW_SPACE_ATOMIC_MASK) ^ HBW_SPACE_ATOMIC_XOR_MASK] ,
                    0);
}

}
}
}
#endif

