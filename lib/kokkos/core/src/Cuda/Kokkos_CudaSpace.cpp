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

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <Kokkos_Macros.hpp>

/* only compile this file if CUDA is enabled for Kokkos */
#ifdef KOKKOS_HAVE_CUDA

#include <Kokkos_Cuda.hpp>
#include <Kokkos_CudaSpace.hpp>

#include <Cuda/Kokkos_Cuda_BasicAllocators.hpp>
#include <Cuda/Kokkos_Cuda_Internal.hpp>
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

namespace {
   cudaStream_t get_deep_copy_stream() {
     static cudaStream_t s = 0;
     if( s == 0) {
       cudaStreamCreate ( &s );
     }
     return s;
   }
}

DeepCopy<CudaSpace,CudaSpace,Cuda>::DeepCopy( void * dst , const void * src , size_t n )
{ CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) ); }

DeepCopy<HostSpace,CudaSpace,Cuda>::DeepCopy( void * dst , const void * src , size_t n )
{ CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) ); }

DeepCopy<CudaSpace,HostSpace,Cuda>::DeepCopy( void * dst , const void * src , size_t n )
{ CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) ); }

DeepCopy<CudaSpace,CudaSpace,Cuda>::DeepCopy( const Cuda & instance , void * dst , const void * src , size_t n )
{ CUDA_SAFE_CALL( cudaMemcpyAsync( dst , src , n , cudaMemcpyDefault , instance.cuda_stream() ) ); }

DeepCopy<HostSpace,CudaSpace,Cuda>::DeepCopy( const Cuda & instance , void * dst , const void * src , size_t n )
{ CUDA_SAFE_CALL( cudaMemcpyAsync( dst , src , n , cudaMemcpyDefault , instance.cuda_stream() ) ); }

DeepCopy<CudaSpace,HostSpace,Cuda>::DeepCopy( const Cuda & instance , void * dst , const void * src , size_t n )
{ CUDA_SAFE_CALL( cudaMemcpyAsync( dst , src , n , cudaMemcpyDefault , instance.cuda_stream() ) ); }

void DeepCopyAsyncCuda( void * dst , const void * src , size_t n) {
  cudaStream_t s = get_deep_copy_stream();
  CUDA_SAFE_CALL( cudaMemcpyAsync( dst , src , n , cudaMemcpyDefault , s ) );
  cudaStreamSynchronize(s);
}

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/


namespace Kokkos {

#if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

namespace {

void texture_object_attach_impl(  Impl::AllocationTracker const & tracker
                                , unsigned type_size
                                , ::cudaChannelFormatDesc const & desc
                               )
{
  enum { TEXTURE_BOUND_1D = 2u << 27 };

  if ( tracker.attribute() == NULL ) {
    // check for correct allocator
    const bool ok_alloc =  tracker.allocator()->support_texture_binding();

    const bool ok_count = (tracker.alloc_size() / type_size) < TEXTURE_BOUND_1D;

    if (ok_alloc && ok_count) {
      Impl::TextureAttribute * attr = new Impl::TextureAttribute( tracker.alloc_ptr(), tracker.alloc_size(), desc );
      tracker.set_attribute( attr );
    }
    else {
      std::ostringstream oss;
      oss << "Error: Cannot attach texture object";
      if (!ok_alloc) {
        oss << ", incompatabile allocator " << tracker.allocator()->name();
      }
      if (!ok_count) {
        oss << ", array " << tracker.label() << " too large";
      }
      oss << ".";
      Kokkos::Impl::throw_runtime_exception( oss.str() );
    }
  }

  if ( NULL == dynamic_cast<Impl::TextureAttribute *>(tracker.attribute()) ) {
    std::ostringstream oss;
    oss << "Error: Allocation " << tracker.label() << " already has an attribute attached.";
    Kokkos::Impl::throw_runtime_exception( oss.str() );
  }

}

} // unnamed namespace

/*--------------------------------------------------------------------------*/

Impl::AllocationTracker CudaSpace::allocate_and_track( const std::string & label, const size_t size )
{
  return Impl::AllocationTracker( allocator(), size, label);
}

void CudaSpace::texture_object_attach(  Impl::AllocationTracker const & tracker
                                      , unsigned type_size
                                      , ::cudaChannelFormatDesc const & desc
                                     )
{
  texture_object_attach_impl( tracker, type_size, desc );
}

#endif /* #if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW ) */

void CudaSpace::access_error()
{
  const std::string msg("Kokkos::CudaSpace::access_error attempt to execute Cuda function from non-Cuda space" );
  Kokkos::Impl::throw_runtime_exception( msg );
}

void CudaSpace::access_error( const void * const )
{
  const std::string msg("Kokkos::CudaSpace::access_error attempt to execute Cuda function from non-Cuda space" );
  Kokkos::Impl::throw_runtime_exception( msg );
}

/*--------------------------------------------------------------------------*/

#if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

Impl::AllocationTracker CudaUVMSpace::allocate_and_track( const std::string & label, const size_t size )
{
  return Impl::AllocationTracker( allocator(), size, label);
}

void CudaUVMSpace::texture_object_attach(  Impl::AllocationTracker const & tracker
                                         , unsigned type_size
                                         , ::cudaChannelFormatDesc const & desc
                                        )
{
  texture_object_attach_impl( tracker, type_size, desc );
}

#endif /* #if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW ) */

bool CudaUVMSpace::available()
{
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION ) && !defined(__APPLE__)
  enum { UVM_available = true };
#else
  enum { UVM_available = false };
#endif
  return UVM_available;
}

/*--------------------------------------------------------------------------*/

#if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

Impl::AllocationTracker CudaHostPinnedSpace::allocate_and_track( const std::string & label, const size_t size )
{
  return Impl::AllocationTracker( allocator(), size, label);
}

#endif /* #if ! defined( KOKKOS_USING_EXPERIMENTAL_VIEW ) */

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

CudaSpace::CudaSpace()
  : m_device( Kokkos::Cuda().cuda_device() )
{
}

CudaUVMSpace::CudaUVMSpace()
  : m_device( Kokkos::Cuda().cuda_device() )
{
}

CudaHostPinnedSpace::CudaHostPinnedSpace()
{
}

void * CudaSpace::allocate( const size_t arg_alloc_size ) const
{
  void * ptr = NULL;

  CUDA_SAFE_CALL( cudaMalloc( &ptr, arg_alloc_size ) );

  return ptr ;
}

void * CudaUVMSpace::allocate( const size_t arg_alloc_size ) const
{
  void * ptr = NULL;

  CUDA_SAFE_CALL( cudaMallocManaged( &ptr, arg_alloc_size , cudaMemAttachGlobal ) );

  return ptr ;
}

void * CudaHostPinnedSpace::allocate( const size_t arg_alloc_size ) const
{
  void * ptr = NULL;

  CUDA_SAFE_CALL( cudaHostAlloc( &ptr, arg_alloc_size , cudaHostAllocDefault ) );

  return ptr ;
}

void CudaSpace::deallocate( void * const arg_alloc_ptr , const size_t /* arg_alloc_size */ ) const
{
  try {
    CUDA_SAFE_CALL( cudaFree( arg_alloc_ptr ) );
  } catch(...) {}
}

void CudaUVMSpace::deallocate( void * const arg_alloc_ptr , const size_t /* arg_alloc_size */ ) const
{
  try {
    CUDA_SAFE_CALL( cudaFree( arg_alloc_ptr ) );
  } catch(...) {}
}

void CudaHostPinnedSpace::deallocate( void * const arg_alloc_ptr , const size_t /* arg_alloc_size */ ) const
{
  try {
    CUDA_SAFE_CALL( cudaFreeHost( arg_alloc_ptr ) );
  } catch(...) {}
}

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

SharedAllocationRecord< void , void >
SharedAllocationRecord< Kokkos::CudaSpace , void >::s_root_record ;

SharedAllocationRecord< void , void >
SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::s_root_record ;

SharedAllocationRecord< void , void >
SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::s_root_record ;

::cudaTextureObject_t
SharedAllocationRecord< Kokkos::CudaSpace , void >::
attach_texture_object( const unsigned sizeof_alias
                     , void *   const alloc_ptr
                     , size_t   const alloc_size )
{
  enum { TEXTURE_BOUND_1D = 1u << 27 };

  if ( ( alloc_ptr == 0 ) || ( sizeof_alias * TEXTURE_BOUND_1D <= alloc_size ) ) {
    std::ostringstream msg ;
    msg << "Kokkos::CudaSpace ERROR: Cannot attach texture object to"
        << " alloc_ptr(" << alloc_ptr << ")"
        << " alloc_size(" << alloc_size << ")"
        << " max_size(" << ( sizeof_alias * TEXTURE_BOUND_1D ) << ")" ;
    std::cerr << msg.str() << std::endl ;
    std::cerr.flush();
    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }

  ::cudaTextureObject_t tex_obj ;

  struct cudaResourceDesc resDesc ;
  struct cudaTextureDesc  texDesc ;

  memset( & resDesc , 0 , sizeof(resDesc) );
  memset( & texDesc , 0 , sizeof(texDesc) );

  resDesc.resType                = cudaResourceTypeLinear ;
  resDesc.res.linear.desc        = ( sizeof_alias ==  4 ?  cudaCreateChannelDesc< int >() :
                                   ( sizeof_alias ==  8 ?  cudaCreateChannelDesc< ::int2 >() :
                                  /* sizeof_alias == 16 */ cudaCreateChannelDesc< ::int4 >() ) );
  resDesc.res.linear.sizeInBytes = alloc_size ;
  resDesc.res.linear.devPtr      = alloc_ptr ;

  CUDA_SAFE_CALL( cudaCreateTextureObject( & tex_obj , & resDesc, & texDesc, NULL ) );

  return tex_obj ;
}

std::string
SharedAllocationRecord< Kokkos::CudaSpace , void >::get_label() const
{
  SharedAllocationHeader header ;

  Kokkos::Impl::DeepCopy< Kokkos::HostSpace , Kokkos::CudaSpace >( & header , RecordBase::head() , sizeof(SharedAllocationHeader) );

  return std::string( header.m_label );
}

std::string
SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::get_label() const
{
  return std::string( RecordBase::head()->m_label );
}

std::string
SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::get_label() const
{
  return std::string( RecordBase::head()->m_label );
}

SharedAllocationRecord< Kokkos::CudaSpace , void > *
SharedAllocationRecord< Kokkos::CudaSpace , void >::
allocate( const Kokkos::CudaSpace &  arg_space
        , const std::string       &  arg_label
        , const size_t               arg_alloc_size
        )
{
  return new SharedAllocationRecord( arg_space , arg_label , arg_alloc_size );
}

SharedAllocationRecord< Kokkos::CudaUVMSpace , void > *
SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::
allocate( const Kokkos::CudaUVMSpace &  arg_space
        , const std::string          &  arg_label
        , const size_t                  arg_alloc_size
        )
{
  return new SharedAllocationRecord( arg_space , arg_label , arg_alloc_size );
}

SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void > *
SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::
allocate( const Kokkos::CudaHostPinnedSpace &  arg_space
        , const std::string                 &  arg_label
        , const size_t                         arg_alloc_size
        )
{
  return new SharedAllocationRecord( arg_space , arg_label , arg_alloc_size );
}

void
SharedAllocationRecord< Kokkos::CudaSpace , void >::
deallocate( SharedAllocationRecord< void , void > * arg_rec )
{
  delete static_cast<SharedAllocationRecord*>(arg_rec);
}

void
SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::
deallocate( SharedAllocationRecord< void , void > * arg_rec )
{
  delete static_cast<SharedAllocationRecord*>(arg_rec);
}

void
SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::
deallocate( SharedAllocationRecord< void , void > * arg_rec )
{
  delete static_cast<SharedAllocationRecord*>(arg_rec);
}

SharedAllocationRecord< Kokkos::CudaSpace , void >::
~SharedAllocationRecord()
{
  m_space.deallocate( SharedAllocationRecord< void , void >::m_alloc_ptr
                    , SharedAllocationRecord< void , void >::m_alloc_size
                    );
}

SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::
~SharedAllocationRecord()
{
  m_space.deallocate( SharedAllocationRecord< void , void >::m_alloc_ptr
                    , SharedAllocationRecord< void , void >::m_alloc_size
                    );
}

SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::
~SharedAllocationRecord()
{
  m_space.deallocate( SharedAllocationRecord< void , void >::m_alloc_ptr
                    , SharedAllocationRecord< void , void >::m_alloc_size
                    );
}

SharedAllocationRecord< Kokkos::CudaSpace , void >::
SharedAllocationRecord( const Kokkos::CudaSpace & arg_space
                      , const std::string       & arg_label
                      , const size_t              arg_alloc_size
                      , const SharedAllocationRecord< void , void >::function_type arg_dealloc
                      )
  // Pass through allocated [ SharedAllocationHeader , user_memory ]
  // Pass through deallocation function
  : SharedAllocationRecord< void , void >
      ( & SharedAllocationRecord< Kokkos::CudaSpace , void >::s_root_record
      , reinterpret_cast<SharedAllocationHeader*>( arg_space.allocate( sizeof(SharedAllocationHeader) + arg_alloc_size ) )
      , sizeof(SharedAllocationHeader) + arg_alloc_size
      , arg_dealloc
      )
  , m_tex_obj( 0 )
  , m_space( arg_space )
{
  SharedAllocationHeader header ;

  // Fill in the Header information
  header.m_record = static_cast< SharedAllocationRecord< void , void > * >( this );

  strncpy( header.m_label
          , arg_label.c_str()
          , SharedAllocationHeader::maximum_label_length
          );

  // Copy to device memory
  Kokkos::Impl::DeepCopy<CudaSpace,HostSpace>::DeepCopy( RecordBase::m_alloc_ptr , & header , sizeof(SharedAllocationHeader) );
}

SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::
SharedAllocationRecord( const Kokkos::CudaUVMSpace & arg_space
                      , const std::string          & arg_label
                      , const size_t                 arg_alloc_size
                      , const SharedAllocationRecord< void , void >::function_type arg_dealloc
                      )
  // Pass through allocated [ SharedAllocationHeader , user_memory ]
  // Pass through deallocation function
  : SharedAllocationRecord< void , void >
      ( & SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::s_root_record
      , reinterpret_cast<SharedAllocationHeader*>( arg_space.allocate( sizeof(SharedAllocationHeader) + arg_alloc_size ) )
      , sizeof(SharedAllocationHeader) + arg_alloc_size
      , arg_dealloc
      )
  , m_tex_obj( 0 )
  , m_space( arg_space )
{
  // Fill in the Header information, directly accessible via UVM

  RecordBase::m_alloc_ptr->m_record = this ;

  strncpy( RecordBase::m_alloc_ptr->m_label
          , arg_label.c_str()
          , SharedAllocationHeader::maximum_label_length
          );
}

SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::
SharedAllocationRecord( const Kokkos::CudaHostPinnedSpace & arg_space
                      , const std::string                 & arg_label
                      , const size_t                        arg_alloc_size
                      , const SharedAllocationRecord< void , void >::function_type arg_dealloc
                      )
  // Pass through allocated [ SharedAllocationHeader , user_memory ]
  // Pass through deallocation function
  : SharedAllocationRecord< void , void >
      ( & SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::s_root_record
      , reinterpret_cast<SharedAllocationHeader*>( arg_space.allocate( sizeof(SharedAllocationHeader) + arg_alloc_size ) )
      , sizeof(SharedAllocationHeader) + arg_alloc_size
      , arg_dealloc
      )
  , m_space( arg_space )
{
  // Fill in the Header information, directly accessible via UVM

  RecordBase::m_alloc_ptr->m_record = this ;

  strncpy( RecordBase::m_alloc_ptr->m_label
          , arg_label.c_str()
          , SharedAllocationHeader::maximum_label_length
          );
}

//----------------------------------------------------------------------------

void * SharedAllocationRecord< Kokkos::CudaSpace , void >::
allocate_tracked( const Kokkos::CudaSpace & arg_space
                , const std::string & arg_alloc_label
                , const size_t arg_alloc_size )
{
  if ( ! arg_alloc_size ) return (void *) 0 ;

  SharedAllocationRecord * const r =
    allocate( arg_space , arg_alloc_label , arg_alloc_size );

  RecordBase::increment( r );

  return r->data();
}

void SharedAllocationRecord< Kokkos::CudaSpace , void >::
deallocate_tracked( void * const arg_alloc_ptr )
{
  if ( arg_alloc_ptr != 0 ) {
    SharedAllocationRecord * const r = get_record( arg_alloc_ptr );

    RecordBase::decrement( r );
  }
}

void * SharedAllocationRecord< Kokkos::CudaSpace , void >::
reallocate_tracked( void * const arg_alloc_ptr
                  , const size_t arg_alloc_size )
{
  SharedAllocationRecord * const r_old = get_record( arg_alloc_ptr );
  SharedAllocationRecord * const r_new = allocate( r_old->m_space , r_old->get_label() , arg_alloc_size );

  Kokkos::Impl::DeepCopy<CudaSpace,CudaSpace>( r_new->data() , r_old->data()
                                             , std::min( r_old->size() , r_new->size() ) );

  RecordBase::increment( r_new );
  RecordBase::decrement( r_old );

  return r_new->data();
}

void * SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::
allocate_tracked( const Kokkos::CudaUVMSpace & arg_space
                , const std::string & arg_alloc_label
                , const size_t arg_alloc_size )
{
  if ( ! arg_alloc_size ) return (void *) 0 ;

  SharedAllocationRecord * const r =
    allocate( arg_space , arg_alloc_label , arg_alloc_size );

  RecordBase::increment( r );

  return r->data();
}

void SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::
deallocate_tracked( void * const arg_alloc_ptr )
{
  if ( arg_alloc_ptr != 0 ) {
    SharedAllocationRecord * const r = get_record( arg_alloc_ptr );

    RecordBase::decrement( r );
  }
}

void * SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::
reallocate_tracked( void * const arg_alloc_ptr
                  , const size_t arg_alloc_size )
{
  SharedAllocationRecord * const r_old = get_record( arg_alloc_ptr );
  SharedAllocationRecord * const r_new = allocate( r_old->m_space , r_old->get_label() , arg_alloc_size );

  Kokkos::Impl::DeepCopy<CudaUVMSpace,CudaUVMSpace>( r_new->data() , r_old->data()
                                             , std::min( r_old->size() , r_new->size() ) );

  RecordBase::increment( r_new );
  RecordBase::decrement( r_old );

  return r_new->data();
}

void * SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::
allocate_tracked( const Kokkos::CudaHostPinnedSpace & arg_space
                , const std::string & arg_alloc_label
                , const size_t arg_alloc_size )
{
  if ( ! arg_alloc_size ) return (void *) 0 ;

  SharedAllocationRecord * const r =
    allocate( arg_space , arg_alloc_label , arg_alloc_size );

  RecordBase::increment( r );

  return r->data();
}

void SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::
deallocate_tracked( void * const arg_alloc_ptr )
{
  if ( arg_alloc_ptr != 0 ) {
    SharedAllocationRecord * const r = get_record( arg_alloc_ptr );

    RecordBase::decrement( r );
  }
}

void * SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::
reallocate_tracked( void * const arg_alloc_ptr
                  , const size_t arg_alloc_size )
{
  SharedAllocationRecord * const r_old = get_record( arg_alloc_ptr );
  SharedAllocationRecord * const r_new = allocate( r_old->m_space , r_old->get_label() , arg_alloc_size );

  Kokkos::Impl::DeepCopy<CudaHostPinnedSpace,CudaHostPinnedSpace>( r_new->data() , r_old->data()
                                             , std::min( r_old->size() , r_new->size() ) );

  RecordBase::increment( r_new );
  RecordBase::decrement( r_old );

  return r_new->data();
}

//----------------------------------------------------------------------------

SharedAllocationRecord< Kokkos::CudaSpace , void > *
SharedAllocationRecord< Kokkos::CudaSpace , void >::get_record( void * alloc_ptr )
{
  using Header     = SharedAllocationHeader ;
  using RecordBase = SharedAllocationRecord< void , void > ;
  using RecordCuda = SharedAllocationRecord< Kokkos::CudaSpace , void > ;

#if 0
  // Copy the header from the allocation
  Header head ;

  Header const * const head_cuda = alloc_ptr ? Header::get_header( alloc_ptr ) : (Header*) 0 ;

  if ( alloc_ptr ) {
    Kokkos::Impl::DeepCopy<HostSpace,CudaSpace>::DeepCopy( & head , head_cuda , sizeof(SharedAllocationHeader) );
  }

  RecordCuda * const record = alloc_ptr ? static_cast< RecordCuda * >( head.m_record ) : (RecordCuda *) 0 ;

  if ( ! alloc_ptr || record->m_alloc_ptr != head_cuda ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::Experimental::Impl::SharedAllocationRecord< Kokkos::CudaSpace , void >::get_record ERROR" ) );
  }

#else

  // Iterate the list to search for the record among all allocations
  // requires obtaining the root of the list and then locking the list.

  RecordCuda * const record = static_cast< RecordCuda * >( RecordBase::find( & s_root_record , alloc_ptr ) );

  if ( record == 0 ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::Experimental::Impl::SharedAllocationRecord< Kokkos::CudaSpace , void >::get_record ERROR" ) );
  }

#endif

  return record ;
}

SharedAllocationRecord< Kokkos::CudaUVMSpace , void > *
SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::get_record( void * alloc_ptr )
{
  using Header     = SharedAllocationHeader ;
  using RecordCuda = SharedAllocationRecord< Kokkos::CudaUVMSpace , void > ;

  Header * const h = alloc_ptr ? reinterpret_cast< Header * >( alloc_ptr ) - 1 : (Header *) 0 ;

  if ( ! alloc_ptr || h->m_record->m_alloc_ptr != h ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::Experimental::Impl::SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::get_record ERROR" ) );
  }

  return static_cast< RecordCuda * >( h->m_record );
}

SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void > *
SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::get_record( void * alloc_ptr )
{
  using Header     = SharedAllocationHeader ;
  using RecordCuda = SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void > ;

  Header * const h = alloc_ptr ? reinterpret_cast< Header * >( alloc_ptr ) - 1 : (Header *) 0 ;

  if ( ! alloc_ptr || h->m_record->m_alloc_ptr != h ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::Experimental::Impl::SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::get_record ERROR" ) );
  }

  return static_cast< RecordCuda * >( h->m_record );
}

// Iterate records to print orphaned memory ...
void
SharedAllocationRecord< Kokkos::CudaSpace , void >::
print_records( std::ostream & s , const Kokkos::CudaSpace & space , bool detail )
{
  SharedAllocationRecord< void , void > * r = & s_root_record ;

  char buffer[256] ;

  SharedAllocationHeader head ;

  if ( detail ) {
    do {
      if ( r->m_alloc_ptr ) {
        Kokkos::Impl::DeepCopy<HostSpace,CudaSpace>::DeepCopy( & head , r->m_alloc_ptr , sizeof(SharedAllocationHeader) );
      }
      else {
        head.m_label[0] = 0 ;
      }

      //Formatting dependent on sizeof(uintptr_t)
      const char * format_string;

      if (sizeof(uintptr_t) == sizeof(unsigned long)) { 
        format_string = "Cuda addr( 0x%.12lx ) list( 0x%.12lx 0x%.12lx ) extent[ 0x%.12lx + %.8ld ] count(%d) dealloc(0x%.12lx) %s\n";
      }
      else if (sizeof(uintptr_t) == sizeof(unsigned long long)) { 
        format_string = "Cuda addr( 0x%.12llx ) list( 0x%.12llx 0x%.12llx ) extent[ 0x%.12llx + %.8ld ] count(%d) dealloc(0x%.12llx) %s\n";
      }

      snprintf( buffer , 256 
              , format_string
              , reinterpret_cast<uintptr_t>( r )
              , reinterpret_cast<uintptr_t>( r->m_prev )
              , reinterpret_cast<uintptr_t>( r->m_next )
              , reinterpret_cast<uintptr_t>( r->m_alloc_ptr )
              , r->m_alloc_size
              , r->m_count
              , reinterpret_cast<uintptr_t>( r->m_dealloc )
              , head.m_label
              );
      std::cout << buffer ;
      r = r->m_next ;
    } while ( r != & s_root_record );
  }
  else {
    do {
      if ( r->m_alloc_ptr ) {

        Kokkos::Impl::DeepCopy<HostSpace,CudaSpace>::DeepCopy( & head , r->m_alloc_ptr , sizeof(SharedAllocationHeader) );

        //Formatting dependent on sizeof(uintptr_t)
        const char * format_string;

        if (sizeof(uintptr_t) == sizeof(unsigned long)) { 
          format_string = "Cuda [ 0x%.12lx + %ld ] %s\n";
        }
        else if (sizeof(uintptr_t) == sizeof(unsigned long long)) { 
          format_string = "Cuda [ 0x%.12llx + %ld ] %s\n";
        }

        snprintf( buffer , 256 
                , format_string
                , reinterpret_cast< uintptr_t >( r->data() )
                , r->size()
                , head.m_label
                );
      }
      else {
        snprintf( buffer , 256 , "Cuda [ 0 + 0 ]\n" );
      }
      std::cout << buffer ;
      r = r->m_next ;
    } while ( r != & s_root_record );
  }
}

void
SharedAllocationRecord< Kokkos::CudaUVMSpace , void >::
print_records( std::ostream & s , const Kokkos::CudaUVMSpace & space , bool detail )
{
  SharedAllocationRecord< void , void >::print_host_accessible_records( s , "CudaUVM" , & s_root_record , detail );
}

void
SharedAllocationRecord< Kokkos::CudaHostPinnedSpace , void >::
print_records( std::ostream & s , const Kokkos::CudaHostPinnedSpace & space , bool detail )
{
  SharedAllocationRecord< void , void >::print_host_accessible_records( s , "CudaHostPinned" , & s_root_record , detail );
}

} // namespace Impl
} // namespace Experimental
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace {
  __global__ void init_lock_array_kernel() {
    unsigned i = blockIdx.x*blockDim.x + threadIdx.x;

    if(i<CUDA_SPACE_ATOMIC_MASK+1)
      kokkos_impl_cuda_atomic_lock_array[i] = 0;
  }
}

namespace Impl {
int* lock_array_cuda_space_ptr(bool deallocate) {
  static int* ptr = NULL;
  if(deallocate) {
    cudaFree(ptr);
    ptr = NULL;
  }

  if(ptr==NULL && !deallocate)
    cudaMalloc(&ptr,sizeof(int)*(CUDA_SPACE_ATOMIC_MASK+1));
  return ptr;
}

void init_lock_array_cuda_space() {
  int is_initialized = 0;
  if(! is_initialized) {
    int* lock_array_ptr = lock_array_cuda_space_ptr();
    cudaMemcpyToSymbol( kokkos_impl_cuda_atomic_lock_array , & lock_array_ptr , sizeof(int*) );
    init_lock_array_kernel<<<(CUDA_SPACE_ATOMIC_MASK+255)/256,256>>>();
  }
}

}
}
#endif // KOKKOS_HAVE_CUDA

