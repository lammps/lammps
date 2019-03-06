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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXPERIMENTAL_CUDA_VIEW_HPP
#define KOKKOS_EXPERIMENTAL_CUDA_VIEW_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_CUDA )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Cuda Texture fetches can be performed for 4, 8 and 16 byte objects (int,int2,int4)
// Via reinterpret_case this can be used to support all scalar types of those sizes.
// Any other scalar type falls back to either normal reads out of global memory,
// or using the __ldg intrinsic on Kepler GPUs or newer (Compute Capability >= 3.0)

template< typename ValueType , typename AliasType >
struct CudaTextureFetch {

  ::cudaTextureObject_t   m_obj ;
  const ValueType       * m_ptr ;
  int                     m_offset ;

  // Deference operator pulls through texture object and returns by value
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  ValueType operator[]( const iType & i ) const
    {
#if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
      AliasType v = tex1Dfetch<AliasType>( m_obj , i + m_offset );
      return  *(reinterpret_cast<ValueType*> (&v));
#else
      return m_ptr[ i ];
#endif
    }

  // Pointer to referenced memory
  KOKKOS_INLINE_FUNCTION
  operator const ValueType * () const { return m_ptr ; }


  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : m_obj() , m_ptr() , m_offset() {}

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs )
    : m_obj(     rhs.m_obj )
    , m_ptr(     rhs.m_ptr )
    , m_offset(  rhs.m_offset )
    {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( CudaTextureFetch && rhs )
    : m_obj(     rhs.m_obj )
    , m_ptr(     rhs.m_ptr )
    , m_offset(  rhs.m_offset )
    {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
    {
      m_obj     = rhs.m_obj ;
      m_ptr     = rhs.m_ptr ;
      m_offset  = rhs.m_offset ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( CudaTextureFetch && rhs )
    {
      m_obj     = rhs.m_obj ;
      m_ptr     = rhs.m_ptr ;
      m_offset  = rhs.m_offset ;
      return *this ;
    }

  // Texture object spans the entire allocation.
  // This handle may view a subset of the allocation, so an offset is required.
  template< class CudaMemorySpace >
  inline explicit
  CudaTextureFetch( const ValueType * const arg_ptr
                  , Kokkos::Impl::SharedAllocationRecord< CudaMemorySpace , void > * record
                  )
    : m_obj( record->template attach_texture_object< AliasType >() )
    , m_ptr( arg_ptr )
    , m_offset( record->attach_texture_object_offset( reinterpret_cast<const AliasType*>( arg_ptr ) ) )
    {}

  // Texture object spans the entire allocation.
  // This handle may view a subset of the allocation, so an offset is required.
  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs , size_t offset )
    : m_obj(     rhs.m_obj )
    , m_ptr(     rhs.m_ptr + offset)
    , m_offset( offset + rhs.m_offset )
    {}
};

#if defined( KOKKOS_ENABLE_CUDA_LDG_INTRINSIC )

template< typename ValueType , typename AliasType >
struct CudaLDGFetch {

  const ValueType * m_ptr ;

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  ValueType operator[]( const iType & i ) const
    {
      #ifdef __CUDA_ARCH__
      AliasType v = __ldg(reinterpret_cast<const AliasType*>(&m_ptr[i]));
      return  *(reinterpret_cast<ValueType*> (&v));
      #else
      return m_ptr[i];
      #endif
    }

  KOKKOS_INLINE_FUNCTION
  operator const ValueType * () const { return m_ptr ; }

  KOKKOS_INLINE_FUNCTION
  CudaLDGFetch() : m_ptr() {}

  KOKKOS_INLINE_FUNCTION
  ~CudaLDGFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaLDGFetch( const CudaLDGFetch & rhs )
    : m_ptr( rhs.m_ptr )
    {}

  KOKKOS_INLINE_FUNCTION
  CudaLDGFetch( CudaLDGFetch && rhs )
    : m_ptr( rhs.m_ptr )
    {}

  KOKKOS_INLINE_FUNCTION
  CudaLDGFetch & operator = ( const CudaLDGFetch & rhs )
    {
      m_ptr = rhs.m_ptr ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION
  CudaLDGFetch & operator = ( CudaLDGFetch && rhs )
    {
      m_ptr = rhs.m_ptr ;
      return *this ;
    }

  template< class CudaMemorySpace >
  inline explicit
  CudaLDGFetch( const ValueType * const arg_ptr
              , Kokkos::Impl::SharedAllocationRecord<CudaMemorySpace,void>*
              )
    : m_ptr( arg_ptr )
    {}

  KOKKOS_INLINE_FUNCTION
  CudaLDGFetch( CudaLDGFetch const rhs ,size_t offset)
    : m_ptr( rhs.m_ptr + offset )
    {}

};

#endif

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Replace Default ViewDataHandle with Cuda texture fetch specialization
 *          if 'const' value type, CudaSpace and random access.
 */
template< class Traits >
class ViewDataHandle< Traits ,
  typename std::enable_if<(
    // Is Cuda memory space
    ( std::is_same< typename Traits::memory_space,Kokkos::CudaSpace>::value ||
      std::is_same< typename Traits::memory_space,Kokkos::CudaUVMSpace>::value )
    &&
    // Is a trivial const value of 4, 8, or 16 bytes
    std::is_trivial<typename Traits::const_value_type>::value
    &&
    std::is_same<typename Traits::const_value_type,typename Traits::value_type>::value
    &&
    ( sizeof(typename Traits::const_value_type) ==  4 ||
      sizeof(typename Traits::const_value_type) ==  8 ||
      sizeof(typename Traits::const_value_type) == 16 )
    &&
    // Random access trait
    ( Traits::memory_traits::RandomAccess != 0 )
  )>::type >
{
public:

  using track_type  = Kokkos::Impl::SharedAllocationTracker ;

  using value_type  = typename Traits::const_value_type ;
  using return_type = typename Traits::const_value_type ; // NOT a reference

  using alias_type = typename std::conditional< ( sizeof(value_type) ==  4 ) , int ,
                     typename std::conditional< ( sizeof(value_type) ==  8 ) , ::int2 ,
                     typename std::conditional< ( sizeof(value_type) == 16 ) , ::int4 , void
                     >::type
                     >::type
                     >::type ;

#if defined( KOKKOS_ENABLE_CUDA_LDG_INTRINSIC )
  using handle_type = Kokkos::Impl::CudaLDGFetch< value_type , alias_type > ;
#else
  using handle_type = Kokkos::Impl::CudaTextureFetch< value_type , alias_type > ;
#endif

  KOKKOS_INLINE_FUNCTION
  static handle_type const & assign( handle_type const & arg_handle , track_type const & /* arg_tracker */ )
    {
      return arg_handle ;
    }

  KOKKOS_INLINE_FUNCTION
  static handle_type const assign( handle_type const & arg_handle , size_t offset )
    {
      return handle_type(arg_handle,offset) ;
    }

  KOKKOS_INLINE_FUNCTION
  static handle_type assign( value_type * arg_data_ptr, track_type const & arg_tracker )
    {
      if(arg_data_ptr == NULL) return handle_type();

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      // Assignment of texture = non-texture requires creation of a texture object
      // which can only occur on the host.  In addition, 'get_record' is only valid
      // if called in a host execution space


      typedef typename Traits::memory_space memory_space ;
      typedef typename Impl::SharedAllocationRecord<memory_space,void> record ;

      record * const r = arg_tracker.template get_record< memory_space >();

#if ! defined( KOKKOS_ENABLE_CUDA_LDG_INTRINSIC )
      if ( 0 == r ) {
        Kokkos::abort("Cuda const random access View using Cuda texture memory requires Kokkos to allocate the View's memory");
      }
#endif

      return handle_type( arg_data_ptr , r );

#else
      Kokkos::Impl::cuda_abort("Cannot create Cuda texture object from within a Cuda kernel");
      return handle_type();
#endif
    }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_CUDA ) */
#endif /* #ifndef KOKKOS_CUDA_VIEW_HPP */

