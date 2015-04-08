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

#ifndef KOKKOS_CUDA_VIEW_HPP
#define KOKKOS_CUDA_VIEW_HPP

#include <cstring>

#include <Kokkos_HostSpace.hpp>
#include <Kokkos_CudaSpace.hpp>
#include <Kokkos_CudaTypes.hpp>
#include <Kokkos_View.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
struct AssertShapeBoundsAbort< CudaSpace >
{
  KOKKOS_INLINE_FUNCTION
  static void apply( const size_t /* rank */ ,
                     const size_t /* n0 */ , const size_t /* n1 */ ,
                     const size_t /* n2 */ , const size_t /* n3 */ ,
                     const size_t /* n4 */ , const size_t /* n5 */ ,
                     const size_t /* n6 */ , const size_t /* n7 */ ,

                     const size_t /* arg_rank */ ,
                     const size_t /* i0 */ , const size_t /* i1 */ ,
                     const size_t /* i2 */ , const size_t /* i3 */ ,
                     const size_t /* i4 */ , const size_t /* i5 */ ,
                     const size_t /* i6 */ , const size_t /* i7 */ )
    {
      Kokkos::abort("Kokkos::View array bounds violation");
    }
};

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Cuda 5.0 <texture_types.h> defines 'cudaTextureObject_t'
// to be an 'unsigned long long'.  This chould change with
// future version of Cuda and this typedef would have to
// change accordingly.

#if defined( CUDA_VERSION ) && ( 5000 <= CUDA_VERSION )

typedef enable_if<
  sizeof(::cudaTextureObject_t) == sizeof(const void *) ,
  ::cudaTextureObject_t >::type cuda_texture_object_type ;

#else

typedef const void * cuda_texture_object_type ;

#endif

//----------------------------------------------------------------------------
// Cuda Texture fetches can be performed for 4, 8 and 16 byte objects (int,int2,int4)
// Via reinterpret_case this can be used to support all scalar types of those sizes.
// Any other scalar type falls back to either normal reads out of global memory,
// or using the __ldg intrinsic on Kepler GPUs or newer (Compute Capability >= 3.0)

template< typename ValueType
        , class MemorySpace
        , class AliasType =
            typename Kokkos::Impl::if_c< ( sizeof(ValueType) ==  4 ) , int ,
            typename Kokkos::Impl::if_c< ( sizeof(ValueType) ==  8 ) , int2 ,
            typename Kokkos::Impl::if_c< ( sizeof(ValueType) == 16 ) , int4 , void
            >::type
            >::type
            >::type
        >
class CudaTextureFetch {
private:

  cuda_texture_object_type  m_obj ;
  const ValueType         * m_alloc_ptr ;
  int                       m_offset ;

public:

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : m_obj( 0 ) , m_alloc_ptr(0) , m_offset(0) {}

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs )
    : m_obj(       rhs.m_obj )
    , m_alloc_ptr( rhs.m_alloc_ptr )
    , m_offset(    rhs.m_offset )
    {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
    {
      m_obj       = rhs.m_obj ;
      m_alloc_ptr = rhs.m_alloc_ptr ;
      m_offset    = rhs.m_offset ;
      return *this ;
    }


  KOKKOS_INLINE_FUNCTION explicit
  CudaTextureFetch( const ValueType * const arg_ptr )
    : m_obj( 0 ) , m_alloc_ptr(0) , m_offset(0)
    {
#if defined( __CUDACC__ ) && ! defined( __CUDA_ARCH__ )
      MemorySpace::texture_object_attach( arg_ptr
                                        , sizeof(ValueType)
                                        , cudaCreateChannelDesc< AliasType >()
                                        , & m_obj
                                        , reinterpret_cast<const void **>( & m_alloc_ptr )
                                        , & m_offset
                                        );
#endif
    }

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const ValueType * arg_ptr )
    {
#if defined( __CUDACC__ ) && ! defined( __CUDA_ARCH__ )
      MemorySpace::texture_object_attach( arg_ptr
                                        , sizeof(ValueType)
                                        , cudaCreateChannelDesc< AliasType >()
                                        , & m_obj
                                        , reinterpret_cast<const void **>( & m_alloc_ptr )
                                        , & m_offset
                                        );
#endif
      return *this ;
    }


  KOKKOS_INLINE_FUNCTION
  operator const ValueType * () const { return m_alloc_ptr + m_offset ; }


  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  ValueType operator[]( const iType & i ) const
    {
#if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
#if defined( KOKKOS_USE_LDG_INTRINSIC )
      // Enable the usage of the _ldg intrinsic even in cases where texture fetches work
      // Currently texture fetches are faster, but that might change in the future
      return _ldg( & m_alloc_ptr[i+m_offset] );
#else /* ! defined( KOKKOS_USE_LDG_INTRINSIC ) */
      AliasType v = tex1Dfetch<AliasType>( m_obj , i + m_offset );

      return  *(reinterpret_cast<ValueType*> (&v));
#endif /* ! defined( KOKKOS_USE_LDG_INTRINSIC ) */
#else  /* ! defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ ) */
      return m_alloc_ptr[ i + m_offset ];
#endif
  }
};

template< typename ValueType >
class CudaTextureFetch< const ValueType, void >
{
private:
  const ValueType * m_ptr ;
public:

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : m_ptr(0) {};

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {
  }

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs ) : m_ptr(rhs.m_ptr) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs ) {
    m_ptr = rhs.m_ptr;
    return *this ;
  }

  explicit KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( ValueType * const base_view_ptr ) {
    m_ptr = base_view_ptr;
  }

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = (const ValueType* base_view_ptr) {
    m_ptr = base_view_ptr;
    return *this;
  }


  KOKKOS_INLINE_FUNCTION
  operator const ValueType * () const { return m_ptr ; }


  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  ValueType operator[]( const iType & i ) const
  {
  #if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
    return _ldg(&m_ptr[i]);
  #else
    return m_ptr[ i ];
  #endif
  }
};

} // namespace Impl
} // namespace Kokkos


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  Replace Default ViewDataHandle with Cuda texture fetch specialization
 *          if 'const' value type, CudaSpace and random access.
 */
template< class ViewTraits >
class ViewDataHandle< ViewTraits ,
  typename enable_if< ( is_same< typename ViewTraits::memory_space,CudaSpace>::value ||
                        is_same< typename ViewTraits::memory_space,CudaUVMSpace>::value )
                      &&
                      is_same<typename ViewTraits::const_value_type,typename ViewTraits::value_type>::value
                      &&
                      ViewTraits::memory_traits::RandomAccess
                    >::type >
{
public:
  enum { ReturnTypeIsReference = false };

  typedef Impl::CudaTextureFetch< typename ViewTraits::value_type
                                , typename ViewTraits::memory_space > handle_type;

  typedef typename ViewTraits::value_type return_type;
};

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CUDA_VIEW_HPP */

