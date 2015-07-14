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

#ifndef KOKKOS_CUDA_VIEW_HPP
#define KOKKOS_CUDA_VIEW_HPP

#include <Kokkos_Macros.hpp>

/* only compile this file if CUDA is enabled for Kokkos */
#ifdef KOKKOS_HAVE_CUDA

#include <cstring>

#include <Kokkos_HostSpace.hpp>
#include <Kokkos_CudaSpace.hpp>
#include <Kokkos_View.hpp>

#include <Cuda/Kokkos_Cuda_BasicAllocators.hpp>

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

//----------------------------------------------------------------------------
// Cuda Texture fetches can be performed for 4, 8 and 16 byte objects (int,int2,int4)
// Via reinterpret_case this can be used to support all scalar types of those sizes.
// Any other scalar type falls back to either normal reads out of global memory,
// or using the __ldg intrinsic on Kepler GPUs or newer (Compute Capability >= 3.0)

template< typename ValueType
        , class MemorySpace
        , class AliasType =
            typename Kokkos::Impl::if_c< ( sizeof(ValueType) ==  4 ) , int ,
            typename Kokkos::Impl::if_c< ( sizeof(ValueType) ==  8 ) , ::int2 ,
            typename Kokkos::Impl::if_c< ( sizeof(ValueType) == 16 ) , ::int4 , void
            >::type
            >::type
            >::type
        >
class CudaTextureFetch {
private:

  cuda_texture_object_type  m_obj ;
  const ValueType         * m_alloc_ptr ;
  int                       m_offset ;

  void attach( const ValueType * const arg_ptr, AllocationTracker const & tracker )
  {
    typedef char const * const byte;

    m_alloc_ptr = reinterpret_cast<ValueType *>(tracker.alloc_ptr());

    size_t byte_offset = reinterpret_cast<byte>(arg_ptr) - reinterpret_cast<byte>(m_alloc_ptr);
    const bool ok_aligned = 0 == byte_offset % sizeof(ValueType);

    const size_t count = tracker.alloc_size() / sizeof(ValueType);
    const bool ok_contains = (m_alloc_ptr <= arg_ptr) && (arg_ptr < (m_alloc_ptr + count));

    if (ok_aligned && ok_contains) {
      if (tracker.attribute() == NULL ) {
        MemorySpace::texture_object_attach(
            tracker
            , sizeof(ValueType)
            , cudaCreateChannelDesc< AliasType >()
            );
      }
      m_obj = dynamic_cast<TextureAttribute*>(tracker.attribute())->m_tex_obj;
      m_offset = arg_ptr - m_alloc_ptr;
    }
    else if( !ok_contains ) {
      throw_runtime_exception("Error: cannot attach a texture object to a tracker which does not bound the pointer.");
    }
    else {
      throw_runtime_exception("Error: cannot attach a texture object to an incorrectly aligned pointer.");
    }
  }

public:

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : m_obj() , m_alloc_ptr() , m_offset() {}

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
  CudaTextureFetch( const ValueType * const arg_ptr, AllocationTracker const & tracker )
    : m_obj( 0 ) , m_alloc_ptr(0) , m_offset(0)
    {
      #if defined( KOKKOS_USE_LDG_INTRINSIC )
        m_alloc_ptr(arg_ptr);
      #elif defined( __CUDACC__ ) && ! defined( __CUDA_ARCH__ )
        if ( arg_ptr != NULL ) {
          if ( tracker.is_valid() ) {
            attach( arg_ptr, tracker );
          }
          else {
            AllocationTracker found_tracker = AllocationTracker::find<typename MemorySpace::allocator>(arg_ptr);
            if ( found_tracker.is_valid() ) {
              attach( arg_ptr, found_tracker );
            } else {
              throw_runtime_exception("Error: cannot attach a texture object to an untracked pointer!");
            }
          }
        }
      #endif
    }

  KOKKOS_INLINE_FUNCTION
  operator const ValueType * () const { return m_alloc_ptr + m_offset ; }


  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  ValueType operator[]( const iType & i ) const
    {
      #if defined( KOKKOS_USE_LDG_INTRINSIC ) && defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
        AliasType v = __ldg(reinterpret_cast<AliasType*>(&m_alloc_ptr[i]));
        return  *(reinterpret_cast<ValueType*> (&v));
      #elif defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
        AliasType v = tex1Dfetch<AliasType>( m_obj , i + m_offset );
        return  *(reinterpret_cast<ValueType*> (&v));
      #else
        return m_alloc_ptr[ i + m_offset ];
      #endif
  }
};

template< typename ValueType, class MemorySpace >
class CudaTextureFetch< const ValueType, MemorySpace, void >
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
  CudaTextureFetch( const ValueType * ptr, const AllocationTracker & ) : m_ptr(ptr) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs ) : m_ptr(rhs.m_ptr) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs ) {
    m_ptr = rhs.m_ptr;
    return *this ;
  }

  explicit KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( ValueType * const base_view_ptr, AllocationTracker const & /*tracker*/ ) {
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
    return m_ptr[ i ];
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
                                , typename ViewTraits::memory_space> handle_type;

  KOKKOS_INLINE_FUNCTION
  static handle_type create_handle( typename ViewTraits::value_type * arg_data_ptr, AllocationTracker const & arg_tracker )
  {
    return handle_type(arg_data_ptr, arg_tracker);
  }

  typedef typename ViewTraits::value_type return_type;
};

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif // KOKKOS_HAVE_CUDA
#endif /* #ifndef KOKKOS_CUDA_VIEW_HPP */

