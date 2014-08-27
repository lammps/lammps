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

#if defined( __CUDACC__ )
#include <cuda_runtime.h>
#endif

#include <Kokkos_View.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_CudaSpace.hpp>
#include <Kokkos_CudaTypes.hpp>
#include <Cuda/Kokkos_Cuda_abort.hpp>
#include <Kokkos_Atomic.hpp>

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
      Kokkos::cuda_abort("Kokkos::View array bounds violation");
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

cuda_texture_object_type
cuda_texture_object_attach(
  const cudaChannelFormatDesc & ,
  const void * const );

int cuda_texture_object_release(
    cuda_texture_object_type obj
    );

int cuda_texture_object_release(
    const void * const
    );

template< typename TextureType >
inline
cuda_texture_object_type
cuda_texture_object_attach( const void * const base_view_ptr )
{
  return cuda_texture_object_attach( cudaCreateChannelDesc<TextureType>() , base_view_ptr );
}

#else

typedef const void * cuda_texture_object_type ;

template< typename TextureType >
inline
cuda_texture_object_type
cuda_texture_object_attach( const void * const )
{ return 0 ; }

int cuda_texture_object_release(
    const void * const
    );

#endif

//----------------------------------------------------------------------------

// Cuda Texture fetches can be performed for 4, 8 and 16 byte objects (int,int2,int4)
// Via reinterpret_case this can be used to support all scalar types of those sizes.
// Any other scalar type falls back to either normal reads out of global memory,
// or using the __ldg intrinsic on Kepler GPUs or newer (Compute Capability >= 3.0)

template< typename T, size_t size = sizeof(T) >
struct alias_type {
  typedef void type;
};

template< typename T >
struct alias_type<T,4> {
  typedef int type;
};

template< typename T >
struct alias_type<T,8> {
  typedef int2 type;
};

template< typename T >
struct alias_type<T,16> {
  typedef int4 type;
};

template< typename ValueType, typename AliasType = typename alias_type<ValueType>::type >
struct CudaTextureFetch {
  private:

    cuda_texture_object_type  obj ;

    int* ref_count ;
  public:

    const ValueType * ptr ;

    KOKKOS_INLINE_FUNCTION
    CudaTextureFetch() : obj( 0 ) , ref_count(0), ptr( 0 ) {}

    KOKKOS_INLINE_FUNCTION
    ~CudaTextureFetch() {
#ifndef __CUDA_ARCH__
      if(ptr!=NULL) {
        //printf("Release D: %p %p %i\n",this,ptr,ref_count[0]);
        int count = Kokkos::atomic_fetch_add(ref_count,-1);
        if(count == 1) {
          cuda_texture_object_release(obj);
          delete [] ref_count;
        }
      }
#endif
    }

    KOKKOS_INLINE_FUNCTION
    CudaTextureFetch( const CudaTextureFetch & rhs ) {
#ifndef __CUDA_ARCH__
      if(rhs.ptr != NULL) {
        obj = rhs.obj;
        ptr = rhs.ptr;
        ref_count = rhs.ref_count;
        Kokkos::atomic_fetch_add(ref_count,1);
        //printf("Attach C: %p %p %i\n",this,ptr,ref_count[0]);
      } else {
        obj = 0;
        ref_count = NULL;
        ptr = NULL;
      }
#else
      obj = rhs.obj;
      ref_count = rhs.ref_count;
      ptr = rhs.ptr;
#endif
}

    KOKKOS_INLINE_FUNCTION
    CudaTextureFetch & operator = ( const CudaTextureFetch & rhs ) {
#ifndef __CUDA_ARCH__
      if(ptr!=NULL) {
        //printf("Release A: %p %p %i\n",this,ptr,ref_count[0]);
        int count = Kokkos::atomic_fetch_add(ref_count,-1);
        if(count == 1) {
          cuda_texture_object_release(obj);
          delete [] ref_count;
        }
      }
      if(rhs.ptr!=NULL) {
        obj = rhs.obj;
        ptr = rhs.ptr;
        ref_count = rhs.ref_count;
        Kokkos::atomic_fetch_add(ref_count,1);
        //printf("Attach A: %p %p %i\n",this,ptr,ref_count[0]);
      } else {
        obj = 0;
        ref_count = NULL;
        ptr = NULL;
      }
#else
      obj = rhs.obj;
      ref_count = rhs.ref_count;
      ptr = rhs.ptr;
#endif
      return *this ;
    }

    explicit KOKKOS_INLINE_FUNCTION
    CudaTextureFetch( ValueType * const base_view_ptr ) {
#ifndef __CUDA_ARCH__
      if( base_view_ptr != NULL ) {
        obj = cuda_texture_object_attach<AliasType>( base_view_ptr );
        ref_count = new int[1];
        ref_count[0] = 1;
        ptr = base_view_ptr;
        //printf("Attach PC: %p %p %i\n",this,ptr,ref_count[0]);
      } else {
        obj = 0;
        ref_count = NULL;
        ptr = NULL;
      }
#else
      cuda_abort("ERROR: Trying to assign a non texture_fetch view to a texture_fetch view in a Device kernel\n.");
#endif
    }

    KOKKOS_INLINE_FUNCTION
    CudaTextureFetch & operator = (const ValueType* base_view_ptr) {
#ifndef __CUDA_ARCH__
      if(ptr!=NULL) {
        //printf("Release P: %p %p %i\n",this,ptr,ref_count[0]);
        int count = Kokkos::atomic_fetch_add(ref_count,-1);
        if(count == 1) {
          cuda_texture_object_release(obj);
          delete [] ref_count;
        }
      }
      if( base_view_ptr != NULL ) {
        obj = cuda_texture_object_attach<AliasType>( base_view_ptr );
        ref_count = new int[1];
        ref_count[0] = 1;
        ptr = base_view_ptr;
        //printf("Attach P: %p %p %i\n",this,ptr,ref_count[0]);
      } else {
        obj = 0;
        ref_count = NULL;
        ptr = NULL;
      }
#else
      cuda_abort("ERROR: Trying to assign a non texture_fetch view to a texture_fetch view in a Device kernel\n.");
#endif
      return *this;
    }

    template< typename iType >
    KOKKOS_INLINE_FUNCTION
    ValueType operator[]( const iType & i ) const
    {
  #if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
  // Enable the usage of the _ldg intrinsic even in cases where texture fetches work
  // Currently texture fetches are faster, but that might change in the future
  #ifdef KOKKOS_USE_LDG_INTRINSIC
      return _ldg(&ptr[i]);
  #else
      AliasType v = tex1Dfetch<AliasType>( obj , i );

      return  *(reinterpret_cast<ValueType*> (&v));
  #endif
  #else
      return ptr[ i ];
  #endif
    }

};

template< typename ValueType >
struct CudaTextureFetch< const ValueType, void > {

  const ValueType * ptr ;

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : ptr(0) {};

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {
  }

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs ) : ptr(rhs.ptr) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs ) {
    ptr = rhs.ptr;
    return *this ;
  }

  explicit KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( ValueType * const base_view_ptr ) {
    ptr = base_view_ptr;
  }

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = (const ValueType* base_view_ptr) {
    ptr = base_view_ptr;
    return *this;
  }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  ValueType operator[]( const iType & i ) const
  {
  #if defined( __CUDA_ARCH__ ) && ( 300 <= __CUDA_ARCH__ )
    return _ldg(&ptr[i]);
  #else
    return ptr[ i ];
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
template<class ViewTraits>
class ViewDataHandle<ViewTraits,
typename enable_if<is_same<typename ViewTraits::memory_space,CudaSpace>::value &&
                   is_same<typename ViewTraits::const_value_type,typename ViewTraits::value_type>::value &&
                   ViewTraits::memory_traits::RandomAccess
                  >::type> {
  typedef ViewDataHandle self_type;
public:
  enum {ReferenceAble = 0};
  typedef Impl::CudaTextureFetch<typename ViewTraits::value_type> type;
  typedef typename ViewTraits::value_type return_type;

  static type allocate(std::string label, size_t count) {
    return type((typename ViewTraits::value_type*)
                typename ViewTraits::memory_space::allocate( label ,
                typeid(typename ViewTraits::value_type) ,
                sizeof(typename ViewTraits::value_type) ,
                count ));
  }

  KOKKOS_INLINE_FUNCTION
  static typename ViewTraits::value_type* get_raw_ptr(type handle) {
    return handle.ptr;
  }
};

}
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CUDA_VIEW_HPP */

