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

  public:

    const ValueType * ptr ;

    KOKKOS_INLINE_FUNCTION
    CudaTextureFetch() : obj( 0 ) , ptr( 0 ) {}

    KOKKOS_INLINE_FUNCTION
    ~CudaTextureFetch() {}

    KOKKOS_INLINE_FUNCTION
    CudaTextureFetch( const CudaTextureFetch & rhs )
      : obj( rhs.obj ) , ptr( rhs.ptr ) {}

    KOKKOS_INLINE_FUNCTION
    CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
      { obj = rhs.obj ; ptr = rhs.ptr ; return *this ; }

    explicit
    CudaTextureFetch( const ValueType * const base_view_ptr )
      : obj( cuda_texture_object_attach<AliasType>( base_view_ptr ) )
      , ptr( base_view_ptr ) {}

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
private:

  cuda_texture_object_type  obj ;

public:

  const ValueType * ptr ;

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch() : obj( 0 ) , ptr( 0 ) {}

  KOKKOS_INLINE_FUNCTION
  ~CudaTextureFetch() {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch( const CudaTextureFetch & rhs )
    : obj( rhs.obj ) , ptr( rhs.ptr ) {}

  KOKKOS_INLINE_FUNCTION
  CudaTextureFetch & operator = ( const CudaTextureFetch & rhs )
    { obj = rhs.obj ; ptr = rhs.ptr ; return *this ; }

  explicit
  CudaTextureFetch( ValueType * const base_view_ptr )
    : obj( cuda_texture_object_attach<ValueType>( base_view_ptr ) )
    , ptr( base_view_ptr ) {}

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

struct ViewCudaTexture {};

#if defined( CUDA_VERSION ) && ( 5000 <= CUDA_VERSION )

/** \brief  Replace ViewDefault specialization with Cuda texture fetch specialization
 *          if 'const' value type and random access.
 */
template< class ValueType , class MemoryTraits >
struct ViewSpecialize< const ValueType , void , LayoutLeft , CudaSpace , MemoryTraits >
{
  typedef typename if_c< MemoryTraits::RandomAccess , ViewCudaTexture , ViewDefault >::type type ;
};

template< class ValueType , class MemoryTraits >
struct ViewSpecialize< const ValueType , void , LayoutRight , CudaSpace , MemoryTraits >
{
  typedef typename if_c< MemoryTraits::RandomAccess , ViewCudaTexture , ViewDefault >::type type ;
};

#endif

//----------------------------------------------------------------------------

template<>
struct ViewAssignment< ViewCudaTexture , ViewCudaTexture , void >
{
  /** \brief Assign compatible views */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,ViewCudaTexture> & dst ,
                  const View<ST,SL,SD,SM,ViewCudaTexture> & src ,
                  const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> , ViewTraits<ST,SL,SD,SM> >::value
                  ) >::type * = 0 )
  {
    dst.m_tracking.decrement( dst.m_texture.ptr );

    dst.m_texture  = src.m_texture ;
    dst.m_offset_map.assign( src.m_offset_map );
    dst.m_tracking = src.m_tracking ;

    dst.m_tracking.increment( dst.m_texture.ptr );
  }
};


template<>
struct ViewAssignment< ViewCudaTexture , ViewDefault , void >
{
  /** \brief Assign compatible views */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  inline
  ViewAssignment(       View<DT,DL,DD,DM,ViewCudaTexture> & dst ,
                  const View<ST,SL,SD,SM,ViewDefault> & src ,
                  const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                  )>::type * = 0 )
  {
    dst.m_tracking.decrement( dst.m_texture.ptr );

    dst.m_texture = CudaTextureFetch< typename ViewTraits<DT,DL,DD,DM>::value_type >( src.m_ptr_on_device );

    dst.m_offset_map.assign( src.m_offset_map );

    dst.m_tracking  = src.m_tracking ;

    dst.m_tracking.increment( dst.m_texture.ptr );
  }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
template< class T , class L, class D , class M >
class View< T , L , D , M , Impl::ViewCudaTexture >
  : public ViewTraits< T , L , D , M >
{
public:

  typedef ViewTraits< T , L , D , M > traits ;

private:

  template< class , class , class > friend struct Impl::ViewAssignment ;

  typedef Impl::ViewOffset< typename traits::shape_type
                          , typename traits::array_layout
                          > offset_map_type ;

  Impl::CudaTextureFetch<typename traits::value_type > m_texture ;
  offset_map_type                                      m_offset_map ;
  Impl::ViewTracking< traits >                         m_tracking ;

public:

  typedef Impl::ViewCudaTexture specialize ;

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type::host_mirror_device_type ,
                void > HostMirror ;

  enum { Rank = traits::rank };

  KOKKOS_INLINE_FUNCTION typename traits::shape_type shape() const { return m_offset_map ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_offset_map.N0 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_offset_map.N1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_2() const { return m_offset_map.N2 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_3() const { return m_offset_map.N3 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_4() const { return m_offset_map.N4 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_5() const { return m_offset_map.N5 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_6() const { return m_offset_map.N6 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_7() const { return m_offset_map.N7 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type size() const
  {
    return   m_offset_map.N0
           * m_offset_map.N1
           * m_offset_map.N2
           * m_offset_map.N3
           * m_offset_map.N4
           * m_offset_map.N5
           * m_offset_map.N6
           * m_offset_map.N7
           ;
  }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type dimension( const iType & i ) const
    { return Impl::dimension( m_offset_map , i ); }

  //------------------------------------

  View() : m_texture()
   { m_offset_map.assign(0,0,0,0,0,0,0,0); }

  KOKKOS_INLINE_FUNCTION
  ~View() { m_tracking.decrement( m_texture.ptr ); }

  View( const View & rhs )
    : m_texture( rhs.m_texture )
    {
      m_offset_map.assign( rhs.m_offset_map );
      m_tracking = rhs.m_tracking ;
      m_tracking.increment( m_texture.ptr );
    }

  View & operator = ( const View & rhs )
    {
      (void)Impl::ViewAssignment< Impl::ViewCudaTexture , Impl::ViewCudaTexture >( *this , rhs );
      return *this ;
    }

  template< class RT , class RL, class RD , class RM , class RS >
  View( const View<RT,RL,RD,RM,RS> & rhs )
    : m_texture(0)
    {
      Impl::ViewAssignment< Impl::ViewCudaTexture , RS >( *this , rhs );
    }

  template< class RT , class RL, class RD, class RM , class RS >
  View & operator = ( const View<RT,RL,RD,RM,RS> & rhs )
    {
      Impl::ViewAssignment< Impl::ViewCudaTexture , RS >( *this , rhs );
      return *this ;
    }

  template< typename TT >
  explicit inline
  View( TT * ptr ,
        const size_t n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        typename Impl::enable_if<(
          Impl::is_same<TT,typename traits::value_type>::value
        ), const size_t >::type n7 = 0 )
    : m_texture( Impl::CudaTextureFetch< typename traits::value_type >(ptr))
    {
      m_offset_map.assign( n0, n1, n2, n3, n4, n5, n6, n7 );
      m_tracking = false ;
    }

  //------------------------------------

  KOKKOS_FORCEINLINE_FUNCTION
  bool is_null() const { return 0 == m_texture.ptr ; }

  //------------------------------------
  // Rank = 1 access operators:

  template < typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type , traits, typename traits::array_layout, 1 , iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_offset_map, i0 );
      return m_texture[ i0 ];
    }

  template < typename iType0 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type , traits , typename traits::array_layout, 1 , iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );
      KOKKOS_ASSERT_SHAPE_BOUNDS_1( m_offset_map, i0 );
      return m_texture[ i0 ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type , traits, typename traits::array_layout, 2, iType0, iType1 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_offset_map, i0,i1 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ m_offset_map(i0,i1) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type ,
                                      traits, typename traits::array_layout, 3, iType0, iType1, iType2 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_3( m_offset_map, i0,i1,i2 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ m_offset_map(i0,i1,i2) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type ,
                                      traits, typename traits::array_layout, 4, iType0, iType1, iType2, iType3 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_4( m_offset_map, i0,i1,i2,i3 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ m_offset_map(i0,i1,i2,i3) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type ,
                                      traits, typename traits::array_layout, 5, iType0, iType1, iType2, iType3, iType4 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_5( m_offset_map, i0,i1,i2,i3,i4 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ m_offset_map(i0,i1,i2,i3,i4) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type ,
                                      traits, typename traits::array_layout, 6, iType0, iType1, iType2, iType3, iType4, iType5 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_6( m_offset_map, i0,i1,i2,i3,i4,i5 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ m_offset_map(i0,i1,i2,i3,i4,i5) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type ,
                                      traits, typename traits::array_layout, 7, iType0, iType1, iType2, iType3, iType4, iType5, iType6 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_7( m_offset_map, i0,i1,i2,i3,i4,i5,i6 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ m_offset_map(i0,i1,i2,i3,i4,i5,i6) ];
    }

  template< typename iType0 , typename iType1 , typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 , typename iType6 , typename iType7 >
  KOKKOS_FORCEINLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< typename traits::value_type ,
                                      traits, typename traits::array_layout, 8, iType0, iType1, iType2, iType3, iType4, iType5, iType6, iType7 >::type
    operator() ( const iType0 & i0 , const iType1 & i1 , const iType2 & i2 , const iType3 & i3 ,
                 const iType4 & i4 , const iType5 & i5 , const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_8( m_offset_map, i0,i1,i2,i3,i4,i5,i6,i7 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_texture.ptr );

      return m_texture[ m_offset_map(i0,i1,i2,i3,i4,i5,i6,i7) ];
    }

  //------------------------------------

  KOKKOS_FORCEINLINE_FUNCTION
  typename traits::value_type * ptr_on_device() const { return m_texture.ptr ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const { m_offset_map.stride(s); }
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_CUDA_VIEW_HPP */

