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

#ifndef KOKKOS_VIEWTILELEFT_HPP
#define KOKKOS_VIEWTILELEFT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class T , unsigned N0 , unsigned N1 , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< T , void , LayoutTileLeft<N0,N1> , MemorySpace , MemoryTraits >
{
  typedef ViewDefault type ;
};

struct ViewTile {};

template< class ShapeType , unsigned N0 , unsigned N1 >
struct ViewOffset< ShapeType
                 , LayoutTileLeft<N0,N1,true> /* Only accept properly shaped tiles */
                 , typename Impl::enable_if<( 2 == ShapeType::rank
                                              &&
                                              2 == ShapeType::rank_dynamic
                                            )>::type >
  : public ShapeType
{
  enum { SHIFT_0 = Impl::power_of_two<N0>::value };
  enum { SHIFT_1 = Impl::power_of_two<N1>::value };
  enum { MASK_0  = N0 - 1 };
  enum { MASK_1  = N1 - 1 };

  typedef size_t                      size_type ;
  typedef ShapeType                   shape_type ;
  typedef LayoutTileLeft<N0,N1,true>  array_layout ;

  enum { has_padding = true };

  size_type tile_N0 ;

  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset & rhs )
    {
      shape_type::N0 = rhs.N0 ;
      shape_type::N1 = rhs.N1 ;
      tile_N0 = ( rhs.N0 + MASK_0 ) >> SHIFT_0 ; // number of tiles in first dimension
    }

  KOKKOS_INLINE_FUNCTION
  void assign( size_t n0 , size_t n1
             , int = 0 , int = 0
             , int = 0 , int = 0
             , int = 0 , int = 0
             , int = 0
             )
    {
      shape_type::N0 = n0 ;
      shape_type::N1 = n1 ;
      tile_N0 = ( n0 + MASK_0 ) >> SHIFT_0 ; // number of tiles in first dimension
    }


  KOKKOS_INLINE_FUNCTION
  void set_padding() {}


  template< typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION
  size_type operator()( I0 const & i0 , I1 const & i1
                      , int = 0 , int = 0
                      , int = 0 , int = 0
                      , int = 0 , int = 0
                      ) const
    {
      return /* ( ( Tile offset                             ) *  ( Tile size       ) ) */
                ( ( (i0>>SHIFT_0) + tile_N0 * (i1>>SHIFT_1) ) << (SHIFT_0 + SHIFT_1) ) +
             /* ( Offset within tile                       ) */
                ( (i0 & MASK_0) + ((i1 & MASK_1)<<SHIFT_0) ) ;
    }

  template< typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION
  size_type tile_begin( I0 const & i_tile0 , I1 const & i_tile1 ) const
    {
      return ( i_tile0 + tile_N0 * i_tile1 ) << ( SHIFT_0 + SHIFT_1 );
    }


  KOKKOS_INLINE_FUNCTION
  size_type capacity() const
    {
      // ( TileDim0 * ( TileDim1 ) ) * TileSize
      return ( tile_N0 * ( ( shape_type::N1 + MASK_1 ) >> SHIFT_1 ) ) << ( SHIFT_0 + SHIFT_1 );
    }
};

template<>
struct ViewAssignment< ViewTile , void , void >
{
  // Some compilers have type-matching issues on the integer values when using:
  //   template< class T , unsigned N0 , unsigned N1 , class A2 , class A3 >
  template< class T , unsigned dN0 , unsigned dN1
          , class A2 , class A3
          , unsigned sN0 , unsigned sN1 >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( View< T[dN0][dN1], LayoutLeft, A2, A3, Impl::ViewDefault > & dst
                , View< T** , LayoutTileLeft<sN0,sN1,true>, A2, A3, Impl::ViewDefault > const & src
                , size_t const i_tile0
                , typename Impl::enable_if< unsigned(dN0) == unsigned(sN0) &&
                                            unsigned(dN1) == unsigned(sN1)
                                          , size_t const
                                          >::type i_tile1
                )
   {
     // Destination is always contiguous but source may be non-contiguous
     // so don't assign the whole view management object.
     // Just query and appropriately set the reference-count state.

     if ( ! src.m_management.is_managed() ) dst.m_management.set_unmanaged();

     dst.m_ptr_on_device = src.m_ptr_on_device + src.m_offset_map.tile_begin(i_tile0,i_tile1);

     dst.m_management.increment( dst.m_ptr_on_device );
   }
};

} /* namespace Impl */
} /* namespace Kokkos */

namespace Kokkos {

template< class T , unsigned N0, unsigned N1, class A2, class A3 >
KOKKOS_INLINE_FUNCTION
View< T[N0][N1], LayoutLeft, A2, A3, Impl::ViewDefault >
tile_subview( const View<T**,LayoutTileLeft<N0,N1,true>,A2,A3,Impl::ViewDefault> & src
            , const size_t i_tile0
            , const size_t i_tile1
            )
{
  View< T[N0][N1], LayoutLeft, A2, A3, Impl::ViewDefault > dst ;

  (void) Impl::ViewAssignment< Impl::ViewTile , void , void >( dst , src , i_tile0 , i_tile1 );

  return dst ;
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEWTILELEFT_HPP */

