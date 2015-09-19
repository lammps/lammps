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

#ifndef KOKKOS_EXPERIMENTAL_VIEWTILE_HPP
#define KOKKOS_EXPERIMENTAL_VIEWTILE_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

// View mapping for rank two tiled array

template< class L >
struct is_layout_tile : public std::false_type {};

template< unsigned N0 , unsigned N1 >
struct is_layout_tile< Kokkos::LayoutTileLeft<N0,N1,true> > : public std::true_type {};

template< class Dimension , class Layout >
struct ViewOffset< Dimension , Layout ,
  typename std::enable_if<(
    ( Dimension::rank == 2 )
    &&
    is_layout_tile< Layout >::value
  )>::type >
{
public:

  enum { SHIFT_0 = Kokkos::Impl::power_of_two<Layout::N0>::value };
  enum { SHIFT_1 = Kokkos::Impl::power_of_two<Layout::N1>::value };
  enum { SHIFT_T = SHIFT_0 + SHIFT_1 };
  enum { MASK_0  = Layout::N0 - 1 };
  enum { MASK_1  = Layout::N1 - 1 };

  // Is an irregular layout that does not have uniform striding for each index.
  using is_mapping_plugin = std::true_type ;
  using is_regular        = std::false_type ;

  typedef size_t     size_type ;
  typedef Dimension  dimension_type ;
  typedef Layout     array_layout ;

  dimension_type m_dim ;
  size_type      m_tile_N0 ;

  //----------------------------------------

  // Only instantiated for rank 2
  template< typename I0 , typename I1 >
  KOKKOS_INLINE_FUNCTION constexpr
  size_type operator()( I0 const & i0 , I1 const & i1
                      , int = 0 , int = 0
                      , int = 0 , int = 0
                      , int = 0 , int = 0
                      ) const
    {
      return /* ( ( Tile offset                               ) * Tile size ) */
                ( ( (i0>>SHIFT_0) + m_tile_N0 * (i1>>SHIFT_1) ) << SHIFT_T) +
             /* ( Offset within tile                       ) */
                ( (i0 & MASK_0) + ((i1 & MASK_1)<<SHIFT_0) ) ;
    }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const { return m_dim.N0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const { return m_dim.N1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION constexpr size_type size() const { return m_dim.N0 * m_dim.N1 ; }

  // Strides are meaningless due to irregularity
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return 0 ; }

  KOKKOS_INLINE_FUNCTION constexpr size_type span() const
    {
      // ( TileDim0 * ( TileDim1 ) ) * TileSize
      return ( m_tile_N0 * ( ( m_dim.N1 + MASK_1 ) >> SHIFT_1 ) ) << SHIFT_T ;
    }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const
    {
      // Only if dimensions align with tile size
      return ( m_dim.N0 & MASK_0 ) == 0 && ( m_dim.N1 & MASK_1 ) == 0 ;
    }

  //----------------------------------------

  ~ViewOffset() = default ;
  ViewOffset() = default ;
  ViewOffset( const ViewOffset & ) = default ;
  ViewOffset & operator = ( const ViewOffset & ) = default ;

  template< unsigned TrivialScalarSize >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( std::integral_constant<unsigned,TrivialScalarSize> const &
                      , size_t aN0   , size_t aN1
                      , unsigned , unsigned , unsigned , unsigned , unsigned , unsigned )
    : m_dim( aN0, aN1, 0, 0, 0, 0, 0, 0 )
    , m_tile_N0( ( aN0 + MASK_0 ) >> SHIFT_0 /* number of tiles in first dimension */ )
    {}
};

} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

namespace Kokkos {
namespace Experimental {

// Using View with an invalid data type to construct the tiling subview.
// View is a friend of View so we use this invalid data type partial specialization
// to access implementation of both source and destination view for constructing
// the tile subview.

template< unsigned N0 , unsigned N1 >
struct View< void , Kokkos::LayoutTileLeft<N0,N1,true> , void , void >
{
  typedef Kokkos::LayoutTileLeft<N0,N1,true>  Layout ;

  template< typename T , class A2 , class A3 >
  KOKKOS_INLINE_FUNCTION static
  Kokkos::Experimental::View< T[N0][N1] , LayoutLeft , A2 , A3 >
  tile_subview( const Kokkos::Experimental::View<T**,Layout,A2,A3> & src
              , const size_t i_tile0
              , const size_t i_tile1
              )
    {
      typedef Kokkos::Experimental::View<T**,Layout,A2,A3>                   SrcView ;
      typedef Kokkos::Experimental::View< T[N0][N1] , LayoutLeft , A2 , A3 > DstView ;

      typedef typename SrcView::map_type::offset_type  src_offset_type ;
      typedef typename DstView::map_type               dst_map_type ;
      typedef typename DstView::map_type::handle_type  dst_handle_type ;
      typedef typename DstView::map_type::offset_type  dst_offset_type ;

      return DstView( src.m_track ,
                      dst_map_type(
                        dst_handle_type( src.m_map.m_handle +
                                         ( ( i_tile0 + src.m_map.m_offset.m_tile_N0 * i_tile1 ) << src_offset_type::SHIFT_T ) ) ,
                        dst_offset_type() )
                    );
    }
};

template< typename T , unsigned N0 , unsigned N1 , class A2 , class A3 >
KOKKOS_INLINE_FUNCTION
Kokkos::Experimental::View< T[N0][N1] , LayoutLeft , A2 , A3 >
tile_subview( const Kokkos::Experimental::View<T**,Kokkos::LayoutTileLeft<N0,N1,true>,A2,A3> & src
            , const size_t i_tile0
            , const size_t i_tile1
            )
{
  return View< void , Kokkos::LayoutTileLeft<N0,N1,true> , void , void >::
    tile_subview( src , i_tile0 , i_tile1 );
}

} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXPERIENTAL_VIEWTILE_HPP */

