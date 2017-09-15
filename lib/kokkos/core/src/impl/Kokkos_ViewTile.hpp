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

  enum { SHIFT_0 = Kokkos::Impl::integral_power_of_two(Layout::N0) };
  enum { SHIFT_1 = Kokkos::Impl::integral_power_of_two(Layout::N1) };
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

  KOKKOS_INLINE_FUNCTION constexpr
  array_layout layout() const
    { return array_layout( m_dim.N0 , m_dim.N1 ); }

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

  KOKKOS_FUNCTION_DEFAULTED ~ViewOffset() = default ;
  KOKKOS_INLINE_FUNCTION ViewOffset() = default ;
  KOKKOS_INLINE_FUNCTION ViewOffset( const ViewOffset & ) = default ;
  KOKKOS_INLINE_FUNCTION ViewOffset & operator = ( const ViewOffset & ) = default ;

  template< unsigned TrivialScalarSize >
  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset( std::integral_constant<unsigned,TrivialScalarSize> const & ,
                        array_layout const arg_layout )
    : m_dim( arg_layout.dimension[0], arg_layout.dimension[1], 0, 0, 0, 0, 0, 0 )
    , m_tile_N0( ( arg_layout.dimension[0] + MASK_0 ) >> SHIFT_0 /* number of tiles in first dimension */ )
    {}
};

template< typename T , unsigned N0 , unsigned N1 , class ... P
        , typename iType0 , typename iType1
        >
struct ViewMapping
  < void
  , Kokkos::ViewTraits<T**,Kokkos::LayoutTileLeft<N0,N1,true>,P...>
  , Kokkos::LayoutTileLeft<N0,N1,true>
  , iType0
  , iType1 >
{
  typedef Kokkos::LayoutTileLeft<N0,N1,true>  src_layout ;
  typedef Kokkos::ViewTraits< T** , src_layout , P... > src_traits ;
  typedef Kokkos::ViewTraits< T[N0][N1] , LayoutLeft , P ... > traits ;
  typedef Kokkos::View< T[N0][N1] , LayoutLeft , P ... > type ;

  KOKKOS_INLINE_FUNCTION static
  void assign( ViewMapping< traits , void > & dst
             , const ViewMapping< src_traits , void > & src
             , const src_layout &
             , const size_t i_tile0
             , const size_t i_tile1
             )
    {
      typedef ViewMapping< traits , void >        dst_map_type ;
      typedef ViewMapping< src_traits , void >    src_map_type ;
      typedef typename dst_map_type::handle_type  dst_handle_type ;
      typedef typename dst_map_type::offset_type  dst_offset_type ;
      typedef typename src_map_type::offset_type  src_offset_type ;

      dst = dst_map_type(
         dst_handle_type( src.m_handle +
                        ( ( i_tile0 + src.m_offset.m_tile_N0 * i_tile1 ) << src_offset_type::SHIFT_T ) ) ,
         dst_offset_type() );
    }
};

} /* namespace Impl */
} /* namespace Kokkos */

namespace Kokkos {

template< typename T , unsigned N0 , unsigned N1 , class ... P >
KOKKOS_INLINE_FUNCTION
Kokkos::View< T[N0][N1] , LayoutLeft , P... >
tile_subview( const Kokkos::View<T**,Kokkos::LayoutTileLeft<N0,N1,true>,P...> & src
            , const size_t i_tile0
            , const size_t i_tile1
            )
{
  // Force the specialized ViewMapping for extracting a tile
  // by using the first subview argument as the layout.
  typedef Kokkos::LayoutTileLeft<N0,N1,true> SrcLayout ;

  return Kokkos::View< T[N0][N1] , LayoutLeft , P... >
    ( src , SrcLayout() , i_tile0 , i_tile1 );
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXPERIENTAL_VIEWTILE_HPP */

