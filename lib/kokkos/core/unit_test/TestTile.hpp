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

#ifndef TEST_TILE_HPP
#define TEST_TILE_HPP

//========================================================================
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_ViewTile.hpp>

namespace TestTile {

template < typename Device, typename TileLayout >
struct ReduceTileErrors
{
  typedef Device execution_space;
  typedef Kokkos::View< ptrdiff_t**, TileLayout, Device >  array_type;
  typedef Kokkos::View< ptrdiff_t[ TileLayout::N0 ][ TileLayout::N1 ], Kokkos::LayoutLeft, Device >  tile_type;
  typedef ptrdiff_t value_type;

  array_type m_array;

  ReduceTileErrors( array_type a ) : m_array( a ) {}

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & errors ) { errors = 0; }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & errors,
                    const volatile value_type & src_errors )
  {
    errors += src_errors;
  }

  // Initialize.
  KOKKOS_INLINE_FUNCTION
  void operator()( size_t iwork ) const
  {
    const size_t i = iwork % m_array.extent(0);
    const size_t j = iwork / m_array.extent(0);

    if ( j < m_array.extent(1) ) {
      m_array( i, j ) = &m_array( i, j ) - &m_array( 0, 0 );

      //printf( "m_array(%d, %d) = %d\n", int( i ), int( j ), int( m_array( i, j ) ) );
    }
  }

  // Verify:
  KOKKOS_INLINE_FUNCTION
  void operator()( size_t iwork, value_type & errors ) const
  {
    const size_t tile_dim0 = ( m_array.extent(0) + TileLayout::N0 - 1 ) / TileLayout::N0;
    const size_t tile_dim1 = ( m_array.extent(1) + TileLayout::N1 - 1 ) / TileLayout::N1;

    const size_t itile = iwork % tile_dim0;
    const size_t jtile = iwork / tile_dim0;

    if ( jtile < tile_dim1 ) {
      tile_type tile = Kokkos::tile_subview( m_array, itile, jtile );

      if ( tile( 0, 0 ) != ptrdiff_t( ( itile + jtile * tile_dim0 ) * TileLayout::N0 * TileLayout::N1 ) ) {
        ++errors;
      }
      else {
        for ( size_t j = 0; j < size_t( TileLayout::N1 ); ++j ) {
          for ( size_t i = 0; i < size_t( TileLayout::N0 ); ++i ) {
            const size_t iglobal = i + itile * TileLayout::N0;
            const size_t jglobal = j + jtile * TileLayout::N1;

            if ( iglobal < m_array.extent(0) && jglobal < m_array.extent(1) ) {
              if ( tile( i, j ) != ptrdiff_t( tile( 0, 0 ) + i + j * TileLayout::N0 ) ) ++errors;

              //printf( "tile(%d, %d)(%d, %d) = %d\n", int( itile ), int( jtile ), int( i ), int( j ), int( tile( i, j ) ) );
            }
          }
        }
      }
    }
  }
};

template< class Space, unsigned N0, unsigned N1 >
void test( const size_t dim0, const size_t dim1 )
{
  typedef Kokkos::LayoutTileLeft< N0, N1 >  array_layout;
  typedef ReduceTileErrors< Space, array_layout > functor_type;

  const size_t tile_dim0 = ( dim0 + N0 - 1 ) / N0;
  const size_t tile_dim1 = ( dim1 + N1 - 1 ) / N1;

  typename functor_type::array_type array( "", dim0, dim1 );

  Kokkos::parallel_for( Kokkos::RangePolicy< Space, size_t >( 0, dim0 * dim1 ), functor_type( array ) );

  ptrdiff_t error = 0;

  Kokkos::parallel_reduce( Kokkos::RangePolicy< Space, size_t >( 0, tile_dim0 * tile_dim1 ), functor_type( array ), error );

  EXPECT_EQ( error, ptrdiff_t( 0 ) );
}

} // namespace TestTile

namespace Test {
TEST_F( TEST_CATEGORY, tile_layout )
{
  TestTile::test< TEST_EXECSPACE, 1, 1 >( 1, 1 );
  TestTile::test< TEST_EXECSPACE, 1, 1 >( 2, 3 );
  TestTile::test< TEST_EXECSPACE, 1, 1 >( 9, 10 );

  TestTile::test< TEST_EXECSPACE, 2, 2 >( 1, 1 );
  TestTile::test< TEST_EXECSPACE, 2, 2 >( 2, 3 );
  TestTile::test< TEST_EXECSPACE, 2, 2 >( 4, 4 );
  TestTile::test< TEST_EXECSPACE, 2, 2 >( 9, 9 );

  TestTile::test< TEST_EXECSPACE, 2, 4 >( 9, 9 );
  TestTile::test< TEST_EXECSPACE, 4, 2 >( 9, 9 );

  TestTile::test< TEST_EXECSPACE, 4, 4 >( 1, 1 );
  TestTile::test< TEST_EXECSPACE, 4, 4 >( 4, 4 );
  TestTile::test< TEST_EXECSPACE, 4, 4 >( 9, 9 );
  TestTile::test< TEST_EXECSPACE, 4, 4 >( 9, 11 );

  TestTile::test< TEST_EXECSPACE, 8, 8 >( 1, 1 );
  TestTile::test< TEST_EXECSPACE, 8, 8 >( 4, 4 );
  TestTile::test< TEST_EXECSPACE, 8, 8 >( 9, 9 );
  TestTile::test< TEST_EXECSPACE, 8, 8 >( 9, 11 );
}

}

#endif // KOKKOS_ENABLE_DEPRECATED_CODE
//=====================================================================

#endif //TEST_TILE_HPP
