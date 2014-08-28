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

#ifndef KOKKOS_SHAPE_HPP
#define KOKKOS_SHAPE_HPP

#include <typeinfo>
#include <utility>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Layout.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  The shape of a Kokkos with dynamic and static dimensions.
 *          Dynamic dimensions are member values and static dimensions are
 *          'static const' values.
 *
 *  The upper bound on the array rank is eight.
 */
template< unsigned ScalarSize ,
          unsigned Rank ,
          unsigned s0  = 1 ,
          unsigned s1  = 1 ,
          unsigned s2  = 1 ,
          unsigned s3  = 1 ,
          unsigned s4  = 1 ,
          unsigned s5  = 1 ,
          unsigned s6  = 1 ,
          unsigned s7  = 1 >
struct Shape ;

template< class ShapeType , class Layout >
struct ShapeMap ;

//----------------------------------------------------------------------------
/** \brief  Shape equality if the value type, layout, and dimensions
 *          are equal.
 */
template< unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          unsigned ySize , unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
KOKKOS_INLINE_FUNCTION
bool operator == ( const Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
                   const Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{
  enum { same_size = xSize == ySize };
  enum { same_rank = xRank == yRank };

  return same_size && same_rank &&
         size_t( x.N0 )   == size_t( y.N0 ) &&
         unsigned( x.N1 ) == unsigned( y.N1 ) &&
         unsigned( x.N2 ) == unsigned( y.N2 ) &&
         unsigned( x.N3 ) == unsigned( y.N3 ) &&
         unsigned( x.N4 ) == unsigned( y.N4 ) &&
         unsigned( x.N5 ) == unsigned( y.N5 ) &&
         unsigned( x.N6 ) == unsigned( y.N6 ) &&
         unsigned( x.N7 ) == unsigned( y.N7 ) ;
}

template< unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          unsigned ySize ,unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
KOKKOS_INLINE_FUNCTION
bool operator != ( const Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
                   const Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{ return ! operator == ( x , y ); }

//----------------------------------------------------------------------------

void assert_counts_are_equal_throw(
  const size_t x_count ,
  const size_t y_count );

inline
void assert_counts_are_equal(
  const size_t x_count ,
  const size_t y_count )
{
  if ( x_count != y_count ) {
    assert_counts_are_equal_throw( x_count , y_count );
  }
}

void assert_shapes_are_equal_throw(
  const unsigned x_scalar_size ,
  const unsigned x_rank ,
  const size_t   x_N0 , const unsigned x_N1 ,
  const unsigned x_N2 , const unsigned x_N3 ,
  const unsigned x_N4 , const unsigned x_N5 ,
  const unsigned x_N6 , const unsigned x_N7 ,

  const unsigned y_scalar_size ,
  const unsigned y_rank ,
  const size_t   y_N0 , const unsigned y_N1 ,
  const unsigned y_N2 , const unsigned y_N3 ,
  const unsigned y_N4 , const unsigned y_N5 ,
  const unsigned y_N6 , const unsigned y_N7 );

template< unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          unsigned ySize , unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
inline
void assert_shapes_are_equal(
  const Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
  const Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{
  typedef Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> x_type ;
  typedef Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> y_type ;

  if ( x != y ) {
    assert_shapes_are_equal_throw(
      x_type::scalar_size, x_type::rank, x.N0, x.N1, x.N2, x.N3, x.N4, x.N5, x.N6, x.N7,
      y_type::scalar_size, y_type::rank, y.N0, y.N1, y.N2, y.N3, y.N4, y.N5, y.N6, y.N7 );
  }
}

template< unsigned xSize , unsigned xRank ,
          unsigned xN0 , unsigned xN1 , unsigned xN2 , unsigned xN3 ,
          unsigned xN4 , unsigned xN5 , unsigned xN6 , unsigned xN7 ,

          unsigned ySize , unsigned yRank ,
          unsigned yN0 , unsigned yN1 , unsigned yN2 , unsigned yN3 ,
          unsigned yN4 , unsigned yN5 , unsigned yN6 , unsigned yN7 >
void assert_shapes_equal_dimension(
  const Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> & x ,
  const Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> & y )
{
  typedef Shape<xSize,xRank,xN0,xN1,xN2,xN3,xN4,xN5,xN6,xN7> x_type ;
  typedef Shape<ySize,yRank,yN0,yN1,yN2,yN3,yN4,yN5,yN6,yN7> y_type ;

  // Omit comparison of scalar_size.
  if ( unsigned( x.rank ) != unsigned( y.rank ) ||
       size_t( x.N0 )   != size_t( y.N0 ) || 
       unsigned( x.N1 ) != unsigned( y.N1 ) || 
       unsigned( x.N2 ) != unsigned( y.N2 ) || 
       unsigned( x.N3 ) != unsigned( y.N3 ) ||
       unsigned( x.N4 ) != unsigned( y.N4 ) || 
       unsigned( x.N5 ) != unsigned( y.N5 ) || 
       unsigned( x.N6 ) != unsigned( y.N6 ) || 
       unsigned( x.N7 ) != unsigned( y.N7 ) ) {
    assert_shapes_are_equal_throw(
      x_type::scalar_size, x_type::rank, x.N0, x.N1, x.N2, x.N3, x.N4, x.N5, x.N6, x.N7,
      y_type::scalar_size, y_type::rank, y.N0, y.N1, y.N2, y.N3, y.N4, y.N5, y.N6, y.N7 );
  }
}

//----------------------------------------------------------------------------

template< class ShapeType > struct assert_shape_is_rank_zero ;
template< class ShapeType > struct assert_shape_is_rank_one ;

template< unsigned Size >
struct assert_shape_is_rank_zero< Shape<Size,0> >
  : public true_type {};

template< unsigned Size , unsigned s0 >
struct assert_shape_is_rank_one< Shape<Size,1,s0> >
  : public true_type {};

//----------------------------------------------------------------------------

/** \brief  Array bounds assertion templated on the execution space
 *          to allow device-specific abort code.
 */
template< class Space >
struct AssertShapeBoundsAbort ;

template<>
struct AssertShapeBoundsAbort< Kokkos::HostSpace >
{
  static void apply( const size_t rank ,
                     const size_t n0 , const size_t n1 ,
                     const size_t n2 , const size_t n3 ,
                     const size_t n4 , const size_t n5 ,
                     const size_t n6 , const size_t n7 ,
                     const size_t arg_rank ,
                     const size_t i0 , const size_t i1 ,
                     const size_t i2 , const size_t i3 ,
                     const size_t i4 , const size_t i5 ,
                     const size_t i6 , const size_t i7 );
};

template< class ExecutionDevice >
struct AssertShapeBoundsAbort
{
  KOKKOS_INLINE_FUNCTION
  static void apply( const size_t rank ,
                     const size_t n0 , const size_t n1 ,
                     const size_t n2 , const size_t n3 ,
                     const size_t n4 , const size_t n5 ,
                     const size_t n6 , const size_t n7 ,
                     const size_t arg_rank ,
                     const size_t i0 , const size_t i1 ,
                     const size_t i2 , const size_t i3 ,
                     const size_t i4 , const size_t i5 ,
                     const size_t i6 , const size_t i7 )
    {
      AssertShapeBoundsAbort< Kokkos::HostSpace >
        ::apply( rank ,    n0 , n1 , n2 , n3 , n4 , n5 , n6 , n7 ,
                 arg_rank, i0 , i1 , i2 , i3 , i4 , i5 , i6 , i7 );
    }
};

template< class ShapeType >
KOKKOS_INLINE_FUNCTION
void assert_shape_bounds( const ShapeType & shape ,
                          const size_t arg_rank ,
                          const size_t i0 ,
                          const size_t i1 = 0 ,
                          const size_t i2 = 0 ,
                          const size_t i3 = 0 ,
                          const size_t i4 = 0 ,
                          const size_t i5 = 0 ,
                          const size_t i6 = 0 ,
                          const size_t i7 = 0 )
{
  // Must supply at least as many indices as ranks.
  // Every index must be within bounds.
  const bool ok = ShapeType::rank <= arg_rank &&
                  i0 < shape.N0 && 
                  i1 < shape.N1 &&
                  i2 < shape.N2 &&
                  i3 < shape.N3 &&
                  i4 < shape.N4 &&
                  i5 < shape.N5 &&
                  i6 < shape.N6 &&
                  i7 < shape.N7 ;

  if ( ! ok ) {
    AssertShapeBoundsAbort< Kokkos::Impl::ActiveExecutionMemorySpace >
      ::apply( ShapeType::rank ,
               shape.N0 , shape.N1 , shape.N2 , shape.N3 ,
               shape.N4 , shape.N5 , shape.N6 , shape.N7 ,
               arg_rank , i0 , i1 , i2 , i3 , i4 , i5 , i6 , i7 );
  }
}

#if defined( KOKKOS_EXPRESSION_CHECK )
#define KOKKOS_ASSERT_SHAPE_BOUNDS_1( S , I0 ) assert_shape_bounds(S,1,I0);
#define KOKKOS_ASSERT_SHAPE_BOUNDS_2( S , I0 , I1 ) assert_shape_bounds(S,2,I0,I1);
#define KOKKOS_ASSERT_SHAPE_BOUNDS_3( S , I0 , I1 , I2 ) assert_shape_bounds(S,3,I0,I1,I2);
#define KOKKOS_ASSERT_SHAPE_BOUNDS_4( S , I0 , I1 , I2 , I3 ) assert_shape_bounds(S,4,I0,I1,I2,I3);
#define KOKKOS_ASSERT_SHAPE_BOUNDS_5( S , I0 , I1 , I2 , I3 , I4 ) assert_shape_bounds(S,5,I0,I1,I2,I3,I4);
#define KOKKOS_ASSERT_SHAPE_BOUNDS_6( S , I0 , I1 , I2 , I3 , I4 , I5 ) assert_shape_bounds(S,6,I0,I1,I2,I3,I4,I5);
#define KOKKOS_ASSERT_SHAPE_BOUNDS_7( S , I0 , I1 , I2 , I3 , I4 , I5 , I6 ) assert_shape_bounds(S,7,I0,I1,I2,I3,I4,I5,I6);
#define KOKKOS_ASSERT_SHAPE_BOUNDS_8( S , I0 , I1 , I2 , I3 , I4 , I5 , I6 , I7 ) assert_shape_bounds(S,8,I0,I1,I2,I3,I4,I5,I6,I7);
#else
#define KOKKOS_ASSERT_SHAPE_BOUNDS_1( S , I0 ) /* */
#define KOKKOS_ASSERT_SHAPE_BOUNDS_2( S , I0 , I1 ) /* */
#define KOKKOS_ASSERT_SHAPE_BOUNDS_3( S , I0 , I1 , I2 ) /* */
#define KOKKOS_ASSERT_SHAPE_BOUNDS_4( S , I0 , I1 , I2 , I3 ) /* */
#define KOKKOS_ASSERT_SHAPE_BOUNDS_5( S , I0 , I1 , I2 , I3 , I4 ) /* */
#define KOKKOS_ASSERT_SHAPE_BOUNDS_6( S , I0 , I1 , I2 , I3 , I4 , I5 ) /* */
#define KOKKOS_ASSERT_SHAPE_BOUNDS_7( S , I0 , I1 , I2 , I3 , I4 , I5 , I6 ) /* */
#define KOKKOS_ASSERT_SHAPE_BOUNDS_8( S , I0 , I1 , I2 , I3 , I4 , I5 , I6 , I7 ) /* */
#endif


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Specialization and optimization for the Rank 0 shape.

template < unsigned ScalarSize >
struct Shape< ScalarSize , 0, 1,1,1,1, 1,1,1,1 >
{
  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 0 };
  enum { rank         = 0 };

  enum { N0 = 1 };
  enum { N1 = 1 };
  enum { N2 = 1 };
  enum { N3 = 1 };
  enum { N4 = 1 };
  enum { N5 = 1 };
  enum { N6 = 1 };
  enum { N7 = 1 };

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  {}
};

//----------------------------------------------------------------------------
// All-static dimension array

template < unsigned ScalarSize ,
           unsigned Rank ,
           unsigned s0 ,
           unsigned s1 ,
           unsigned s2 ,
           unsigned s3 ,
           unsigned s4 ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape {

  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 0 };
  enum { rank         = Rank };

  enum { N0 = s0 };
  enum { N1 = s1 };
  enum { N2 = s2 };
  enum { N3 = s3 };
  enum { N4 = s4 };
  enum { N5 = s5 };
  enum { N6 = s6 };
  enum { N7 = s7 };

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  {}
};

// 1 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize ,
           unsigned Rank ,
           unsigned s1 ,
           unsigned s2 ,
           unsigned s3 ,
           unsigned s4 ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,s1,s2,s3, s4,s5,s6,s7 >
{
  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 1 };
  enum { rank         = Rank };

  size_t N0 ; // For 1 == dynamic_rank allow  N0 > 2^32

  enum { N1 = s1 };
  enum { N2 = s2 };
  enum { N3 = s3 };
  enum { N4 = s4 };
  enum { N5 = s5 };
  enum { N6 = s6 };
  enum { N7 = s7 };

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               size_t n0 , unsigned = 0 , unsigned = 0 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; }
};

// 2 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s2 ,
           unsigned s3 ,
           unsigned s4 ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,s2,s3, s4,s5,s6,s7 >
{
  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 2 };
  enum { rank         = Rank };

  unsigned N0 ;
  unsigned N1 ;

  enum { N2 = s2 };
  enum { N3 = s3 };
  enum { N4 = s4 };
  enum { N5 = s5 };
  enum { N6 = s6 };
  enum { N7 = s7 };

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned = 0 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; s.N1 = n1 ; }
};

// 3 == dynamic_rank <= rank <= 8
template < unsigned Rank , unsigned ScalarSize ,
           unsigned s3 ,
           unsigned s4 ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,0,s3, s4,s5,s6,s7>
{
  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 3 };
  enum { rank         = Rank };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;

  enum { N3 = s3 };
  enum { N4 = s4 };
  enum { N5 = s5 };
  enum { N6 = s6 };
  enum { N7 = s7 };

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned = 0 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; }
};

// 4 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s4 ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank, 0,0,0,0, s4,s5,s6,s7 >
{
  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 4 };
  enum { rank         = Rank };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;

  enum { N4 = s4 };
  enum { N5 = s5 };
  enum { N6 = s6 };
  enum { N7 = s7 };

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned = 0 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ; }
};

// 5 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s5 ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,0,0, 0,s5,s6,s7 >
{
  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 5 };
  enum { rank         = Rank };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;

  enum { N5 = s5 };
  enum { N6 = s6 };
  enum { N7 = s7 };

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned n4 , unsigned = 0 , unsigned = 0 , unsigned = 0 )
  { s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ; s.N4 = n4 ; }
};

// 6 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s6 ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,0,0, 0,0,s6,s7 >
{
  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 6 };
  enum { rank         = Rank };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;

  enum { N6 = s6 };
  enum { N7 = s7 };

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned n4 , unsigned n5 = 0 , unsigned = 0 , unsigned = 0 )
  {
    s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ;
    s.N4 = n4 ; s.N5 = n5 ;
  }
};

// 7 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize , unsigned Rank ,
           unsigned s7 >
struct Shape< ScalarSize , Rank , 0,0,0,0, 0,0,0,s7 >
{
  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 7 };
  enum { rank         = Rank };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;

  enum { N7 = s7 };

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned n4 , unsigned n5 , unsigned n6 , unsigned = 0 )
  {
    s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ;
    s.N4 = n4 ; s.N5 = n5 ; s.N6 = n6 ;
  }
};

// 8 == dynamic_rank <= rank <= 8
template < unsigned ScalarSize >
struct Shape< ScalarSize , 8 , 0,0,0,0, 0,0,0,0 >
{
  enum { scalar_size   = ScalarSize };
  enum { rank_dynamic = 8 };
  enum { rank         = 8 };

  unsigned N0 ;
  unsigned N1 ;
  unsigned N2 ;
  unsigned N3 ;
  unsigned N4 ;
  unsigned N5 ;
  unsigned N6 ;
  unsigned N7 ;

  KOKKOS_INLINE_FUNCTION
  static
  void assign( Shape & s ,
               unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3 ,
               unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7 )
  {
    s.N0 = n0 ; s.N1 = n1 ; s.N2 = n2 ; s.N3 = n3 ;
    s.N4 = n4 ; s.N5 = n5 ; s.N6 = n6 ; s.N7 = n7 ;
  }
};

//----------------------------------------------------------------------------

template< class ShapeType , unsigned N ,
          unsigned R = ShapeType::rank_dynamic >
struct ShapeInsert ;

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 0 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 N ,
                 ShapeType::N0 ,
                 ShapeType::N1 ,
                 ShapeType::N2 ,
                 ShapeType::N3 ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 1 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 N ,
                 ShapeType::N1 ,
                 ShapeType::N2 ,
                 ShapeType::N3 ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 2 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N2 ,
                 ShapeType::N3 ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 3 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N3 ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 4 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N4 ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 5 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N5 ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 6 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 N ,
                 ShapeType::N6 > type ;
};

template< class ShapeType , unsigned N >
struct ShapeInsert< ShapeType , N , 7 >
{
  typedef Shape< ShapeType::scalar_size ,
                 ShapeType::rank + 1 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 0 ,
                 N > type ;
};

//----------------------------------------------------------------------------

template< class DstShape , class SrcShape ,
          unsigned DstRankDynamic   = DstShape::rank_dynamic ,
          bool     DstRankDynamicOK = unsigned(DstShape::rank_dynamic) >= unsigned(SrcShape::rank_dynamic) >
struct ShapeCompatible { enum { value = false }; };

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 8 , true >
{
  enum { value = unsigned(DstShape::scalar_size) == unsigned(SrcShape::scalar_size) };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 7 , true >
{
  enum { value = unsigned(DstShape::scalar_size) == unsigned(SrcShape::scalar_size) &&
                 unsigned(DstShape::N7) == unsigned(SrcShape::N7) };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 6 , true >
{
  enum { value = unsigned(DstShape::scalar_size) == unsigned(SrcShape::scalar_size) &&
                 unsigned(DstShape::N6) == unsigned(SrcShape::N6) &&
                 unsigned(DstShape::N7) == unsigned(SrcShape::N7) };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 5 , true >
{
  enum { value = unsigned(DstShape::scalar_size) == unsigned(SrcShape::scalar_size) &&
                 unsigned(DstShape::N5) == unsigned(SrcShape::N5) &&
                 unsigned(DstShape::N6) == unsigned(SrcShape::N6) &&
                 unsigned(DstShape::N7) == unsigned(SrcShape::N7) };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 4 , true >
{
  enum { value = unsigned(DstShape::scalar_size) == unsigned(SrcShape::scalar_size) &&
                 unsigned(DstShape::N4) == unsigned(SrcShape::N4) &&
                 unsigned(DstShape::N5) == unsigned(SrcShape::N5) &&
                 unsigned(DstShape::N6) == unsigned(SrcShape::N6) &&
                 unsigned(DstShape::N7) == unsigned(SrcShape::N7) };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 3 , true >
{
  enum { value = unsigned(DstShape::scalar_size) == unsigned(SrcShape::scalar_size) &&
                 unsigned(DstShape::N3) == unsigned(SrcShape::N3) &&
                 unsigned(DstShape::N4) == unsigned(SrcShape::N4) &&
                 unsigned(DstShape::N5) == unsigned(SrcShape::N5) &&
                 unsigned(DstShape::N6) == unsigned(SrcShape::N6) &&
                 unsigned(DstShape::N7) == unsigned(SrcShape::N7) };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 2 , true >
{
  enum { value = unsigned(DstShape::scalar_size) == unsigned(SrcShape::scalar_size) &&
                 unsigned(DstShape::N2) == unsigned(SrcShape::N2) &&
                 unsigned(DstShape::N3) == unsigned(SrcShape::N3) &&
                 unsigned(DstShape::N4) == unsigned(SrcShape::N4) &&
                 unsigned(DstShape::N5) == unsigned(SrcShape::N5) &&
                 unsigned(DstShape::N6) == unsigned(SrcShape::N6) &&
                 unsigned(DstShape::N7) == unsigned(SrcShape::N7) };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 1 , true >
{
  enum { value = unsigned(DstShape::scalar_size) == unsigned(SrcShape::scalar_size) &&
                 unsigned(DstShape::N1) == unsigned(SrcShape::N1) &&
                 unsigned(DstShape::N2) == unsigned(SrcShape::N2) &&
                 unsigned(DstShape::N3) == unsigned(SrcShape::N3) &&
                 unsigned(DstShape::N4) == unsigned(SrcShape::N4) &&
                 unsigned(DstShape::N5) == unsigned(SrcShape::N5) &&
                 unsigned(DstShape::N6) == unsigned(SrcShape::N6) &&
                 unsigned(DstShape::N7) == unsigned(SrcShape::N7) };
};

template< class DstShape , class SrcShape >
struct ShapeCompatible< DstShape , SrcShape , 0 , true >
{
  enum { value = unsigned(DstShape::scalar_size) == unsigned(SrcShape::scalar_size) &&
                 unsigned(DstShape::N0) == unsigned(SrcShape::N0) &&
                 unsigned(DstShape::N1) == unsigned(SrcShape::N1) &&
                 unsigned(DstShape::N2) == unsigned(SrcShape::N2) &&
                 unsigned(DstShape::N3) == unsigned(SrcShape::N3) &&
                 unsigned(DstShape::N4) == unsigned(SrcShape::N4) &&
                 unsigned(DstShape::N5) == unsigned(SrcShape::N5) &&
                 unsigned(DstShape::N6) == unsigned(SrcShape::N6) &&
                 unsigned(DstShape::N7) == unsigned(SrcShape::N7) };
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< unsigned ScalarSize , unsigned Rank ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 , unsigned s7 ,
          typename iType >
KOKKOS_INLINE_FUNCTION
size_t dimension( 
  const Shape<ScalarSize,Rank,s0,s1,s2,s3,s4,s5,s6,s7> & shape ,
  const iType & r )
{
  return 0 == r ? shape.N0 : (
         1 == r ? shape.N1 : (
         2 == r ? shape.N2 : (
         3 == r ? shape.N3 : (
         4 == r ? shape.N4 : (
         5 == r ? shape.N5 : (
         6 == r ? shape.N6 : (
         7 == r ? shape.N7 : 1 )))))));
}

template< unsigned ScalarSize , unsigned Rank ,
          unsigned s0 , unsigned s1 , unsigned s2 , unsigned s3 ,
          unsigned s4 , unsigned s5 , unsigned s6 , unsigned s7 >
KOKKOS_INLINE_FUNCTION
size_t cardinality_count(
  const Shape<ScalarSize,Rank,s0,s1,s2,s3,s4,s5,s6,s7> & shape )
{
  return size_t(shape.N0) * shape.N1 * shape.N2 * shape.N3 *
         shape.N4 * shape.N5 * shape.N6 * shape.N7 ;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_CORESHAPE_HPP */

