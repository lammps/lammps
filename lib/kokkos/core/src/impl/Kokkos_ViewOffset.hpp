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

#ifndef KOKKOS_VIEWOFFSET_HPP
#define KOKKOS_VIEWOFFSET_HPP

#include <Kokkos_Pair.hpp>
#include <Kokkos_Layout.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
struct ALL ;
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos { namespace Impl {

template < class ShapeType , class LayoutType , typename Enable = void >
struct ViewOffset ;

//----------------------------------------------------------------------------
// LayoutLeft AND ( 1 >= rank OR 0 == rank_dynamic ) : no padding / striding
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutLeft
                 , typename enable_if<( 1 >= ShapeType::rank
                                        ||
                                        0 == ShapeType::rank_dynamic
                                      )>::type >
  : public ShapeType
{
  typedef size_t     size_type ;
  typedef ShapeType  shape_type ;
  typedef LayoutLeft array_layout ;

  enum { has_padding = false };

  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( size_t n )
    { assign_shape_dimension<R>( *this , n ); }

  // Return whether the subview introduced noncontiguity
  template< class S , class L >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if<( 0 == shape_type::rank &&
                             Impl::is_same<L,LayoutLeft>::value
                           ), bool >::type
  assign_subview( const ViewOffset<S,L,void> &
                , const size_t n0
                , const size_t n1
                , const size_t n2
                , const size_t n3
                , const size_t n4
                , const size_t n5
                , const size_t n6
                , const size_t n7
                )
    {
      return false ; // did not introduce noncontiguity
    }

  // This subview must be 1 == rank and 1 == rank_dynamic.
  // The source dimension #0 must be non-zero and all other dimensions are zero.
  // Return whether the subview introduced noncontiguity
  template< class S , class L >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if<( 1 == shape_type::rank &&
                             1 == shape_type::rank_dynamic &&
                             1 <= S::rank &&
                             Impl::is_same<L,LayoutLeft>::value
                           ), bool >::type
  assign_subview( const ViewOffset<S,L,void> &
                , const size_t n0
                , const size_t n1
                , const size_t n2
                , const size_t n3
                , const size_t n4
                , const size_t n5
                , const size_t n6
                , const size_t n7
                )
    {
      // n1 .. n7 must be zero
      shape_type::N0 = n0 ;
      return false ; // did not introduce noncontiguity
    }


  KOKKOS_INLINE_FUNCTION
  void assign( size_t n0 , unsigned n1 , unsigned n2 , unsigned n3
             , unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7
             , unsigned = 0 )
    { shape_type::assign( *this , n0, n1, n2, n3, n4, n5, n6, n7 ); }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutLeft > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                  )>::type * = 0 )
    { shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 ); }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutRight > & rhs
             , typename enable_if<( 1 == int(ShapeRHS::rank)
                                    &&
                                    1 == int(shape_type::rank)
                                    &&
                                    1 == int(shape_type::rank_dynamic)
                                  )>::type * = 0 )
    { shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 ); }

  KOKKOS_INLINE_FUNCTION
  void set_padding() {}

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const
    { return size_type(shape_type::N0) * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const
    { return size_type(shape_type::N0) * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  // Stride with [ rank ] value is the total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      s[0] = 1 ;
      if ( 0 < shape_type::rank ) { s[1] = shape_type::N0 ; }
      if ( 1 < shape_type::rank ) { s[2] = s[1] * shape_type::N1 ; }
      if ( 2 < shape_type::rank ) { s[3] = s[2] * shape_type::N2 ; }
      if ( 3 < shape_type::rank ) { s[4] = s[3] * shape_type::N3 ; }
      if ( 4 < shape_type::rank ) { s[5] = s[4] * shape_type::N4 ; }
      if ( 5 < shape_type::rank ) { s[6] = s[5] * shape_type::N5 ; }
      if ( 6 < shape_type::rank ) { s[7] = s[6] * shape_type::N6 ; }
      if ( 7 < shape_type::rank ) { s[8] = s[7] * shape_type::N7 ; }
    }

  KOKKOS_INLINE_FUNCTION size_type stride_0() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_1() const { return shape_type::N0 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_2() const { return shape_type::N0 * shape_type::N1 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_3() const { return shape_type::N0 * shape_type::N1 * shape_type::N2 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_4() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_5() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_6() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_7() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 ; }

  // rank 1
  template< typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const & i0 ) const { return i0 ; }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const & i0 , I1 const & i1 ) const
    { return i0 + shape_type::N0 * i1 ; }

  //rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0
                      , I1 const& i1
                      , I2 const& i2
                      ) const
    {
      return i0 + shape_type::N0 * (
             i1 + shape_type::N1 * i2 );
    }

  //rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3 ) const
    {
      return i0 + shape_type::N0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * i3 ));
    }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4 ) const
    {
      return i0 + shape_type::N0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * (
             i3 + shape_type::N3 * i4 )));
    }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5 ) const
    {
      return i0 + shape_type::N0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * (
             i3 + shape_type::N3 * (
             i4 + shape_type::N4 * i5 ))));
    }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6) const
  {
    return i0 + shape_type::N0 * (
           i1 + shape_type::N1 * (
           i2 + shape_type::N2 * (
           i3 + shape_type::N3 * (
           i4 + shape_type::N4 * (
           i5 + shape_type::N5 * i6 )))));
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7) const
  {
    return i0 + shape_type::N0 * (
           i1 + shape_type::N1 * (
           i2 + shape_type::N2 * (
           i3 + shape_type::N3 * (
           i4 + shape_type::N4 * (
           i5 + shape_type::N5 * (
           i6 + shape_type::N6 * i7 ))))));
  }
};

//----------------------------------------------------------------------------
// LayoutLeft AND ( 1 < rank AND 0 < rank_dynamic ) : has padding / striding
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutLeft
                 , typename enable_if<( 1 < ShapeType::rank
                                        &&
                                        0 < ShapeType::rank_dynamic
                                      )>::type >
  : public ShapeType
{
  typedef size_t     size_type ;
  typedef ShapeType  shape_type ;
  typedef LayoutLeft array_layout ;

  enum { has_padding = true };

  size_type S0 ;

  // This subview must be 2 == rank and 2 == rank_dynamic
  // due to only having stride #0.
  // The source dimension #0 must be non-zero for stride-one leading dimension.
  // If source is rank deficient then set to zero.
  // Return whether the subview introduced noncontiguity
  template< class S , class L >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if<( 2 == shape_type::rank &&
                             2 == shape_type::rank_dynamic &&
                             2 <= S::rank &&
                             Impl::is_same<L,LayoutLeft>::value
                           ), bool >::type
  assign_subview( const ViewOffset<S,L,void> & rhs
                , const size_t n0
                , const size_t n1
                , const size_t n2
                , const size_t n3
                , const size_t n4
                , const size_t n5
                , const size_t n6
                , const size_t n7
                )
    {
      // N0 = n0 ;
      // N1 = second non-zero dimension
      // S0 = stride for second non-zero dimension
      shape_type::N0 = 0 ;
      shape_type::N1 = 0 ;
      S0 = 0 ;

      if ( 0 == n0 ) {}
      else if (                n1 ) { shape_type::N0 = n0 ; shape_type::N1 = n1 ; S0 = rhs.stride_1(); }
      else if ( 2 < S::rank && n2 ) { shape_type::N0 = n0 ; shape_type::N1 = n2 ; S0 = rhs.stride_2(); }
      else if ( 3 < S::rank && n3 ) { shape_type::N0 = n0 ; shape_type::N1 = n3 ; S0 = rhs.stride_3(); }
      else if ( 4 < S::rank && n4 ) { shape_type::N0 = n0 ; shape_type::N1 = n4 ; S0 = rhs.stride_4(); }
      else if ( 5 < S::rank && n5 ) { shape_type::N0 = n0 ; shape_type::N1 = n5 ; S0 = rhs.stride_5(); }
      else if ( 6 < S::rank && n6 ) { shape_type::N0 = n0 ; shape_type::N1 = n6 ; S0 = rhs.stride_6(); }
      else if ( 7 < S::rank && n7 ) { shape_type::N0 = n0 ; shape_type::N1 = n7 ; S0 = rhs.stride_7(); }

      // Introduce noncontiguity if change the first dimension
      // or took a range of a dimension after the second.
      return ( size_t(shape_type::N0) != size_t(rhs.N0) ) || ( 0 == n1 );
    }


  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( size_t n )
    { assign_shape_dimension<R>( *this , n ); }


  KOKKOS_INLINE_FUNCTION
  void assign( size_t n0 , unsigned n1 , unsigned n2 , unsigned n3
             , unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7
             , unsigned = 0 )
    { shape_type::assign( *this , n0, n1, n2, n3, n4, n5, n6, n7 ); S0 = shape_type::N0 ; }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutLeft > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                    &&
                                    int(ShapeRHS::rank_dynamic) == 0
                                  )>::type * = 0 )
    {
      shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 );
      S0 = shape_type::N0 ; // No padding when dynamic_rank == 0
    }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutLeft > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                    &&
                                    int(ShapeRHS::rank_dynamic) > 0
                                  )>::type * = 0 )
    {
      shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 );
      S0 = rhs.S0 ; // possibly padding when dynamic rank > 0
    }

  KOKKOS_INLINE_FUNCTION
  void set_padding()
    {
      enum { div   = MEMORY_ALIGNMENT / shape_type::scalar_size };
      enum { mod   = MEMORY_ALIGNMENT % shape_type::scalar_size };
      enum { align = 0 == mod ? div : 0 };

      if ( align && MEMORY_ALIGNMENT_THRESHOLD * align < S0 ) {

        const size_type count_mod = S0 % ( div ? div : 1 );

        if ( count_mod ) { S0 += align - count_mod ; }
      }
    }

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const
    { return size_type(shape_type::N0) * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const
    { return size_type(S0) * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  // Stride with [ rank ] as total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      s[0] = 1 ;
      if ( 0 < shape_type::rank ) { s[1] = S0 ; }
      if ( 1 < shape_type::rank ) { s[2] = s[1] * shape_type::N1 ; }
      if ( 2 < shape_type::rank ) { s[3] = s[2] * shape_type::N2 ; }
      if ( 3 < shape_type::rank ) { s[4] = s[3] * shape_type::N3 ; }
      if ( 4 < shape_type::rank ) { s[5] = s[4] * shape_type::N4 ; }
      if ( 5 < shape_type::rank ) { s[6] = s[5] * shape_type::N5 ; }
      if ( 6 < shape_type::rank ) { s[7] = s[6] * shape_type::N6 ; }
      if ( 7 < shape_type::rank ) { s[8] = s[7] * shape_type::N6 ; }
    }

  KOKKOS_INLINE_FUNCTION size_type stride_0() const { return 1 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_1() const { return S0 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_2() const { return S0 * shape_type::N1 ; }
  KOKKOS_INLINE_FUNCTION size_type stride_3() const { return S0 * shape_type::N1 * shape_type::N2 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_4() const
    { return S0 * shape_type::N1 * shape_type::N2 * shape_type::N3 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_5() const
    { return S0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_6() const
    { return S0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_7() const
    { return S0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 ; }

  // rank 2
  template < typename I0 , typename I1 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const & i0 , I1 const & i1) const
    { return i0 + S0 * i1 ; }

  //rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 ) const
    {
      return i0 + S0 * (
             i1 + shape_type::N1 * i2 );
    }

  //rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3 ) const
    {
      return i0 + S0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * i3 ));
    }

  //rank 5
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4 ) const
    {
      return i0 + S0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * (
             i3 + shape_type::N3 * i4 )));
    }

  //rank 6
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5 ) const
    {
      return i0 + S0 * (
             i1 + shape_type::N1 * (
             i2 + shape_type::N2 * (
             i3 + shape_type::N3 * (
             i4 + shape_type::N4 * i5 ))));
    }

  //rank 7
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6 ) const
  {
    return i0 + S0 * (
           i1 + shape_type::N1 * (
           i2 + shape_type::N2 * (
           i3 + shape_type::N3 * (
           i4 + shape_type::N4 * (
           i5 + shape_type::N5 * i6 )))));
  }

  //rank 8
  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7 ) const
  {
    return i0 + S0 * (
           i1 + shape_type::N1 * (
           i2 + shape_type::N2 * (
           i3 + shape_type::N3 * (
           i4 + shape_type::N4 * (
           i5 + shape_type::N5 * (
           i6 + shape_type::N6 * i7 ))))));
  }
};

//----------------------------------------------------------------------------
// LayoutRight AND ( 1 >= rank OR 1 >= rank_dynamic ) : no padding / striding
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutRight
                 , typename enable_if<( 1 >= ShapeType::rank
                                        ||
                                        1 >= ShapeType::rank_dynamic
                                      )>::type >
  : public ShapeType
{
  typedef size_t       size_type;
  typedef ShapeType    shape_type;
  typedef LayoutRight  array_layout ;

  enum { has_padding = false };

  // This subview must be 1 == rank and 1 == rank_dynamic
  // The source view's last dimension must be non-zero
  // Return whether the subview introduced noncontiguity
  template< class S , class L >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if<( 0 == shape_type::rank &&
                             Impl::is_same<L,LayoutRight>::value
                           ), bool >::type
  assign_subview( const ViewOffset<S,L,void> &
                , const size_t n0
                , const size_t n1
                , const size_t n2
                , const size_t n3
                , const size_t n4
                , const size_t n5
                , const size_t n6
                , const size_t n7
                )
    { return false ; }

  // This subview must be 1 == rank and 1 == rank_dynamic
  // The source view's last dimension must be non-zero
  // Return whether the subview introduced noncontiguity
  template< class S , class L >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if<( 1 == shape_type::rank &&
                             1 == shape_type::rank_dynamic &&
                             1 <= S::rank &&
                             Impl::is_same<L,LayoutRight>::value
                           ), bool >::type
  assign_subview( const ViewOffset<S,L,void> &
                , const size_t n0
                , const size_t n1
                , const size_t n2
                , const size_t n3
                , const size_t n4
                , const size_t n5
                , const size_t n6
                , const size_t n7
                )
    {
      shape_type::N0 = S::rank == 1 ? n0 : (
                       S::rank == 2 ? n1 : (
                       S::rank == 3 ? n2 : (
                       S::rank == 4 ? n3 : (
                       S::rank == 5 ? n4 : (
                       S::rank == 6 ? n5 : (
                       S::rank == 7 ? n6 : n7 ))))));
      // should have n0 .. n_(rank-2) equal zero
      return false ;
    }

  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n )
    { assign_shape_dimension<R>( *this , n ); }

  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3
             , unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7
             , unsigned = 0 )
    { shape_type::assign( *this , n0, n1, n2, n3, n4, n5, n6, n7 ); }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutRight > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                  )>::type * = 0 )
    { shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 ); }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutLeft > & rhs
             , typename enable_if<( 1 == int(ShapeRHS::rank)
                                    &&
                                    1 == int(shape_type::rank)
                                    &&
                                    1 == int(shape_type::rank_dynamic)
                                  )>::type * = 0 )
    { shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 ); }

  KOKKOS_INLINE_FUNCTION
  void set_padding() {}

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const
    { return size_type(shape_type::N0) * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const
    { return size_type(shape_type::N0) * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  size_type stride_R() const
    {
      return size_type(shape_type::N1) * shape_type::N2 * shape_type::N3 *
             shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ;
    };

  // Stride with [rank] as total length
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      size_type n = 1 ;
      if ( 7 < shape_type::rank ) { s[7] = n ; n *= shape_type::N7 ; }
      if ( 6 < shape_type::rank ) { s[6] = n ; n *= shape_type::N6 ; }
      if ( 5 < shape_type::rank ) { s[5] = n ; n *= shape_type::N5 ; }
      if ( 4 < shape_type::rank ) { s[4] = n ; n *= shape_type::N4 ; }
      if ( 3 < shape_type::rank ) { s[3] = n ; n *= shape_type::N3 ; }
      if ( 2 < shape_type::rank ) { s[2] = n ; n *= shape_type::N2 ; }
      if ( 1 < shape_type::rank ) { s[1] = n ; n *= shape_type::N1 ; }
      if ( 0 < shape_type::rank ) { s[0] = n ; }
      s[shape_type::rank] = n * shape_type::N0 ;
    }

  KOKKOS_INLINE_FUNCTION
  size_type stride_7() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_6() const { return shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_5() const { return shape_type::N7 * shape_type::N6 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_4() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_3() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_2() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_1() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 * shape_type::N2 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_0() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 * shape_type::N2 * shape_type::N1 ; }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1 ) const
    {
      return i1 + shape_type::N1 * i0 ;
    }

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 ) const
    {
      return i2 + shape_type::N2 * (
             i1 + shape_type::N1 * ( i0 ));
    }

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3 ) const
    {
      return i3 + shape_type::N3 * (
             i2 + shape_type::N2 * (
             i1 + shape_type::N1 * ( i0 )));
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4 ) const
    {
      return i4 + shape_type::N4 * (
             i3 + shape_type::N3 * (
             i2 + shape_type::N2 * (
             i1 + shape_type::N1 * ( i0 ))));
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5 ) const
  {
    return i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * (
           i1 + shape_type::N1 * ( i0 )))));
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6 ) const
  {
    return i6 + shape_type::N6 * (
           i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * (
           i1 + shape_type::N1 * ( i0 ))))));
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7 ) const
  {
    return i7 + shape_type::N7 * (
           i6 + shape_type::N6 * (
           i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * (
           i1 + shape_type::N1 * ( i0 )))))));
  }
};

//----------------------------------------------------------------------------
// LayoutRight AND ( 1 < rank AND 1 < rank_dynamic ) : has padding / striding
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutRight
                 , typename enable_if<( 1 < ShapeType::rank
                                        &&
                                        1 < ShapeType::rank_dynamic
                                      )>::type >
  : public ShapeType
{
  typedef size_t       size_type;
  typedef ShapeType    shape_type;
  typedef LayoutRight  array_layout ;

  enum { has_padding = true };

  size_type SR ;

  // This subview must be 2 == rank and 2 == rank_dynamic
  // due to only having stride #(rank-1).
  // The source dimension #(rank-1) must be non-zero for stride-one leading dimension.
  // If source is rank deficient then set to zero.
  // Return whether the subview introduced noncontiguity
  template< class S , class L >
  KOKKOS_INLINE_FUNCTION
  typename Impl::enable_if<( 2 == shape_type::rank &&
                             2 == shape_type::rank_dynamic &&
                             2 <= S::rank &&
                             Impl::is_same<L,LayoutRight>::value
                           ), bool >::type
  assign_subview( const ViewOffset<S,L,void> & rhs
                , const size_t n0
                , const size_t n1
                , const size_t n2
                , const size_t n3
                , const size_t n4
                , const size_t n5
                , const size_t n6
                , const size_t n7
                )
    {
      const size_type nR = S::rank == 2 ? n1 : (
                           S::rank == 3 ? n2 : (
                           S::rank == 4 ? n3 : (
                           S::rank == 5 ? n4 : (
                           S::rank == 6 ? n5 : (
                           S::rank == 7 ? n6 : n7 )))));

      // N0 = first non-zero-dimension
      // N1 = last non-zero dimension
      // SR = stride for second non-zero dimension
      shape_type::N0 = 0 ;
      shape_type::N1 = 0 ;
      SR = 0 ;

      if ( 0 == nR ) {}
      else if (                n0 ) { shape_type::N0 = n0 ; shape_type::N1 = nR ; SR = rhs.stride_0(); }
      else if ( 2 < S::rank && n1 ) { shape_type::N0 = n1 ; shape_type::N1 = nR ; SR = rhs.stride_1(); }
      else if ( 3 < S::rank && n2 ) { shape_type::N0 = n2 ; shape_type::N1 = nR ; SR = rhs.stride_2(); }
      else if ( 4 < S::rank && n3 ) { shape_type::N0 = n3 ; shape_type::N1 = nR ; SR = rhs.stride_3(); }
      else if ( 5 < S::rank && n4 ) { shape_type::N0 = n4 ; shape_type::N1 = nR ; SR = rhs.stride_4(); }
      else if ( 6 < S::rank && n5 ) { shape_type::N0 = n5 ; shape_type::N1 = nR ; SR = rhs.stride_5(); }
      else if ( 7 < S::rank && n6 ) { shape_type::N0 = n6 ; shape_type::N1 = nR ; SR = rhs.stride_6(); }

      // Introduce noncontiguous if change the last dimension
      // or take a range of a dimension other than the second-to-last dimension.

      return 2 == S::rank ? ( size_t(shape_type::N1) != size_t(rhs.N1) || 0 == n0 ) : (
             3 == S::rank ? ( size_t(shape_type::N1) != size_t(rhs.N2) || 0 == n1 ) : (
             4 == S::rank ? ( size_t(shape_type::N1) != size_t(rhs.N3) || 0 == n2 ) : (
             5 == S::rank ? ( size_t(shape_type::N1) != size_t(rhs.N4) || 0 == n3 ) : (
             6 == S::rank ? ( size_t(shape_type::N1) != size_t(rhs.N5) || 0 == n4 ) : (
             7 == S::rank ? ( size_t(shape_type::N1) != size_t(rhs.N6) || 0 == n5 ) : (
                            ( size_t(shape_type::N1) != size_t(rhs.N7) || 0 == n6 ) ))))));
    }

  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n )
    { assign_shape_dimension<R>( *this , n ); }

  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n0 , unsigned n1 , unsigned n2 , unsigned n3
             , unsigned n4 , unsigned n5 , unsigned n6 , unsigned n7
             , unsigned = 0 )
    {
      shape_type::assign( *this , n0, n1, n2, n3, n4, n5, n6, n7 );
      SR = size_type(shape_type::N1) * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ;
    }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutRight > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= 1
                                  )>::type * = 0 )
    {
      shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 );
      SR = shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ;
    }

  template< class ShapeRHS >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset< ShapeRHS , LayoutRight > & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank)
                                    &&
                                    int(ShapeRHS::rank_dynamic) <= int(shape_type::rank_dynamic)
                                    &&
                                    int(ShapeRHS::rank_dynamic) > 1
                                  )>::type * = 0 )
    {
      shape_type::assign( *this , rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 );
      SR = rhs.SR ;
    }

  KOKKOS_INLINE_FUNCTION
  void set_padding()
    {
      enum { div   = MEMORY_ALIGNMENT / shape_type::scalar_size };
      enum { mod   = MEMORY_ALIGNMENT % shape_type::scalar_size };
      enum { align = 0 == mod ? div : 0 };

      if ( align && MEMORY_ALIGNMENT_THRESHOLD * align < SR ) {

        const size_type count_mod = SR % ( div ? div : 1 );

        if ( count_mod ) { SR += align - count_mod ; }
      }
    }

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const
    { return size_type(shape_type::N0) * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const { return shape_type::N0 * SR ; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    {
      size_type n = 1 ;
      if ( 7 < shape_type::rank ) { s[7] = n ; n *= shape_type::N7 ; }
      if ( 6 < shape_type::rank ) { s[6] = n ; n *= shape_type::N6 ; }
      if ( 5 < shape_type::rank ) { s[5] = n ; n *= shape_type::N5 ; }
      if ( 4 < shape_type::rank ) { s[4] = n ; n *= shape_type::N4 ; }
      if ( 3 < shape_type::rank ) { s[3] = n ; n *= shape_type::N3 ; }
      if ( 2 < shape_type::rank ) { s[2] = n ; n *= shape_type::N2 ; }
      if ( 1 < shape_type::rank ) { s[1] = n ; n *= shape_type::N1 ; }
      if ( 0 < shape_type::rank ) { s[0] = SR ; }
      s[shape_type::rank] = SR * shape_type::N0 ;
    }

  KOKKOS_INLINE_FUNCTION
  size_type stride_7() const { return 1 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_6() const { return shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_5() const { return shape_type::N7 * shape_type::N6 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_4() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_3() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_2() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_1() const { return shape_type::N7 * shape_type::N6 * shape_type::N5 * shape_type::N4 * shape_type::N3 * shape_type::N2 ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_0() const { return SR ; }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1 ) const
    {
      return i1 + i0 * SR ;
    }

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 ) const
    {
      return i2 + shape_type::N2 * ( i1 ) +
             i0 * SR ;
    }

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3 ) const
    {
      return i3 + shape_type::N3 * (
             i2 + shape_type::N2 * ( i1 )) +
             i0 * SR ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4 ) const
    {
      return i4 + shape_type::N4 * (
             i3 + shape_type::N3 * (
             i2 + shape_type::N2 * ( i1 ))) +
             i0 * SR ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5 ) const
  {
    return i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * ( i1 )))) +
           i0 * SR ;
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6 ) const
  {
    return i6 + shape_type::N6 * (
           i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * ( i1 ))))) +
           i0 * SR ;
  }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7 ) const
  {
    return i7 + shape_type::N7 * (
           i6 + shape_type::N6 * (
           i5 + shape_type::N5 * (
           i4 + shape_type::N4 * (
           i3 + shape_type::N3 * (
           i2 + shape_type::N2 * ( i1 )))))) +
           i0 * SR ;
  }
};

//----------------------------------------------------------------------------
// LayoutStride : 
template < class ShapeType >
struct ViewOffset< ShapeType , LayoutStride
                 , typename enable_if<( 0 < ShapeType::rank )>::type >
  : public ShapeType
{
  typedef size_t        size_type;
  typedef ShapeType     shape_type;
  typedef LayoutStride  array_layout ;

  size_type S[ shape_type::rank + 1 ];

  template< class SType , class L >
  KOKKOS_INLINE_FUNCTION
  bool assign_subview( const ViewOffset<SType,L,void> & rhs
                     , const size_type n0
                     , const size_type n1
                     , const size_type n2
                     , const size_type n3
                     , const size_type n4
                     , const size_type n5
                     , const size_type n6
                     , const size_type n7
                     )
    {
      shape_type::assign( *this, 0,0,0,0, 0,0,0,0 );

      for ( int i = 0 ; i < int(shape_type::rank+1) ; ++i ) { S[i] = 0 ; }

      // preconditions:
      //  shape_type::rank <= rhs.rank
      //  shape_type::rank == count of nonzero( rhs_dim[i] )
      size_type dim[8] = { n0 , n1 , n2 , n3 , n4 , n5 , n6 , n7 };
      size_type str[ SType::rank + 1 ];

      rhs.stride( str );

      // contract the zero-dimensions
      int r = 0 ;
      for ( int i = 0 ; i < int(SType::rank) ; ++i ) {
        if ( 0 != dim[i] ) {
          dim[r] = dim[i] ;
          str[r] = str[i] ;
          ++r ;
        }
      }

      if ( int(shape_type::rank) == r ) {
        // The shape is non-zero
        for ( int i = 0 ; i < int(shape_type::rank) ; ++i ) {
          const size_type cap = dim[i] * ( S[i] = str[i] );
          if ( S[ shape_type::rank ] < cap ) S[ shape_type::rank ] = cap ;
        }
        // set the contracted nonzero dimensions
        shape_type::assign( *this, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7] );
      }

      return true ; // definitely noncontiguous
    }

  template< unsigned R >
  KOKKOS_INLINE_FUNCTION
  void assign( unsigned n )
    { assign_shape_dimension<R>( *this , n ); }

  template< class ShapeRHS , class Layout >
  KOKKOS_INLINE_FUNCTION
  void assign( const ViewOffset<ShapeRHS,Layout> & rhs
             , typename enable_if<( int(ShapeRHS::rank) == int(shape_type::rank) )>::type * = 0 )
    {
      rhs.stride(S);
      shape_type::assign( *this, rhs.N0, rhs.N1, rhs.N2, rhs.N3, rhs.N4, rhs.N5, rhs.N6, rhs.N7 );
    }

  KOKKOS_INLINE_FUNCTION
  void assign( const LayoutStride & layout )
  {
    size_type max = 0 ;
    for ( int i = 0 ; i < shape_type::rank ; ++i ) {
      S[i] = layout.stride[i] ;
      const size_type m = layout.dimension[i] * S[i] ;
      if ( max < m ) { max = m ; }
    }
    S[ shape_type::rank ] = max ;
    shape_type::assign( *this, layout.dimension[0], layout.dimension[1],
                               layout.dimension[2], layout.dimension[3],
                               layout.dimension[4], layout.dimension[5],
                               layout.dimension[6], layout.dimension[7] );
  }

  KOKKOS_INLINE_FUNCTION
  void assign( size_t s0 , size_t s1 , size_t s2 , size_t s3
             , size_t s4 , size_t s5 , size_t s6 , size_t s7
             , size_t s8 )
    {
      const size_t str[9] = { s0, s1, s2, s3, s4, s5, s6, s7, s8 };

      // Last argument is the total length.
      // Total length must be non-zero.
      // All strides must be non-zero and less than total length.
      bool ok = 0 < str[ shape_type::rank ] ;

      for ( int i = 0 ; ( i < shape_type::rank ) &&
                        ( ok = 0 < str[i] && str[i] < str[ shape_type::rank ] ); ++i );

      if ( ok ) {
        size_t dim[8] = { 1,1,1,1,1,1,1,1 }; 
        int iorder[9] = { 0,0,0,0,0,0,0,0,0 }; 

        // Ordering of strides smallest to largest.
        for ( int i = 1 ; i < shape_type::rank ; ++i ) {
          int j = i ;
          for ( ; 0 < j && str[i] < str[ iorder[j-1] ] ; --j ) {
            iorder[j] = iorder[j-1] ;
          }
          iorder[j] = i ;
        }

        // Last argument is the total length.
        iorder[ shape_type::rank ] = shape_type::rank ;

        // Determine dimension associated with each stride.
        // Guarantees non-overlap by truncating dimension
        // if ( 0 != str[ iorder[i+1] ] % str[ iorder[i] ] )
        for ( int i = 0 ; i < shape_type::rank ; ++i ) {
          dim[ iorder[i] ] = str[ iorder[i+1] ] / str[ iorder[i] ] ;
        }

        // Assign dimensions and strides:
        shape_type::assign( *this, dim[0], dim[1], dim[2], dim[3], dim[4], dim[5], dim[6], dim[7] );
        for ( int i = 0 ; i <= shape_type::rank ; ++i ) { S[i] = str[i] ; }
      }
      else {
        shape_type::assign(*this,0,0,0,0,0,0,0,0);
        for ( int i = 0 ; i <= shape_type::rank ; ++i ) { S[i] = 0 ; }
      }
    }

  KOKKOS_INLINE_FUNCTION
  void set_padding() {}

  KOKKOS_INLINE_FUNCTION
  size_type cardinality() const
    { return shape_type::N0 * shape_type::N1 * shape_type::N2 * shape_type::N3 * shape_type::N4 * shape_type::N5 * shape_type::N6 * shape_type::N7 ; }

  KOKKOS_INLINE_FUNCTION
  size_type capacity() const { return S[ shape_type::rank ]; }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
    { for ( int i = 0 ; i <= shape_type::rank ; ++i ) { s[i] = S[i] ; } }

  KOKKOS_INLINE_FUNCTION
  size_type stride_0() const { return S[0] ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_1() const { return S[1] ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_2() const { return S[2] ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_3() const { return S[3] ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_4() const { return S[4] ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_5() const { return S[5] ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_6() const { return S[6] ; }

  KOKKOS_INLINE_FUNCTION
  size_type stride_7() const { return S[7] ; }

  // rank 1
  template <typename I0 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0 ) const
    {
      return i0 * S[0] ;
    }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1 ) const
    {
      return i0 * S[0] + i1 * S[1] ;
    }

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] ;
    }

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] + i4 * S[4] ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] + i4 * S[4] + i5 * S[5] ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] + i4 * S[4] + i5 * S[5] + i6 * S[6] ;
    }

  template < typename I0, typename I1, typename I2, typename I3
            ,typename I4, typename I5, typename I6, typename I7 >
  KOKKOS_FORCEINLINE_FUNCTION
  size_type operator()( I0 const& i0, I1 const& i1, I2 const& i2 , I3 const& i3, I4 const& i4, I5 const& i5, I6 const& i6, I7 const& i7 ) const
    {
      return i0 * S[0] + i1 * S[1] + i2 * S[2] + i3 * S[3] + i4 * S[4] + i5 * S[5] + i6 * S[6] + i7 * S[7] ;
    }
};

//----------------------------------------------------------------------------

template< class T >
struct ViewOffsetRange {

  enum { OK_integral_type = Impl::StaticAssert< Impl::is_integral<T>::value >::value };

  enum { is_range = false };

  KOKKOS_INLINE_FUNCTION static
  size_t dimension( size_t const , T const & ) { return 0 ; }

  KOKKOS_INLINE_FUNCTION static
  size_t begin( T const & i ) { return size_t(i) ; }
};

template<>
struct ViewOffsetRange<void> {
  enum { is_range = false };
};

template<>
struct ViewOffsetRange< Kokkos::ALL > {
  enum { is_range = true };

  KOKKOS_INLINE_FUNCTION static
  size_t dimension( size_t const n , ALL const & ) { return n ; }

  KOKKOS_INLINE_FUNCTION static
  size_t begin( ALL const & ) { return 0 ; }
};

template< typename iType >
struct ViewOffsetRange< std::pair<iType,iType> > {

  enum { OK_integral_type = Impl::StaticAssert< Impl::is_integral<iType>::value >::value };

  enum { is_range = true };

  KOKKOS_INLINE_FUNCTION static
  size_t dimension( size_t const n , std::pair<iType,iType> const & r )
    { return ( size_t(r.first) < size_t(r.second) && size_t(r.second) <= n ) ? size_t(r.second) - size_t(r.first) : 0 ; }

  KOKKOS_INLINE_FUNCTION static
  size_t begin( std::pair<iType,iType> const & r ) { return size_t(r.first) ; }
};

template< typename iType >
struct ViewOffsetRange< Kokkos::pair<iType,iType> > {

  enum { OK_integral_type = Impl::StaticAssert< Impl::is_integral<iType>::value >::value };

  enum { is_range = true };

  KOKKOS_INLINE_FUNCTION static
  size_t dimension( size_t const n , Kokkos::pair<iType,iType> const & r )
    { return ( size_t(r.first) < size_t(r.second) && size_t(r.second) <= n ) ? size_t(r.second) - size_t(r.first) : 0 ; }

  KOKKOS_INLINE_FUNCTION static
  size_t begin( Kokkos::pair<iType,iType> const & r ) { return size_t(r.first) ; }
};

}} // namespace Kokkos::Impl

#endif //KOKKOS_VIEWOFFSET_HPP

