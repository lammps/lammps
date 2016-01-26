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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <sstream>
#include <iostream>

/*--------------------------------------------------------------------------*/

#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

namespace Test {

template< typename T, class DeviceType >
class TestViewAPI {
public:
  TestViewAPI() {}
};

}

#else

/*--------------------------------------------------------------------------*/

namespace Test {

template< class T , class L , class D , class M , class S >
size_t allocation_count( const Kokkos::View<T,L,D,M,S> & view )
{
  const size_t card  = Kokkos::Impl::cardinality_count( view.shape() );
  const size_t alloc = view.capacity();

  return card <= alloc ? alloc : 0 ;
}

/*--------------------------------------------------------------------------*/

template< typename T, class DeviceType>
struct TestViewOperator
{
  typedef DeviceType  execution_space ;

  static const unsigned N = 100 ;
  static const unsigned D = 3 ;

  typedef Kokkos::View< T*[D] , execution_space > view_type ;

  const view_type v1 ;
  const view_type v2 ;

  TestViewOperator()
    : v1( "v1" , N )
    , v2( "v2" , N )
    {}

  static void testit()
  {
    Kokkos::parallel_for( N , TestViewOperator() );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned i ) const
  {
    const unsigned X = 0 ;
    const unsigned Y = 1 ;
    const unsigned Z = 2 ;

    v2(i,X) = v1(i,X);
    v2(i,Y) = v1(i,Y);
    v2(i,Z) = v1(i,Z);
  }
};

/*--------------------------------------------------------------------------*/

template< class DataType >
struct rank {
private:
  typedef typename Kokkos::Impl::AnalyzeShape<DataType>::shape shape ;
public:
  static const unsigned value = shape::rank ;
};

template< class DataType ,
          class DeviceType ,
          unsigned Rank = rank< DataType >::value >
struct TestViewOperator_LeftAndRight ;

template< class DataType , class DeviceType >
struct TestViewOperator_LeftAndRight< DataType , DeviceType , 8 >
{
  typedef DeviceType                          execution_space ;
  typedef typename execution_space::memory_space  memory_space ;
  typedef typename execution_space::size_type     size_type ;

  typedef int value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update |= input ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }


  typedef Kokkos::
    View< DataType, Kokkos::LayoutLeft, execution_space > left_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutRight, execution_space > right_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutStride, execution_space > stride_view ;

  typedef typename left_view ::shape_type  left_shape ;
  typedef typename right_view::shape_type  right_shape ;

  left_shape   lsh ;
  right_shape  rsh ;
  left_view    left ;
  right_view   right ;
  stride_view  left_stride ;
  stride_view  right_stride ;
  long         left_alloc ;
  long         right_alloc ;

  TestViewOperator_LeftAndRight()
    : lsh()
    , rsh()
    , left(  "left" )
    , right( "right" )
    , left_stride( left )
    , right_stride( right )
    , left_alloc( allocation_count( left ) )
    , right_alloc( allocation_count( right ) )
    {}

  static void testit()
  {
    TestViewOperator_LeftAndRight driver ;

    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.lsh ) <= driver.left_alloc );
    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.rsh ) <= driver.right_alloc );

    int error_flag = 0 ;

    Kokkos::parallel_reduce( 1 , driver , error_flag );

    ASSERT_EQ( error_flag , 0 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type , value_type & update ) const
  {
    long offset ;

    offset = -1 ;
    for ( unsigned i7 = 0 ; i7 < unsigned(lsh.N7) ; ++i7 )
    for ( unsigned i6 = 0 ; i6 < unsigned(lsh.N6) ; ++i6 )
    for ( unsigned i5 = 0 ; i5 < unsigned(lsh.N5) ; ++i5 )
    for ( unsigned i4 = 0 ; i4 < unsigned(lsh.N4) ; ++i4 )
    for ( unsigned i3 = 0 ; i3 < unsigned(lsh.N3) ; ++i3 )
    for ( unsigned i2 = 0 ; i2 < unsigned(lsh.N2) ; ++i2 )
    for ( unsigned i1 = 0 ; i1 < unsigned(lsh.N1) ; ++i1 )
    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    {
      const long j = & left( i0, i1, i2, i3, i4, i5, i6, i7 ) -
                     & left(  0,  0,  0,  0,  0,  0,  0,  0 );
      if ( j <= offset || left_alloc <= j ) { update |= 1 ; }
      offset = j ;

      if ( & left(i0,i1,i2,i3,i4,i5,i6,i7) !=
           & left_stride(i0,i1,i2,i3,i4,i5,i6,i7) ) {
        update |= 4 ;
      }
    }

    offset = -1 ;
    for ( unsigned i0 = 0 ; i0 < unsigned(rsh.N0) ; ++i0 )
    for ( unsigned i1 = 0 ; i1 < unsigned(rsh.N1) ; ++i1 )
    for ( unsigned i2 = 0 ; i2 < unsigned(rsh.N2) ; ++i2 )
    for ( unsigned i3 = 0 ; i3 < unsigned(rsh.N3) ; ++i3 )
    for ( unsigned i4 = 0 ; i4 < unsigned(rsh.N4) ; ++i4 )
    for ( unsigned i5 = 0 ; i5 < unsigned(rsh.N5) ; ++i5 )
    for ( unsigned i6 = 0 ; i6 < unsigned(rsh.N6) ; ++i6 )
    for ( unsigned i7 = 0 ; i7 < unsigned(rsh.N7) ; ++i7 )
    {
      const long j = & right( i0, i1, i2, i3, i4, i5, i6, i7 ) -
                     & right(  0,  0,  0,  0,  0,  0,  0,  0 );
      if ( j <= offset || right_alloc <= j ) { update |= 2 ; }
      offset = j ;

      if ( & right(i0,i1,i2,i3,i4,i5,i6,i7) !=
           & right_stride(i0,i1,i2,i3,i4,i5,i6,i7) ) {
        update |= 8 ;
      }
    }
  }
};

template< class DataType , class DeviceType >
struct TestViewOperator_LeftAndRight< DataType , DeviceType , 7 >
{
  typedef DeviceType                          execution_space ;
  typedef typename execution_space::memory_space  memory_space ;
  typedef typename execution_space::size_type     size_type ;

  typedef int value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update |= input ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }


  typedef Kokkos::
    View< DataType, Kokkos::LayoutLeft, execution_space > left_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutRight, execution_space > right_view ;

  typedef typename left_view ::shape_type  left_shape ;
  typedef typename right_view::shape_type  right_shape ;

  left_shape   lsh ;
  right_shape  rsh ;
  left_view    left ;
  right_view   right ;
  long         left_alloc ;
  long         right_alloc ;

  TestViewOperator_LeftAndRight()
    : lsh()
    , rsh()
    , left(  "left" )
    , right( "right" )
    , left_alloc( allocation_count( left ) )
    , right_alloc( allocation_count( right ) )
    {}

  static void testit()
  {
    TestViewOperator_LeftAndRight driver ;

    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.lsh ) <= driver.left_alloc );
    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.rsh ) <= driver.right_alloc );

    int error_flag = 0 ;

    Kokkos::parallel_reduce( 1 , driver , error_flag );

    ASSERT_EQ( error_flag , 0 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type , value_type & update ) const
  {
    long offset ;

    offset = -1 ;
    for ( unsigned i6 = 0 ; i6 < unsigned(lsh.N6) ; ++i6 )
    for ( unsigned i5 = 0 ; i5 < unsigned(lsh.N5) ; ++i5 )
    for ( unsigned i4 = 0 ; i4 < unsigned(lsh.N4) ; ++i4 )
    for ( unsigned i3 = 0 ; i3 < unsigned(lsh.N3) ; ++i3 )
    for ( unsigned i2 = 0 ; i2 < unsigned(lsh.N2) ; ++i2 )
    for ( unsigned i1 = 0 ; i1 < unsigned(lsh.N1) ; ++i1 )
    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    {
      const long j = & left( i0, i1, i2, i3, i4, i5, i6 ) -
                     & left(  0,  0,  0,  0,  0,  0,  0 );
      if ( j <= offset || left_alloc <= j ) { update |= 1 ; }
      offset = j ;
    }

    offset = -1 ;
    for ( unsigned i0 = 0 ; i0 < unsigned(rsh.N0) ; ++i0 )
    for ( unsigned i1 = 0 ; i1 < unsigned(rsh.N1) ; ++i1 )
    for ( unsigned i2 = 0 ; i2 < unsigned(rsh.N2) ; ++i2 )
    for ( unsigned i3 = 0 ; i3 < unsigned(rsh.N3) ; ++i3 )
    for ( unsigned i4 = 0 ; i4 < unsigned(rsh.N4) ; ++i4 )
    for ( unsigned i5 = 0 ; i5 < unsigned(rsh.N5) ; ++i5 )
    for ( unsigned i6 = 0 ; i6 < unsigned(rsh.N6) ; ++i6 )
    {
      const long j = & right( i0, i1, i2, i3, i4, i5, i6 ) -
                     & right(  0,  0,  0,  0,  0,  0,  0 );
      if ( j <= offset || right_alloc <= j ) { update |= 2 ; }
      offset = j ;
    }
  }
};

template< class DataType , class DeviceType >
struct TestViewOperator_LeftAndRight< DataType , DeviceType , 6 >
{
  typedef DeviceType                          execution_space ;
  typedef typename execution_space::memory_space  memory_space ;
  typedef typename execution_space::size_type     size_type ;

  typedef int value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update |= input ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }


  typedef Kokkos::
    View< DataType, Kokkos::LayoutLeft, execution_space > left_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutRight, execution_space > right_view ;

  typedef typename left_view ::shape_type  left_shape ;
  typedef typename right_view::shape_type  right_shape ;

  left_shape   lsh ;
  right_shape  rsh ;
  left_view    left ;
  right_view   right ;
  long         left_alloc ;
  long         right_alloc ;

  TestViewOperator_LeftAndRight()
    : lsh()
    , rsh()
    , left(  "left" )
    , right( "right" )
    , left_alloc( allocation_count( left ) )
    , right_alloc( allocation_count( right ) )
    {}

  static void testit()
  {
    TestViewOperator_LeftAndRight driver ;

    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.lsh ) <= driver.left_alloc );
    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.rsh ) <= driver.right_alloc );

    int error_flag = 0 ;

    Kokkos::parallel_reduce( 1 , driver , error_flag );

    ASSERT_EQ( error_flag , 0 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type , value_type & update ) const
  {
    long offset ;

    offset = -1 ;
    for ( unsigned i5 = 0 ; i5 < unsigned(lsh.N5) ; ++i5 )
    for ( unsigned i4 = 0 ; i4 < unsigned(lsh.N4) ; ++i4 )
    for ( unsigned i3 = 0 ; i3 < unsigned(lsh.N3) ; ++i3 )
    for ( unsigned i2 = 0 ; i2 < unsigned(lsh.N2) ; ++i2 )
    for ( unsigned i1 = 0 ; i1 < unsigned(lsh.N1) ; ++i1 )
    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    {
      const long j = & left( i0, i1, i2, i3, i4, i5 ) -
                     & left(  0,  0,  0,  0,  0,  0 );
      if ( j <= offset || left_alloc <= j ) { update |= 1 ; }
      offset = j ;
    }

    offset = -1 ;
    for ( unsigned i0 = 0 ; i0 < unsigned(rsh.N0) ; ++i0 )
    for ( unsigned i1 = 0 ; i1 < unsigned(rsh.N1) ; ++i1 )
    for ( unsigned i2 = 0 ; i2 < unsigned(rsh.N2) ; ++i2 )
    for ( unsigned i3 = 0 ; i3 < unsigned(rsh.N3) ; ++i3 )
    for ( unsigned i4 = 0 ; i4 < unsigned(rsh.N4) ; ++i4 )
    for ( unsigned i5 = 0 ; i5 < unsigned(rsh.N5) ; ++i5 )
    {
      const long j = & right( i0, i1, i2, i3, i4, i5 ) -
                     & right(  0,  0,  0,  0,  0,  0 );
      if ( j <= offset || right_alloc <= j ) { update |= 2 ; }
      offset = j ;
    }
  }
};

template< class DataType , class DeviceType >
struct TestViewOperator_LeftAndRight< DataType , DeviceType , 5 >
{
  typedef DeviceType                          execution_space ;
  typedef typename execution_space::memory_space  memory_space ;
  typedef typename execution_space::size_type     size_type ;

  typedef int value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update |= input ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }


  typedef Kokkos::
    View< DataType, Kokkos::LayoutLeft, execution_space > left_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutRight, execution_space > right_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutStride, execution_space > stride_view ;

  typedef typename left_view ::shape_type  left_shape ;
  typedef typename right_view::shape_type  right_shape ;

  left_shape   lsh ;
  right_shape  rsh ;
  left_view    left ;
  right_view   right ;
  stride_view  left_stride ;
  stride_view  right_stride ;
  long         left_alloc ;
  long         right_alloc ;

  TestViewOperator_LeftAndRight()
    : lsh()
    , rsh()
    , left(  "left" )
    , right( "right" )
    , left_stride( left )
    , right_stride( right )
    , left_alloc( allocation_count( left ) )
    , right_alloc( allocation_count( right ) )
    {}

  static void testit()
  {
    TestViewOperator_LeftAndRight driver ;

    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.lsh ) <= driver.left_alloc );
    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.rsh ) <= driver.right_alloc );

    int error_flag = 0 ;

    Kokkos::parallel_reduce( 1 , driver , error_flag );

    ASSERT_EQ( error_flag , 0 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type , value_type & update ) const
  {
    long offset ;

    offset = -1 ;
    for ( unsigned i4 = 0 ; i4 < unsigned(lsh.N4) ; ++i4 )
    for ( unsigned i3 = 0 ; i3 < unsigned(lsh.N3) ; ++i3 )
    for ( unsigned i2 = 0 ; i2 < unsigned(lsh.N2) ; ++i2 )
    for ( unsigned i1 = 0 ; i1 < unsigned(lsh.N1) ; ++i1 )
    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    {
      const long j = & left( i0, i1, i2, i3, i4 ) -
                     & left(  0,  0,  0,  0,  0 );
      if ( j <= offset || left_alloc <= j ) { update |= 1 ; }
      offset = j ;

      if ( & left( i0, i1, i2, i3, i4 ) !=
           & left_stride( i0, i1, i2, i3, i4 ) ) { update |= 4 ; }
    }

    offset = -1 ;
    for ( unsigned i0 = 0 ; i0 < unsigned(rsh.N0) ; ++i0 )
    for ( unsigned i1 = 0 ; i1 < unsigned(rsh.N1) ; ++i1 )
    for ( unsigned i2 = 0 ; i2 < unsigned(rsh.N2) ; ++i2 )
    for ( unsigned i3 = 0 ; i3 < unsigned(rsh.N3) ; ++i3 )
    for ( unsigned i4 = 0 ; i4 < unsigned(rsh.N4) ; ++i4 )
    {
      const long j = & right( i0, i1, i2, i3, i4 ) -
                     & right(  0,  0,  0,  0,  0 );
      if ( j <= offset || right_alloc <= j ) { update |= 2 ; }
      offset = j ;

      if ( & right( i0, i1, i2, i3, i4 ) !=
           & right_stride( i0, i1, i2, i3, i4 ) ) { update |= 8 ; }
    }
  }
};

template< class DataType , class DeviceType >
struct TestViewOperator_LeftAndRight< DataType , DeviceType , 4 >
{
  typedef DeviceType                          execution_space ;
  typedef typename execution_space::memory_space  memory_space ;
  typedef typename execution_space::size_type     size_type ;

  typedef int value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update |= input ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }


  typedef Kokkos::
    View< DataType, Kokkos::LayoutLeft, execution_space > left_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutRight, execution_space > right_view ;

  typedef typename left_view ::shape_type  left_shape ;
  typedef typename right_view::shape_type  right_shape ;

  left_shape   lsh ;
  right_shape  rsh ;
  left_view    left ;
  right_view   right ;
  long         left_alloc ;
  long         right_alloc ;

  TestViewOperator_LeftAndRight()
    : lsh()
    , rsh()
    , left(  "left" )
    , right( "right" )
    , left_alloc( allocation_count( left ) )
    , right_alloc( allocation_count( right ) )
    {}

  static void testit()
  {
    TestViewOperator_LeftAndRight driver ;

    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.lsh ) <= driver.left_alloc );
    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.rsh ) <= driver.right_alloc );

    int error_flag = 0 ;

    Kokkos::parallel_reduce( 1 , driver , error_flag );

    ASSERT_EQ( error_flag , 0 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type , value_type & update ) const
  {
    long offset ;

    offset = -1 ;
    for ( unsigned i3 = 0 ; i3 < unsigned(lsh.N3) ; ++i3 )
    for ( unsigned i2 = 0 ; i2 < unsigned(lsh.N2) ; ++i2 )
    for ( unsigned i1 = 0 ; i1 < unsigned(lsh.N1) ; ++i1 )
    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    {
      const long j = & left( i0, i1, i2, i3 ) -
                     & left(  0,  0,  0,  0 );
      if ( j <= offset || left_alloc <= j ) { update |= 1 ; }
      offset = j ;
    }

    offset = -1 ;
    for ( unsigned i0 = 0 ; i0 < unsigned(rsh.N0) ; ++i0 )
    for ( unsigned i1 = 0 ; i1 < unsigned(rsh.N1) ; ++i1 )
    for ( unsigned i2 = 0 ; i2 < unsigned(rsh.N2) ; ++i2 )
    for ( unsigned i3 = 0 ; i3 < unsigned(rsh.N3) ; ++i3 )
    {
      const long j = & right( i0, i1, i2, i3 ) -
                     & right(  0,  0,  0,  0 );
      if ( j <= offset || right_alloc <= j ) { update |= 2 ; }
      offset = j ;
    }
  }
};

template< class DataType , class DeviceType >
struct TestViewOperator_LeftAndRight< DataType , DeviceType , 3 >
{
  typedef DeviceType                          execution_space ;
  typedef typename execution_space::memory_space  memory_space ;
  typedef typename execution_space::size_type     size_type ;

  typedef int value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update |= input ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }


  typedef Kokkos::
    View< DataType, Kokkos::LayoutLeft, execution_space > left_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutRight, execution_space > right_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutStride, execution_space > stride_view ;

  typedef typename left_view ::shape_type  left_shape ;
  typedef typename right_view::shape_type  right_shape ;

  left_shape   lsh ;
  right_shape  rsh ;
  left_view    left ;
  right_view   right ;
  stride_view  left_stride ;
  stride_view  right_stride ;
  long         left_alloc ;
  long         right_alloc ;

  TestViewOperator_LeftAndRight()
    : lsh()
    , rsh()
    , left(  std::string("left") )
    , right( std::string("right") )
    , left_stride( left )
    , right_stride( right )
    , left_alloc( allocation_count( left ) )
    , right_alloc( allocation_count( right ) )
    {}

  static void testit()
  {
    TestViewOperator_LeftAndRight driver ;

    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.lsh ) <= driver.left_alloc );
    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.rsh ) <= driver.right_alloc );

    int error_flag = 0 ;

    Kokkos::parallel_reduce( 1 , driver , error_flag );

    ASSERT_EQ( error_flag , 0 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type , value_type & update ) const
  {
    long offset ;

    offset = -1 ;
    for ( unsigned i2 = 0 ; i2 < unsigned(lsh.N2) ; ++i2 )
    for ( unsigned i1 = 0 ; i1 < unsigned(lsh.N1) ; ++i1 )
    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    {
      const long j = & left( i0, i1, i2 ) -
                     & left(  0,  0,  0 );
      if ( j <= offset || left_alloc <= j ) { update |= 1 ; }
      offset = j ;

      if ( & left(i0,i1,i2) != & left_stride(i0,i1,i2) ) { update |= 4 ; }
    }

    offset = -1 ;
    for ( unsigned i0 = 0 ; i0 < unsigned(rsh.N0) ; ++i0 )
    for ( unsigned i1 = 0 ; i1 < unsigned(rsh.N1) ; ++i1 )
    for ( unsigned i2 = 0 ; i2 < unsigned(rsh.N2) ; ++i2 )
    {
      const long j = & right( i0, i1, i2 ) -
                     & right(  0,  0,  0 );
      if ( j <= offset || right_alloc <= j ) { update |= 2 ; }
      offset = j ;

      if ( & right(i0,i1,i2) != & right_stride(i0,i1,i2) ) { update |= 8 ; }
    }

    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    for ( unsigned i1 = 0 ; i1 < unsigned(lsh.N1) ; ++i1 )
    for ( unsigned i2 = 0 ; i2 < unsigned(lsh.N2) ; ++i2 )
    {
      if ( & left(i0,i1,i2)  != & left.at(i0,i1,i2,0,0,0,0,0) )  { update |= 3 ; }
      if ( & right(i0,i1,i2) != & right.at(i0,i1,i2,0,0,0,0,0) ) { update |= 3 ; }
    }
  }
};

template< class DataType , class DeviceType >
struct TestViewOperator_LeftAndRight< DataType , DeviceType , 2 >
{
  typedef DeviceType                          execution_space ;
  typedef typename execution_space::memory_space  memory_space ;
  typedef typename execution_space::size_type     size_type ;

  typedef int value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update |= input ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }


  typedef Kokkos::
    View< DataType, Kokkos::LayoutLeft, execution_space > left_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutRight, execution_space > right_view ;

  typedef typename left_view ::shape_type  left_shape ;
  typedef typename right_view::shape_type  right_shape ;

  left_shape   lsh ;
  right_shape  rsh ;
  left_view    left ;
  right_view   right ;
  long         left_alloc ;
  long         right_alloc ;

  TestViewOperator_LeftAndRight()
    : lsh()
    , rsh()
    , left(  Kokkos::ViewAllocate("left") )
    , right( Kokkos::ViewAllocate("right") )
    , left_alloc( allocation_count( left ) )
    , right_alloc( allocation_count( right ) )
    {}

  static void testit()
  {
    TestViewOperator_LeftAndRight driver ;

    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.lsh ) <= driver.left_alloc );
    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.rsh ) <= driver.right_alloc );

    int error_flag = 0 ;

    Kokkos::parallel_reduce( 1 , driver , error_flag );

    ASSERT_EQ( error_flag , 0 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type , value_type & update ) const
  {
    long offset ;

    offset = -1 ;
    for ( unsigned i1 = 0 ; i1 < unsigned(lsh.N1) ; ++i1 )
    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    {
      const long j = & left( i0, i1 ) -
                     & left(  0,  0 );
      if ( j <= offset || left_alloc <= j ) { update |= 1 ; }
      offset = j ;
    }

    offset = -1 ;
    for ( unsigned i0 = 0 ; i0 < unsigned(rsh.N0) ; ++i0 )
    for ( unsigned i1 = 0 ; i1 < unsigned(rsh.N1) ; ++i1 )
    {
      const long j = & right( i0, i1 ) -
                     & right(  0,  0 );
      if ( j <= offset || right_alloc <= j ) { update |= 2 ; }
      offset = j ;
    }

    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    for ( unsigned i1 = 0 ; i1 < unsigned(lsh.N1) ; ++i1 )
    {
      if ( & left(i0,i1)  != & left.at(i0,i1,0,0,0,0,0,0) )  { update |= 3 ; }
      if ( & right(i0,i1) != & right.at(i0,i1,0,0,0,0,0,0) ) { update |= 3 ; }
    }
  }
};

template< class DataType , class DeviceType >
struct TestViewOperator_LeftAndRight< DataType , DeviceType , 1 >
{
  typedef DeviceType                          execution_space ;
  typedef typename execution_space::memory_space  memory_space ;
  typedef typename execution_space::size_type     size_type ;

  typedef int value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    const volatile value_type & input )
    { update |= input ; }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & update )
    { update = 0 ; }


  typedef Kokkos::
    View< DataType, Kokkos::LayoutLeft, execution_space > left_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutRight, execution_space > right_view ;

  typedef Kokkos::
    View< DataType, Kokkos::LayoutStride, execution_space > stride_view ;

  typedef typename left_view ::shape_type  left_shape ;
  typedef typename right_view::shape_type  right_shape ;

  left_shape   lsh ;
  right_shape  rsh ;
  left_view    left ;
  right_view   right ;
  stride_view  left_stride ;
  stride_view  right_stride ;
  long         left_alloc ;
  long         right_alloc ;

  TestViewOperator_LeftAndRight()
    : lsh()
    , rsh()
    , left(  Kokkos::ViewAllocate() )
    , right( Kokkos::ViewAllocate() )
    , left_stride( left )
    , right_stride( right )
    , left_alloc( allocation_count( left ) )
    , right_alloc( allocation_count( right ) )
    {}

  static void testit()
  {
    TestViewOperator_LeftAndRight driver ;

    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.lsh ) <= driver.left_alloc );
    ASSERT_TRUE( (long) Kokkos::Impl::cardinality_count( driver.rsh ) <= driver.right_alloc );

    int error_flag = 0 ;

    Kokkos::parallel_reduce( 1 , driver , error_flag );

    ASSERT_EQ( error_flag , 0 );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const size_type , value_type & update ) const
  {
    for ( unsigned i0 = 0 ; i0 < unsigned(lsh.N0) ; ++i0 )
    {
      if ( & left(i0)  != & left.at(i0,0,0,0,0,0,0,0) )  { update |= 3 ; }
      if ( & right(i0) != & right.at(i0,0,0,0,0,0,0,0) ) { update |= 3 ; }
      if ( & left(i0)  != & left_stride(i0) ) { update |= 4 ; }
      if ( & right(i0) != & right_stride(i0) ) { update |= 8 ; }
    }
  }
};

/*--------------------------------------------------------------------------*/

template< typename T, class DeviceType >
class TestViewAPI
{
public:
  typedef DeviceType        device ;

  enum { N0 = 1000 ,
         N1 = 3 ,
         N2 = 5 ,
         N3 = 7 };

  typedef Kokkos::View< T , device > dView0 ;
  typedef Kokkos::View< T* , device > dView1 ;
  typedef Kokkos::View< T*[N1] , device > dView2 ;
  typedef Kokkos::View< T*[N1][N2] , device > dView3 ;
  typedef Kokkos::View< T*[N1][N2][N3] , device > dView4 ;
  typedef Kokkos::View< const T*[N1][N2][N3] , device > const_dView4 ;

  typedef Kokkos::View< T****, device, Kokkos::MemoryUnmanaged > dView4_unmanaged ;

  typedef typename dView0::host_mirror_space host ;

  TestViewAPI()
  {
    run_test_mirror();
    run_test();
    run_test_scalar();
    run_test_const();
    run_test_subview();
    run_test_subview_strided();
    run_test_vector();

    TestViewOperator< T , device >::testit();
    TestViewOperator_LeftAndRight< int[2][3][4][2][3][4][2][3] , device >::testit();
    TestViewOperator_LeftAndRight< int[2][3][4][2][3][4][2] , device >::testit();
    TestViewOperator_LeftAndRight< int[2][3][4][2][3][4] , device >::testit();
    TestViewOperator_LeftAndRight< int[2][3][4][2][3] , device >::testit();
    TestViewOperator_LeftAndRight< int[2][3][4][2] , device >::testit();
    TestViewOperator_LeftAndRight< int[2][3][4] , device >::testit();
    TestViewOperator_LeftAndRight< int[2][3] , device >::testit();
    TestViewOperator_LeftAndRight< int[2] , device >::testit();
  }

  static void run_test_mirror()
  {
    typedef Kokkos::View< int , host > view_type ;
    typedef typename view_type::HostMirror mirror_type ;
    view_type a("a");
    mirror_type am = Kokkos::create_mirror_view(a);
    mirror_type ax = Kokkos::create_mirror(a);
    ASSERT_EQ( & a() , & am() );
  }

  static void run_test_scalar()
  {
    typedef typename dView0::HostMirror  hView0 ;

    dView0 dx , dy ;
    hView0 hx , hy ;

    dx = dView0( "dx" );
    dy = dView0( "dy" );

    hx = Kokkos::create_mirror( dx );
    hy = Kokkos::create_mirror( dy );

    hx = 1 ;

    Kokkos::deep_copy( dx , hx );
    Kokkos::deep_copy( dy , dx );
    Kokkos::deep_copy( hy , dy );

    ASSERT_EQ( hx(), hy() );
  }

  static void run_test()
  {
    // mfh 14 Feb 2014: This test doesn't actually create instances of
    // these types.  In order to avoid "declared but unused typedef"
    // warnings, we declare empty instances of these types, with the
    // usual "(void)" marker to avoid compiler warnings for unused
    // variables.

    typedef typename dView0::HostMirror  hView0 ;
    typedef typename dView1::HostMirror  hView1 ;
    typedef typename dView2::HostMirror  hView2 ;
    typedef typename dView3::HostMirror  hView3 ;
    typedef typename dView4::HostMirror  hView4 ;

    {
      hView0 thing;
      (void) thing;
    }
    {
      hView1 thing;
      (void) thing;
    }
    {
      hView2 thing;
      (void) thing;
    }
    {
      hView3 thing;
      (void) thing;
    }
    {
      hView4 thing;
      (void) thing;
    }

    dView4 dx , dy , dz ;
    hView4 hx , hy , hz ;

    ASSERT_TRUE( dx.is_null() );
    ASSERT_TRUE( dy.is_null() );
    ASSERT_TRUE( dz.is_null() );
    ASSERT_TRUE( hx.is_null() );
    ASSERT_TRUE( hy.is_null() );
    ASSERT_TRUE( hz.is_null() );
    ASSERT_EQ( dx.dimension_0() , 0u );
    ASSERT_EQ( dy.dimension_0() , 0u );
    ASSERT_EQ( dz.dimension_0() , 0u );
    ASSERT_EQ( hx.dimension_0() , 0u );
    ASSERT_EQ( hy.dimension_0() , 0u );
    ASSERT_EQ( hz.dimension_0() , 0u );
    ASSERT_EQ( dx.dimension_1() , unsigned(N1) );
    ASSERT_EQ( dy.dimension_1() , unsigned(N1) );
    ASSERT_EQ( dz.dimension_1() , unsigned(N1) );
    ASSERT_EQ( hx.dimension_1() , unsigned(N1) );
    ASSERT_EQ( hy.dimension_1() , unsigned(N1) );
    ASSERT_EQ( hz.dimension_1() , unsigned(N1) );

    dx = dView4( "dx" , N0 );
    dy = dView4( "dy" , N0 );

    #ifndef KOKKOS_USING_EXPERIMENTAL_VIEW
    ASSERT_EQ( dx.tracker().ref_count() , size_t(1) );
    #endif

    dView4_unmanaged unmanaged_dx = dx;
    #ifndef KOKKOS_USING_EXPERIMENTAL_VIEW
    ASSERT_EQ( dx.tracker().ref_count() , size_t(1) );
    #endif

    dView4_unmanaged unmanaged_from_ptr_dx = dView4_unmanaged(dx.ptr_on_device(),
                                                              dx.dimension_0(),
                                                              dx.dimension_1(),
                                                              dx.dimension_2(),
                                                              dx.dimension_3());

    {
      // Destruction of this view should be harmless
      const_dView4 unmanaged_from_ptr_const_dx( dx.ptr_on_device() ,
                                                dx.dimension_0() ,
                                                dx.dimension_1() ,
                                                dx.dimension_2() ,
                                                dx.dimension_3() );
    }

    const_dView4 const_dx = dx ;
    #ifndef KOKKOS_USING_EXPERIMENTAL_VIEW
    ASSERT_EQ( dx.tracker().ref_count() , size_t(2) );
    #endif

    {
      const_dView4 const_dx2;
      const_dx2 = const_dx;
      #ifndef KOKKOS_USING_EXPERIMENTAL_VIEW
      ASSERT_EQ( dx.tracker().ref_count() , size_t(3) );
      #endif

      const_dx2 = dy;
      #ifndef KOKKOS_USING_EXPERIMENTAL_VIEW
      ASSERT_EQ( dx.tracker().ref_count() , size_t(2) );
      #endif

      const_dView4 const_dx3(dx);
      #ifndef KOKKOS_USING_EXPERIMENTAL_VIEW
      ASSERT_EQ( dx.tracker().ref_count() , size_t(3) );
      #endif
      
      dView4_unmanaged dx4_unmanaged(dx);
      #ifndef KOKKOS_USING_EXPERIMENTAL_VIEW
      ASSERT_EQ( dx.tracker().ref_count() , size_t(3) );
      #endif
    }

    #ifndef KOKKOS_USING_EXPERIMENTAL_VIEW
    ASSERT_EQ( dx.tracker().ref_count() , size_t(2) );
    #endif


    ASSERT_FALSE( dx.is_null() );
    ASSERT_FALSE( const_dx.is_null() );
    ASSERT_FALSE( unmanaged_dx.is_null() );
    ASSERT_FALSE( unmanaged_from_ptr_dx.is_null() );
    ASSERT_FALSE( dy.is_null() );
    ASSERT_NE( dx , dy );

    ASSERT_EQ( dx.dimension_0() , unsigned(N0) );
    ASSERT_EQ( dx.dimension_1() , unsigned(N1) );
    ASSERT_EQ( dx.dimension_2() , unsigned(N2) );
    ASSERT_EQ( dx.dimension_3() , unsigned(N3) );

    ASSERT_EQ( dy.dimension_0() , unsigned(N0) );
    ASSERT_EQ( dy.dimension_1() , unsigned(N1) );
    ASSERT_EQ( dy.dimension_2() , unsigned(N2) );
    ASSERT_EQ( dy.dimension_3() , unsigned(N3) );

    ASSERT_EQ( unmanaged_from_ptr_dx.capacity(),unsigned(N0)*unsigned(N1)*unsigned(N2)*unsigned(N3) );

    hx = Kokkos::create_mirror( dx );
    hy = Kokkos::create_mirror( dy );

    // T v1 = hx() ;    // Generates compile error as intended
    // T v2 = hx(0,0) ; // Generates compile error as intended
    // hx(0,0) = v2 ;   // Generates compile error as intended

    // Testing with asynchronous deep copy with respect to device
    {
      size_t count = 0 ;
      for ( size_t ip = 0 ; ip < N0 ; ++ip ) {
      for ( size_t i1 = 0 ; i1 < hx.dimension_1() ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < hx.dimension_2() ; ++i2 ) {
      for ( size_t i3 = 0 ; i3 < hx.dimension_3() ; ++i3 ) {
        hx(ip,i1,i2,i3) = ++count ;
      }}}}

      Kokkos::deep_copy(typename hView4::execution_space(), dx , hx );
      Kokkos::deep_copy(typename hView4::execution_space(), dy , dx );
      Kokkos::deep_copy(typename hView4::execution_space(), hy , dy );

      for ( size_t ip = 0 ; ip < N0 ; ++ip ) {
      for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
      for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
        { ASSERT_EQ( hx(ip,i1,i2,i3) , hy(ip,i1,i2,i3) ); }
      }}}}

      Kokkos::deep_copy(typename hView4::execution_space(), dx , T(0) );
      Kokkos::deep_copy(typename hView4::execution_space(), hx , dx );

      for ( size_t ip = 0 ; ip < N0 ; ++ip ) {
      for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
      for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
        { ASSERT_EQ( hx(ip,i1,i2,i3) , T(0) ); }
      }}}}
    }
    // Testing with asynchronous deep copy with respect to host
    {
      size_t count = 0 ;
      for ( size_t ip = 0 ; ip < N0 ; ++ip ) {
      for ( size_t i1 = 0 ; i1 < hx.dimension_1() ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < hx.dimension_2() ; ++i2 ) {
      for ( size_t i3 = 0 ; i3 < hx.dimension_3() ; ++i3 ) {
        hx(ip,i1,i2,i3) = ++count ;
      }}}}

      Kokkos::deep_copy(typename dView4::execution_space(), dx , hx );
      Kokkos::deep_copy(typename dView4::execution_space(), dy , dx );
      Kokkos::deep_copy(typename dView4::execution_space(), hy , dy );

      for ( size_t ip = 0 ; ip < N0 ; ++ip ) {
      for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
      for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
        { ASSERT_EQ( hx(ip,i1,i2,i3) , hy(ip,i1,i2,i3) ); }
      }}}}

      Kokkos::deep_copy(typename dView4::execution_space(), dx , T(0) );
      Kokkos::deep_copy(typename dView4::execution_space(), hx , dx );

      for ( size_t ip = 0 ; ip < N0 ; ++ip ) {
      for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
      for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
        { ASSERT_EQ( hx(ip,i1,i2,i3) , T(0) ); }
      }}}}
    }

    // Testing with synchronous deep copy
    {
      size_t count = 0 ;
      for ( size_t ip = 0 ; ip < N0 ; ++ip ) {
      for ( size_t i1 = 0 ; i1 < hx.dimension_1() ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < hx.dimension_2() ; ++i2 ) {
      for ( size_t i3 = 0 ; i3 < hx.dimension_3() ; ++i3 ) {
        hx(ip,i1,i2,i3) = ++count ;
      }}}}

      Kokkos::deep_copy( dx , hx );
      Kokkos::deep_copy( dy , dx );
      Kokkos::deep_copy( hy , dy );

      for ( size_t ip = 0 ; ip < N0 ; ++ip ) {
      for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
      for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
        { ASSERT_EQ( hx(ip,i1,i2,i3) , hy(ip,i1,i2,i3) ); }
      }}}}

      Kokkos::deep_copy( dx , T(0) );
      Kokkos::deep_copy( hx , dx );

      for ( size_t ip = 0 ; ip < N0 ; ++ip ) {
      for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
      for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
      for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
        { ASSERT_EQ( hx(ip,i1,i2,i3) , T(0) ); }
      }}}}
    }
    dz = dx ; ASSERT_EQ( dx, dz); ASSERT_NE( dy, dz);
    dz = dy ; ASSERT_EQ( dy, dz); ASSERT_NE( dx, dz);

    dx = dView4();
    ASSERT_TRUE( dx.is_null() );
    ASSERT_FALSE( dy.is_null() );
    ASSERT_FALSE( dz.is_null() );
    dy = dView4();
    ASSERT_TRUE( dx.is_null() );
    ASSERT_TRUE( dy.is_null() );
    ASSERT_FALSE( dz.is_null() );
    dz = dView4();
    ASSERT_TRUE( dx.is_null() );
    ASSERT_TRUE( dy.is_null() );
    ASSERT_TRUE( dz.is_null() );
  }

  typedef T DataType[2] ;

  static void
  check_auto_conversion_to_const(
     const Kokkos::View< const DataType , device > & arg_const ,
     const Kokkos::View< DataType , device > & arg )
  {
    ASSERT_TRUE( arg_const == arg );
  }

  static void run_test_const()
  {
    typedef Kokkos::View< DataType , device > typeX ;
    typedef Kokkos::View< const DataType , device > const_typeX ;
    typedef Kokkos::View< const DataType , device , Kokkos::MemoryRandomAccess > const_typeR ;
    typeX x( "X" );
    const_typeX xc = x ;
    const_typeR xr = x ;

    ASSERT_TRUE( xc == x );
    ASSERT_TRUE( x == xc );
    ASSERT_TRUE( x.ptr_on_device() == xr.ptr_on_device() );

    // typeX xf = xc ; // setting non-const from const must not compile

    check_auto_conversion_to_const( x , x );
  }

  static void run_test_subview()
  {
    typedef Kokkos::View< const T , device > sView ;

    dView0 d0( "d0" );
    dView1 d1( "d1" , N0 );
    dView2 d2( "d2" , N0 );
    dView3 d3( "d3" , N0 );
    dView4 d4( "d4" , N0 );

    sView s0 = d0 ;
    sView s1 = Kokkos::subview( d1 , 1 );
    sView s2 = Kokkos::subview( d2 , 1 , 1 );
    sView s3 = Kokkos::subview( d3 , 1 , 1 , 1 );
    sView s4 = Kokkos::subview( d4 , 1 , 1 , 1 , 1 );
  }

  static void run_test_subview_strided()
  {
    typedef Kokkos::View< int **** , Kokkos::LayoutLeft  , host >  view_left_4 ;
    typedef Kokkos::View< int **** , Kokkos::LayoutRight , host >  view_right_4 ;
    typedef Kokkos::View< int **   , Kokkos::LayoutLeft  , host >  view_left_2 ;
    typedef Kokkos::View< int **   , Kokkos::LayoutRight , host >  view_right_2 ;

    typedef Kokkos::View< int * ,  Kokkos::LayoutStride , host >  view_stride_1 ;
    typedef Kokkos::View< int ** ,  Kokkos::LayoutStride , host >  view_stride_2 ;

    view_left_2  xl2("xl2", 100 , 200 );
    view_right_2 xr2("xr2", 100 , 200 );
    view_stride_1  yl1 = Kokkos::subview( xl2 , 0 , Kokkos::ALL() );
    view_stride_1  yl2 = Kokkos::subview( xl2 , 1 , Kokkos::ALL() );
    view_stride_1  yr1 = Kokkos::subview( xr2 , 0 , Kokkos::ALL() );
    view_stride_1  yr2 = Kokkos::subview( xr2 , 1 , Kokkos::ALL() );

    ASSERT_EQ( yl1.dimension_0() , xl2.dimension_1() );
    ASSERT_EQ( yl2.dimension_0() , xl2.dimension_1() );
    ASSERT_EQ( yr1.dimension_0() , xr2.dimension_1() );
    ASSERT_EQ( yr2.dimension_0() , xr2.dimension_1() );

    ASSERT_EQ( & yl1(0) - & xl2(0,0) , 0 );
    ASSERT_EQ( & yl2(0) - & xl2(1,0) , 0 );
    ASSERT_EQ( & yr1(0) - & xr2(0,0) , 0 );
    ASSERT_EQ( & yr2(0) - & xr2(1,0) , 0 );

    view_left_4 xl4( "xl4" , 10 , 20 , 30 , 40 );
    view_right_4 xr4( "xr4" , 10 , 20 , 30 , 40 );

    view_stride_2 yl4 = Kokkos::subview( xl4 , 1 , Kokkos::ALL() , 2 , Kokkos::ALL() );
    view_stride_2 yr4 = Kokkos::subview( xr4 , 1 , Kokkos::ALL() , 2 , Kokkos::ALL() );

    ASSERT_EQ( yl4.dimension_0() , xl4.dimension_1() );
    ASSERT_EQ( yl4.dimension_1() , xl4.dimension_3() );
    ASSERT_EQ( yr4.dimension_0() , xr4.dimension_1() );
    ASSERT_EQ( yr4.dimension_1() , xr4.dimension_3() );

    ASSERT_EQ( & yl4(4,4) - & xl4(1,4,2,4) , 0 );
    ASSERT_EQ( & yr4(4,4) - & xr4(1,4,2,4) , 0 );
  }

  static void run_test_vector()
  {
    static const unsigned Length = 1000 , Count = 8 ;

    typedef Kokkos::View< T* ,  Kokkos::LayoutLeft , host > vector_type ;
    typedef Kokkos::View< T** , Kokkos::LayoutLeft , host > multivector_type ;

    typedef Kokkos::View< T* ,  Kokkos::LayoutRight , host > vector_right_type ;
    typedef Kokkos::View< T** , Kokkos::LayoutRight , host > multivector_right_type ;

    typedef Kokkos::View< const T* , Kokkos::LayoutRight, host > const_vector_right_type ;
    typedef Kokkos::View< const T* , Kokkos::LayoutLeft , host > const_vector_type ;
    typedef Kokkos::View< const T** , Kokkos::LayoutLeft , host > const_multivector_type ;

    multivector_type mv = multivector_type( "mv" , Length , Count );
    multivector_right_type mv_right = multivector_right_type( "mv" , Length , Count );

    vector_type v1 = Kokkos::subview( mv , Kokkos::ALL() , 0 );
    vector_type v2 = Kokkos::subview( mv , Kokkos::ALL() , 1 );
    vector_type v3 = Kokkos::subview( mv , Kokkos::ALL() , 2 );

    vector_type rv1 = Kokkos::subview( mv_right , 0 , Kokkos::ALL() );
    vector_type rv2 = Kokkos::subview( mv_right , 1 , Kokkos::ALL() );
    vector_type rv3 = Kokkos::subview( mv_right , 2 , Kokkos::ALL() );

    multivector_type mv1 = Kokkos::subview( mv , std::make_pair( 1 , 998 ) ,
                                                 std::make_pair( 2 , 5 ) );

    multivector_right_type mvr1 =
      Kokkos::subview( mv_right ,
                       std::make_pair( 1 , 998 ) ,
                       std::make_pair( 2 , 5 ) );

    const_vector_type cv1 = Kokkos::subview( mv , Kokkos::ALL(), 0 );
    const_vector_type cv2 = Kokkos::subview( mv , Kokkos::ALL(), 1 );
    const_vector_type cv3 = Kokkos::subview( mv , Kokkos::ALL(), 2 );

    vector_right_type vr1 = Kokkos::subview( mv , Kokkos::ALL() , 0 );
    vector_right_type vr2 = Kokkos::subview( mv , Kokkos::ALL() , 1 );
    vector_right_type vr3 = Kokkos::subview( mv , Kokkos::ALL() , 2 );

    const_vector_right_type cvr1 = Kokkos::subview( mv , Kokkos::ALL() , 0 );
    const_vector_right_type cvr2 = Kokkos::subview( mv , Kokkos::ALL() , 1 );
    const_vector_right_type cvr3 = Kokkos::subview( mv , Kokkos::ALL() , 2 );

    ASSERT_TRUE( & v1[0] == & v1(0) );
    ASSERT_TRUE( & v1[0] == & mv(0,0) );
    ASSERT_TRUE( & v2[0] == & mv(0,1) );
    ASSERT_TRUE( & v3[0] == & mv(0,2) );

    ASSERT_TRUE( & cv1[0] == & mv(0,0) );
    ASSERT_TRUE( & cv2[0] == & mv(0,1) );
    ASSERT_TRUE( & cv3[0] == & mv(0,2) );

    ASSERT_TRUE( & vr1[0] == & mv(0,0) );
    ASSERT_TRUE( & vr2[0] == & mv(0,1) );
    ASSERT_TRUE( & vr3[0] == & mv(0,2) );

    ASSERT_TRUE( & cvr1[0] == & mv(0,0) );
    ASSERT_TRUE( & cvr2[0] == & mv(0,1) );
    ASSERT_TRUE( & cvr3[0] == & mv(0,2) );

    ASSERT_TRUE( & mv1(0,0) == & mv( 1 , 2 ) );
    ASSERT_TRUE( & mv1(1,1) == & mv( 2 , 3 ) );
    ASSERT_TRUE( & mv1(3,2) == & mv( 4 , 4 ) );
    ASSERT_TRUE( & mvr1(0,0) == & mv_right( 1 , 2 ) );
    ASSERT_TRUE( & mvr1(1,1) == & mv_right( 2 , 3 ) );
    ASSERT_TRUE( & mvr1(3,2) == & mv_right( 4 , 4 ) );

    const_vector_type c_cv1( v1 );
    typename vector_type::const_type c_cv2( v2 );
    typename const_vector_type::const_type c_ccv2( v2 );

    const_multivector_type cmv( mv );
    typename multivector_type::const_type cmvX( cmv );
    typename const_multivector_type::const_type ccmvX( cmv );
  }
};

} // namespace Test

#endif

/*--------------------------------------------------------------------------*/

