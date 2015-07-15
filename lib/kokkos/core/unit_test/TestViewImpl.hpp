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

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

/*--------------------------------------------------------------------------*/

#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )

namespace Test {

template < class Device >
void test_view_impl() {}

}

#else

/*--------------------------------------------------------------------------*/

namespace Test {

struct DummyMemorySpace
{
  typedef DummyMemorySpace memory_space ;
  typedef unsigned size_type ;
};

/*--------------------------------------------------------------------------*/

template< class Type >
struct DefineShape {
  typedef typename Kokkos::Impl::AnalyzeShape<Type>::shape type ;
};

template< class Type >
struct ExtractValueType {
  typedef typename Kokkos::Impl::AnalyzeShape<Type>::value_type type ;
};

template< class Type >
struct ArrayType { typedef Type type ; };

template < class Device >
void test_view_impl()
{
  //typedef typename Device::memory_space memory_space ; // unused

  typedef ArrayType< int[100]                >::type type_01 ;
  typedef ArrayType< int*                    >::type type_11 ;
  typedef ArrayType< int[5][6][700]          >::type type_03 ;
  typedef ArrayType< double*[8][9][900]      >::type type_14 ;
  typedef ArrayType< long**                  >::type type_22 ;
  typedef ArrayType< short **[5][6][7]       >::type type_25 ;
  typedef ArrayType< const short **[5][6][7] >::type const_type_25 ;
  typedef ArrayType< short***[5][6][7]       >::type type_36 ;
  typedef ArrayType< const short***[5][6][7] >::type const_type_36 ;

  // mfh 14 Feb 2014: With gcc 4.8.2 -Wall, this emits a warning:
  //
  // typedef ‘ok_const_25’ locally defined but not used [-Wunused-local-typedefs]
  //
  // It's unfortunate that this is the case, because the typedef is
  // being used for a compile-time check!  We deal with this by
  // declaring an instance of ok_const_25, and marking it with
  // "(void)" so that instance doesn't emit an "unused variable"
  // warning.
  //
  // typedef typename Kokkos::Impl::StaticAssertSame<
  //    typename Kokkos::Impl::AnalyzeShape<type_25>::const_type ,
  //    typename Kokkos::Impl::AnalyzeShape<const_type_25>::type
  //      > ok_const_25 ;

  typedef typename Kokkos::Impl::StaticAssertSame<
    typename Kokkos::Impl::AnalyzeShape<type_25>::const_type,
    typename Kokkos::Impl::AnalyzeShape<const_type_25>::type
      > ok_const_25 ;

  typedef typename Kokkos::Impl::StaticAssertSame<
    typename Kokkos::Impl::AnalyzeShape<type_36>::const_type,
    typename Kokkos::Impl::AnalyzeShape<const_type_36>::type
      > ok_const_36 ;
  {
    ok_const_25 thing_25 ;
    ok_const_36 thing_36 ;
    (void) thing_25 ; // silence warning for unused variable
    (void) thing_36 ; // silence warning for unused variable
  }

  ASSERT_TRUE( ( Kokkos::Impl::is_same< ExtractValueType<type_03>::type , int >::value ) );
  ASSERT_TRUE( ( Kokkos::Impl::is_same< ExtractValueType<type_14>::type , double >::value ) );
  ASSERT_TRUE( ( Kokkos::Impl::is_same< ExtractValueType<type_22>::type , long >::value ) );
  ASSERT_TRUE( ( Kokkos::Impl::is_same< ExtractValueType<type_36>::type , short >::value ) );

  ASSERT_FALSE( ( Kokkos::Impl::is_same< ExtractValueType<type_36>::type , int >::value ) );

  typedef typename DefineShape< type_01 >::type  shape_01_type ;
  typedef typename DefineShape< type_11 >::type  shape_11_type ;
  typedef typename DefineShape< type_03 >::type  shape_03_type ;
  typedef typename DefineShape< type_14 >::type  shape_14_type ;
  typedef typename DefineShape< type_22 >::type  shape_22_type ;
  typedef typename DefineShape< type_36 >::type  shape_36_type ;

  ASSERT_TRUE( ( Kokkos::Impl::StaticAssert< shape_36_type::rank == 6 >::value ) );
  ASSERT_TRUE( ( Kokkos::Impl::StaticAssert< shape_03_type::rank == 3 >::value ) );

  shape_01_type shape_01 ; shape_01_type::assign( shape_01 );
  shape_11_type shape_11 ; shape_11_type::assign( shape_11, 1000 );
  shape_03_type shape_03 ; shape_03_type::assign( shape_03 );
  shape_14_type shape_14 ; shape_14_type::assign( shape_14 , 0 );
  shape_22_type shape_22 ; shape_22_type::assign( shape_22 , 0 , 0 );
  shape_36_type shape_36 ; shape_36_type::assign( shape_36 , 10 , 20 , 30 );

  ASSERT_TRUE( shape_01.rank_dynamic == 0u );
  ASSERT_TRUE( shape_01.rank         == 1u );
  ASSERT_TRUE( shape_01.N0           == 100u );

  ASSERT_TRUE( shape_11.rank_dynamic == 1u );
  ASSERT_TRUE( shape_11.rank         == 1u );
  ASSERT_TRUE( shape_11.N0           == 1000u );

  ASSERT_TRUE( shape_03.rank_dynamic == 0u );
  ASSERT_TRUE( shape_03.rank         == 3u );
  ASSERT_TRUE( shape_03.N0           == 5u );
  ASSERT_TRUE( shape_03.N1           == 6u );
  ASSERT_TRUE( shape_03.N2           == 700u );

  ASSERT_TRUE( shape_14.rank_dynamic == 1u );
  ASSERT_TRUE( shape_14.rank         == 4u );
  ASSERT_TRUE( shape_14.N0           == 0u );
  ASSERT_TRUE( shape_14.N1           == 8u );
  ASSERT_TRUE( shape_14.N2           == 9u );
  ASSERT_TRUE( shape_14.N3           == 900u );

  ASSERT_TRUE( shape_22.rank_dynamic == 2u );
  ASSERT_TRUE( shape_22.rank         == 2u );
  ASSERT_TRUE( shape_22.N0           == 0u );
  ASSERT_TRUE( shape_22.N1           == 0u );

  ASSERT_TRUE( shape_36.rank_dynamic == 3u );
  ASSERT_TRUE( shape_36.rank         == 6u );
  ASSERT_TRUE( shape_36.N0           == 10u );
  ASSERT_TRUE( shape_36.N1           == 20u );
  ASSERT_TRUE( shape_36.N2           == 30u );
  ASSERT_TRUE( shape_36.N3           == 5u  );
  ASSERT_TRUE( shape_36.N4           == 6u  );
  ASSERT_TRUE( shape_36.N5           == 7u  );


  ASSERT_TRUE( shape_01 == shape_01 );
  ASSERT_TRUE( shape_11 == shape_11 );
  ASSERT_TRUE( shape_36 == shape_36 );
  ASSERT_TRUE( shape_01 != shape_36 );
  ASSERT_TRUE( shape_22 != shape_36 );

  //------------------------------------------------------------------------

  typedef Kokkos::Impl::ViewOffset< shape_01_type , Kokkos::LayoutLeft > shape_01_left_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_11_type , Kokkos::LayoutLeft > shape_11_left_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_03_type , Kokkos::LayoutLeft > shape_03_left_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_14_type , Kokkos::LayoutLeft > shape_14_left_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_22_type , Kokkos::LayoutLeft > shape_22_left_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_36_type , Kokkos::LayoutLeft > shape_36_left_offset ;

  typedef Kokkos::Impl::ViewOffset< shape_01_type , Kokkos::LayoutRight > shape_01_right_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_11_type , Kokkos::LayoutRight > shape_11_right_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_03_type , Kokkos::LayoutRight > shape_03_right_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_14_type , Kokkos::LayoutRight > shape_14_right_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_22_type , Kokkos::LayoutRight > shape_22_right_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_36_type , Kokkos::LayoutRight > shape_36_right_offset ;

  ASSERT_TRUE( ! shape_01_left_offset::has_padding );
  ASSERT_TRUE( ! shape_11_left_offset::has_padding );
  ASSERT_TRUE( ! shape_03_left_offset::has_padding );
  ASSERT_TRUE(   shape_14_left_offset::has_padding );
  ASSERT_TRUE(   shape_22_left_offset::has_padding );
  ASSERT_TRUE(   shape_36_left_offset::has_padding );

  ASSERT_TRUE( ! shape_01_right_offset::has_padding );
  ASSERT_TRUE( ! shape_11_right_offset::has_padding );
  ASSERT_TRUE( ! shape_03_right_offset::has_padding );
  ASSERT_TRUE( ! shape_14_right_offset::has_padding );
  ASSERT_TRUE(   shape_22_right_offset::has_padding );
  ASSERT_TRUE(   shape_36_right_offset::has_padding );

  //------------------------------------------------------------------------

  typedef Kokkos::Impl::ViewOffset< shape_01_type , Kokkos::LayoutStride > shape_01_stride_offset ;
  typedef Kokkos::Impl::ViewOffset< shape_36_type , Kokkos::LayoutStride > shape_36_stride_offset ;

  {
    shape_01_stride_offset stride_offset_01 ;

    stride_offset_01.assign( 1, stride_offset_01.N0, 0,0,0,0,0,0,0 );

    ASSERT_EQ( int(stride_offset_01.S[0]) , int(1) );
    ASSERT_EQ( int(stride_offset_01.S[1]) , int(stride_offset_01.N0) );
  }

  {
    shape_36_stride_offset stride_offset_36 ;

    size_t str[7] ;
    str[5] = 1 ;
    str[4] = str[5] * stride_offset_36.N5 ;
    str[3] = str[4] * stride_offset_36.N4 ;
    str[2] = str[3] * stride_offset_36.N3 ;
    str[1] = str[2] * 100 ;
    str[0] = str[1] * 200 ;
    str[6] = str[0] * 300 ;

    stride_offset_36.assign( str[0] , str[1] , str[2] , str[3] , str[4] , str[5] , str[6] , 0 , 0 );

    ASSERT_EQ( size_t(stride_offset_36.S[6]) , size_t(str[6]) );
    ASSERT_EQ( size_t(stride_offset_36.N2)   , size_t(100) );
    ASSERT_EQ( size_t(stride_offset_36.N1)   , size_t(200) );
    ASSERT_EQ( size_t(stride_offset_36.N0)   , size_t(300) );
  }

  //------------------------------------------------------------------------

  {
    const int rank = 6 ;
    const int order[] = { 5 , 3 , 1 , 0 , 2 , 4 };
    const unsigned dim[] = { 2 , 3 , 5 , 7 , 11 , 13 };
    Kokkos::LayoutStride stride_6 = Kokkos::LayoutStride::order_dimensions( rank , order , dim );
    size_t n = 1 ;
    for ( int i = 0 ; i < rank ; ++i ) {
      ASSERT_EQ( size_t(dim[i]) , size_t( stride_6.dimension[i] ) );
      ASSERT_EQ( size_t(n) , size_t( stride_6.stride[ order[i] ] ) );
      n *= dim[order[i]] ;
    }
  }

  //------------------------------------------------------------------------
}

} /* namespace Test */

#endif

/*--------------------------------------------------------------------------*/

