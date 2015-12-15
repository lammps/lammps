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

namespace Test {

template< class ExecSpace >
void test_view_mapping()
{
  typedef Kokkos::Experimental::Impl::ViewDimension<>  dim_0 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<2> dim_s2 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<2,3> dim_s2_s3 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<2,3,4> dim_s2_s3_s4 ;

  typedef Kokkos::Experimental::Impl::ViewDimension<0> dim_s0 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<0,3> dim_s0_s3 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<0,3,4> dim_s0_s3_s4 ;

  typedef Kokkos::Experimental::Impl::ViewDimension<0,0> dim_s0_s0 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<0,0,4> dim_s0_s0_s4 ;

  typedef Kokkos::Experimental::Impl::ViewDimension<0,0,0> dim_s0_s0_s0 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<0,0,0,0> dim_s0_s0_s0_s0 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<0,0,0,0,0> dim_s0_s0_s0_s0_s0 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<0,0,0,0,0,0> dim_s0_s0_s0_s0_s0_s0 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<0,0,0,0,0,0,0> dim_s0_s0_s0_s0_s0_s0_s0 ;
  typedef Kokkos::Experimental::Impl::ViewDimension<0,0,0,0,0,0,0,0> dim_s0_s0_s0_s0_s0_s0_s0_s0 ;

  // Fully static dimensions should not be larger than an int
  ASSERT_LE( sizeof(dim_0) , sizeof(int) );
  ASSERT_LE( sizeof(dim_s2) , sizeof(int) );
  ASSERT_LE( sizeof(dim_s2_s3) , sizeof(int) );
  ASSERT_LE( sizeof(dim_s2_s3_s4) , sizeof(int) );

  // Rank 1 is size_t
  ASSERT_EQ( sizeof(dim_s0) , sizeof(size_t) );
  ASSERT_EQ( sizeof(dim_s0_s3) , sizeof(size_t) );
  ASSERT_EQ( sizeof(dim_s0_s3_s4) , sizeof(size_t) );

  // Allow for padding
  ASSERT_LE( sizeof(dim_s0_s0) , 2 * sizeof(size_t) );
  ASSERT_LE( sizeof(dim_s0_s0_s4) , 2 * sizeof(size_t) );

  ASSERT_LE( sizeof(dim_s0_s0_s0) , 4 * sizeof(size_t) );
  ASSERT_EQ( sizeof(dim_s0_s0_s0_s0) , 4 * sizeof(unsigned) );
  ASSERT_LE( sizeof(dim_s0_s0_s0_s0_s0) , 6 * sizeof(unsigned) );
  ASSERT_EQ( sizeof(dim_s0_s0_s0_s0_s0_s0) , 6 * sizeof(unsigned) );
  ASSERT_LE( sizeof(dim_s0_s0_s0_s0_s0_s0_s0) , 8 * sizeof(unsigned) );
  ASSERT_EQ( sizeof(dim_s0_s0_s0_s0_s0_s0_s0_s0) , 8 * sizeof(unsigned) );

  ASSERT_EQ( int(dim_0::rank) , int(0) );
  ASSERT_EQ( int(dim_0::rank_dynamic) , int(0) );

  ASSERT_EQ( int(dim_s2::rank) , int(1) );
  ASSERT_EQ( int(dim_s2::rank_dynamic) , int(0) );

  ASSERT_EQ( int(dim_s2_s3::rank) , int(2) );
  ASSERT_EQ( int(dim_s2_s3::rank_dynamic) , int(0) );

  ASSERT_EQ( int(dim_s2_s3_s4::rank) , int(3) );
  ASSERT_EQ( int(dim_s2_s3_s4::rank_dynamic) , int(0) );

  ASSERT_EQ( int(dim_s0::rank) , int(1) );
  ASSERT_EQ( int(dim_s0::rank_dynamic) , int(1) );

  ASSERT_EQ( int(dim_s0_s3::rank) , int(2) );
  ASSERT_EQ( int(dim_s0_s3::rank_dynamic) , int(1) );

  ASSERT_EQ( int(dim_s0_s3_s4::rank) , int(3) );
  ASSERT_EQ( int(dim_s0_s3_s4::rank_dynamic) , int(1) );

  ASSERT_EQ( int(dim_s0_s0_s4::rank) , int(3) );
  ASSERT_EQ( int(dim_s0_s0_s4::rank_dynamic) , int(2) );

  ASSERT_EQ( int(dim_s0_s0_s0::rank) , int(3) );
  ASSERT_EQ( int(dim_s0_s0_s0::rank_dynamic) , int(3) );

  ASSERT_EQ( int(dim_s0_s0_s0_s0::rank) , int(4) );
  ASSERT_EQ( int(dim_s0_s0_s0_s0::rank_dynamic) , int(4) );

  ASSERT_EQ( int(dim_s0_s0_s0_s0_s0::rank) , int(5) );
  ASSERT_EQ( int(dim_s0_s0_s0_s0_s0::rank_dynamic) , int(5) );

  ASSERT_EQ( int(dim_s0_s0_s0_s0_s0_s0::rank) , int(6) );
  ASSERT_EQ( int(dim_s0_s0_s0_s0_s0_s0::rank_dynamic) , int(6) );

  ASSERT_EQ( int(dim_s0_s0_s0_s0_s0_s0_s0::rank) , int(7) );
  ASSERT_EQ( int(dim_s0_s0_s0_s0_s0_s0_s0::rank_dynamic) , int(7) );

  ASSERT_EQ( int(dim_s0_s0_s0_s0_s0_s0_s0_s0::rank) , int(8) );
  ASSERT_EQ( int(dim_s0_s0_s0_s0_s0_s0_s0_s0::rank_dynamic) , int(8) );

  dim_s0          d1( 2, 3, 4, 5, 6, 7, 8, 9 ); 
  dim_s0_s0       d2( 2, 3, 4, 5, 6, 7, 8, 9 );
  dim_s0_s0_s0    d3( 2, 3, 4, 5, 6, 7, 8, 9 );
  dim_s0_s0_s0_s0 d4( 2, 3, 4, 5, 6, 7, 8, 9 );

  ASSERT_EQ( d1.N0 , 2 );
  ASSERT_EQ( d2.N0 , 2 );
  ASSERT_EQ( d3.N0 , 2 );
  ASSERT_EQ( d4.N0 , 2 );

  ASSERT_EQ( d1.N1 , 1 );
  ASSERT_EQ( d2.N1 , 3 );
  ASSERT_EQ( d3.N1 , 3 );
  ASSERT_EQ( d4.N1 , 3 );

  ASSERT_EQ( d1.N2 , 1 );
  ASSERT_EQ( d2.N2 , 1 );
  ASSERT_EQ( d3.N2 , 4 );
  ASSERT_EQ( d4.N2 , 4 );

  ASSERT_EQ( d1.N3 , 1 );
  ASSERT_EQ( d2.N3 , 1 );
  ASSERT_EQ( d3.N3 , 1 );
  ASSERT_EQ( d4.N3 , 5 );

  //----------------------------------------

  typedef Kokkos::Experimental::Impl::ViewOffset< dim_s0_s0_s0 , Kokkos::LayoutStride >  stride_s0_s0_s0 ;

  //----------------------------------------
  // Static dimension
  {
    typedef Kokkos::Experimental::Impl::ViewOffset< dim_s2_s3_s4 , Kokkos::LayoutLeft > left_s2_s3_s4 ;

    ASSERT_EQ( sizeof(left_s2_s3_s4) , sizeof(dim_s2_s3_s4) );

    left_s2_s3_s4 off3 ;

    stride_s0_s0_s0  stride3( off3 );

    ASSERT_EQ( off3.stride_0() , 1 );
    ASSERT_EQ( off3.stride_1() , 2 );
    ASSERT_EQ( off3.stride_2() , 6 );
    ASSERT_EQ( off3.span() , 24 );

    ASSERT_EQ( off3.stride_0() , stride3.stride_0() );
    ASSERT_EQ( off3.stride_1() , stride3.stride_1() );
    ASSERT_EQ( off3.stride_2() , stride3.stride_2() );
    ASSERT_EQ( off3.span() , stride3.span() );

    int offset = 0 ;

    for ( int k = 0 ; k < 4 ; ++k ){
    for ( int j = 0 ; j < 3 ; ++j ){
    for ( int i = 0 ; i < 2 ; ++i , ++offset ){
      ASSERT_EQ( off3(i,j,k) , offset );
      ASSERT_EQ( stride3(i,j,k) , off3(i,j,k) );
    }}}
  }

  //----------------------------------------
  // Small dimension is unpadded
  {
    typedef Kokkos::Experimental::Impl::ViewOffset< dim_s0_s0_s4 , Kokkos::LayoutLeft > left_s0_s0_s4 ;

    left_s0_s0_s4 dyn_off3( std::integral_constant<unsigned,sizeof(int)>(), 2, 3, 0, 0, 0, 0, 0, 0 );

    stride_s0_s0_s0  stride3( dyn_off3 );

    ASSERT_EQ( dyn_off3.m_dim.rank , 3 );
    ASSERT_EQ( dyn_off3.m_dim.N0 , 2 );
    ASSERT_EQ( dyn_off3.m_dim.N1 , 3 );
    ASSERT_EQ( dyn_off3.m_dim.N2 , 4 );
    ASSERT_EQ( dyn_off3.m_dim.N3 , 1 );
    ASSERT_EQ( dyn_off3.size() , 2 * 3 * 4 );

    ASSERT_EQ( stride3.m_dim.rank , 3 );
    ASSERT_EQ( stride3.m_dim.N0 , 2 );
    ASSERT_EQ( stride3.m_dim.N1 , 3 );
    ASSERT_EQ( stride3.m_dim.N2 , 4 );
    ASSERT_EQ( stride3.m_dim.N3 , 1 );
    ASSERT_EQ( stride3.size() , 2 * 3 * 4 );

    int offset = 0 ;

    for ( int k = 0 ; k < 4 ; ++k ){
    for ( int j = 0 ; j < 3 ; ++j ){
    for ( int i = 0 ; i < 2 ; ++i , ++offset ){
      ASSERT_EQ( offset , dyn_off3(i,j,k) );
      ASSERT_EQ( stride3(i,j,k) , dyn_off3(i,j,k) );
    }}}

    ASSERT_EQ( dyn_off3.span() , offset );
    ASSERT_EQ( stride3.span() , dyn_off3.span() );
  }

  // Large dimension is likely padded
  {
    constexpr int N0 = 2000 ;
    constexpr int N1 = 300 ;

    typedef Kokkos::Experimental::Impl::ViewOffset< dim_s0_s0_s4 , Kokkos::LayoutLeft > left_s0_s0_s4 ;

    left_s0_s0_s4 dyn_off3( std::integral_constant<unsigned,sizeof(int)>(), N0, N1, 0, 0, 0, 0, 0, 0 );

    stride_s0_s0_s0  stride3( dyn_off3 );

    ASSERT_EQ( dyn_off3.m_dim.rank , 3 );
    ASSERT_EQ( dyn_off3.m_dim.N0 , N0 );
    ASSERT_EQ( dyn_off3.m_dim.N1 , N1 );
    ASSERT_EQ( dyn_off3.m_dim.N2 , 4 );
    ASSERT_EQ( dyn_off3.m_dim.N3 , 1 );
    ASSERT_EQ( dyn_off3.size() , N0 * N1 * 4 );

    ASSERT_EQ( stride3.m_dim.rank , 3 );
    ASSERT_EQ( stride3.m_dim.N0 , N0 );
    ASSERT_EQ( stride3.m_dim.N1 , N1 );
    ASSERT_EQ( stride3.m_dim.N2 , 4 );
    ASSERT_EQ( stride3.m_dim.N3 , 1 );
    ASSERT_EQ( stride3.size() , N0 * N1 * 4 );
    ASSERT_EQ( stride3.span() , dyn_off3.span() );

    int offset = 0 ;

    for ( int k = 0 ; k < 4 ; ++k ){
    for ( int j = 0 ; j < N1 ; ++j ){
    for ( int i = 0 ; i < N0 ; ++i ){
      ASSERT_LE( offset , dyn_off3(i,j,k) );
      ASSERT_EQ( stride3(i,j,k) , dyn_off3(i,j,k) );
      offset = dyn_off3(i,j,k) + 1 ;
    }}}

    ASSERT_LE( offset , dyn_off3.span() );
  }

  //----------------------------------------
  // Static dimension
  {
    typedef Kokkos::Experimental::Impl::ViewOffset< dim_s2_s3_s4 , Kokkos::LayoutRight > right_s2_s3_s4 ;

    ASSERT_EQ( sizeof(right_s2_s3_s4) , sizeof(dim_s2_s3_s4) );

    right_s2_s3_s4 off3 ;

    stride_s0_s0_s0  stride3( off3 );

    ASSERT_EQ( off3.stride_0() , 12 );
    ASSERT_EQ( off3.stride_1() , 4 );
    ASSERT_EQ( off3.stride_2() , 1 );

    ASSERT_EQ( off3.dimension_0() , stride3.dimension_0() );
    ASSERT_EQ( off3.dimension_1() , stride3.dimension_1() );
    ASSERT_EQ( off3.dimension_2() , stride3.dimension_2() );
    ASSERT_EQ( off3.stride_0() , stride3.stride_0() );
    ASSERT_EQ( off3.stride_1() , stride3.stride_1() );
    ASSERT_EQ( off3.stride_2() , stride3.stride_2() );
    ASSERT_EQ( off3.span() , stride3.span() );

    int offset = 0 ;

    for ( int i = 0 ; i < 2 ; ++i ){
    for ( int j = 0 ; j < 3 ; ++j ){
    for ( int k = 0 ; k < 4 ; ++k , ++offset ){
      ASSERT_EQ( off3(i,j,k) , offset );
      ASSERT_EQ( off3(i,j,k) , stride3(i,j,k) );
    }}}

    ASSERT_EQ( off3.span() , offset );
  }

  //----------------------------------------
  // Small dimension is unpadded
  {
    typedef Kokkos::Experimental::Impl::ViewOffset< dim_s0_s0_s4 , Kokkos::LayoutRight > right_s0_s0_s4 ;

    right_s0_s0_s4 dyn_off3( std::integral_constant<unsigned,sizeof(int)>(), 2, 3, 0, 0, 0, 0, 0, 0 );

    stride_s0_s0_s0  stride3( dyn_off3 );

    ASSERT_EQ( dyn_off3.m_dim.rank , 3 );
    ASSERT_EQ( dyn_off3.m_dim.N0 , 2 );
    ASSERT_EQ( dyn_off3.m_dim.N1 , 3 );
    ASSERT_EQ( dyn_off3.m_dim.N2 , 4 );
    ASSERT_EQ( dyn_off3.m_dim.N3 , 1 );
    ASSERT_EQ( dyn_off3.size() , 2 * 3 * 4 );

    ASSERT_EQ( dyn_off3.dimension_0() , stride3.dimension_0() );
    ASSERT_EQ( dyn_off3.dimension_1() , stride3.dimension_1() );
    ASSERT_EQ( dyn_off3.dimension_2() , stride3.dimension_2() );
    ASSERT_EQ( dyn_off3.stride_0() , stride3.stride_0() );
    ASSERT_EQ( dyn_off3.stride_1() , stride3.stride_1() );
    ASSERT_EQ( dyn_off3.stride_2() , stride3.stride_2() );
    ASSERT_EQ( dyn_off3.span() , stride3.span() );

    int offset = 0 ;

    for ( int i = 0 ; i < 2 ; ++i ){
    for ( int j = 0 ; j < 3 ; ++j ){
    for ( int k = 0 ; k < 4 ; ++k , ++offset ){
      ASSERT_EQ( offset , dyn_off3(i,j,k) );
      ASSERT_EQ( dyn_off3(i,j,k) , stride3(i,j,k) );
    }}}

    ASSERT_EQ( dyn_off3.span() , offset );
  }

  // Large dimension is likely padded
  {
    constexpr int N0 = 2000 ;
    constexpr int N1 = 300 ;

    typedef Kokkos::Experimental::Impl::ViewOffset< dim_s0_s0_s4 , Kokkos::LayoutRight > right_s0_s0_s4 ;

    right_s0_s0_s4 dyn_off3( std::integral_constant<unsigned,sizeof(int)>(), N0, N1, 0, 0, 0, 0, 0, 0 );

    stride_s0_s0_s0  stride3( dyn_off3 );

    ASSERT_EQ( dyn_off3.m_dim.rank , 3 );
    ASSERT_EQ( dyn_off3.m_dim.N0 , N0 );
    ASSERT_EQ( dyn_off3.m_dim.N1 , N1 );
    ASSERT_EQ( dyn_off3.m_dim.N2 , 4 );
    ASSERT_EQ( dyn_off3.m_dim.N3 , 1 );
    ASSERT_EQ( dyn_off3.size() , N0 * N1 * 4 );

    ASSERT_EQ( dyn_off3.dimension_0() , stride3.dimension_0() );
    ASSERT_EQ( dyn_off3.dimension_1() , stride3.dimension_1() );
    ASSERT_EQ( dyn_off3.dimension_2() , stride3.dimension_2() );
    ASSERT_EQ( dyn_off3.stride_0() , stride3.stride_0() );
    ASSERT_EQ( dyn_off3.stride_1() , stride3.stride_1() );
    ASSERT_EQ( dyn_off3.stride_2() , stride3.stride_2() );
    ASSERT_EQ( dyn_off3.span() , stride3.span() );

    int offset = 0 ;

    for ( int i = 0 ; i < N0 ; ++i ){
    for ( int j = 0 ; j < N1 ; ++j ){
    for ( int k = 0 ; k < 4 ; ++k ){
      ASSERT_LE( offset , dyn_off3(i,j,k) );
      ASSERT_EQ( dyn_off3(i,j,k) , stride3(i,j,k) );
      offset = dyn_off3(i,j,k) + 1 ;
    }}}

    ASSERT_LE( offset , dyn_off3.span() );
  }

  //----------------------------------------
  // Subview
  {
    // Mapping rank 4 to rank 3
    typedef Kokkos::Experimental::Impl::SubviewExtents<4,3> SubviewExtents ;

    constexpr int N0 = 1000 ;
    constexpr int N1 = 2000 ;
    constexpr int N2 = 3000 ;
    constexpr int N3 = 4000 ;

    Kokkos::Experimental::Impl::ViewDimension<N0,N1,N2,N3> dim ;

    SubviewExtents tmp( dim
                      , N0 / 2
                      , Kokkos::Experimental::ALL
                      , std::pair<int,int>( N2 / 4 , 10 + N2 / 4 )
                      , Kokkos::pair<int,int>( N3 / 4 , 20 + N3 / 4 )
                      );

    ASSERT_EQ( tmp.domain_offset(0) , N0 / 2 );
    ASSERT_EQ( tmp.domain_offset(1) , 0 );
    ASSERT_EQ( tmp.domain_offset(2) , N2 / 4 );
    ASSERT_EQ( tmp.domain_offset(3) , N3 / 4 );

    ASSERT_EQ( tmp.range_index(0) , 1 );
    ASSERT_EQ( tmp.range_index(1) , 2 );
    ASSERT_EQ( tmp.range_index(2) , 3 );

    ASSERT_EQ( tmp.range_extent(0) , N1 );
    ASSERT_EQ( tmp.range_extent(1) , 10 );
    ASSERT_EQ( tmp.range_extent(2) , 20 );
  }
  //----------------------------------------
  {
    constexpr int N0 = 2000 ;
    constexpr int N1 = 300 ;

    constexpr int sub_N0 = 1000 ;
    constexpr int sub_N1 = 200 ;
    constexpr int sub_N2 = 4 ;

    typedef Kokkos::Experimental::Impl::ViewOffset< dim_s0_s0_s4 , Kokkos::LayoutLeft > left_s0_s0_s4 ;

    left_s0_s0_s4 dyn_off3( std::integral_constant<unsigned,sizeof(int)>(), N0, N1, 0, 0, 0, 0, 0, 0 );

    Kokkos::Experimental::Impl::SubviewExtents< 3 , 3 >
      sub( dyn_off3.m_dim
         , Kokkos::pair<int,int>(0,sub_N0)
         , Kokkos::pair<int,int>(0,sub_N1)
         , Kokkos::pair<int,int>(0,sub_N2)
         );

    stride_s0_s0_s0  stride3( dyn_off3 , sub );

    ASSERT_EQ( stride3.dimension_0() , sub_N0 );
    ASSERT_EQ( stride3.dimension_1() , sub_N1 );
    ASSERT_EQ( stride3.dimension_2() , sub_N2 );
    ASSERT_EQ( stride3.size() , sub_N0 * sub_N1 * sub_N2 );

    ASSERT_EQ( dyn_off3.stride_0() , stride3.stride_0() );
    ASSERT_EQ( dyn_off3.stride_1() , stride3.stride_1() );
    ASSERT_EQ( dyn_off3.stride_2() , stride3.stride_2() );
    ASSERT_GE( dyn_off3.span()   , stride3.span() );

    for ( int k = 0 ; k < sub_N2 ; ++k ){
    for ( int j = 0 ; j < sub_N1 ; ++j ){
    for ( int i = 0 ; i < sub_N0 ; ++i ){
      ASSERT_EQ( stride3(i,j,k) , dyn_off3(i,j,k) );
    }}}
  }

  {
    constexpr int N0 = 2000 ;
    constexpr int N1 = 300 ;

    constexpr int sub_N0 = 1000 ;
    constexpr int sub_N1 = 200 ;
    constexpr int sub_N2 = 4 ;

    typedef Kokkos::Experimental::Impl::ViewOffset< dim_s0_s0_s4 , Kokkos::LayoutRight > right_s0_s0_s4 ;

    right_s0_s0_s4 dyn_off3( std::integral_constant<unsigned,sizeof(int)>(), N0, N1, 0, 0, 0, 0, 0, 0 );

    Kokkos::Experimental::Impl::SubviewExtents< 3 , 3 >
      sub( dyn_off3.m_dim
         , Kokkos::pair<int,int>(0,sub_N0)
         , Kokkos::pair<int,int>(0,sub_N1)
         , Kokkos::pair<int,int>(0,sub_N2)
         );

    stride_s0_s0_s0  stride3( dyn_off3 , sub );

    ASSERT_EQ( stride3.dimension_0() , sub_N0 );
    ASSERT_EQ( stride3.dimension_1() , sub_N1 );
    ASSERT_EQ( stride3.dimension_2() , sub_N2 );
    ASSERT_EQ( stride3.size() , sub_N0 * sub_N1 * sub_N2 );

    ASSERT_EQ( dyn_off3.stride_0() , stride3.stride_0() );
    ASSERT_EQ( dyn_off3.stride_1() , stride3.stride_1() );
    ASSERT_EQ( dyn_off3.stride_2() , stride3.stride_2() );
    ASSERT_GE( dyn_off3.span()   , stride3.span() );

    for ( int i = 0 ; i < sub_N0 ; ++i ){
    for ( int j = 0 ; j < sub_N1 ; ++j ){
    for ( int k = 0 ; k < sub_N2 ; ++k ){
      ASSERT_EQ( stride3(i,j,k) , dyn_off3(i,j,k) );
    }}}
  }

  //----------------------------------------
  // view data analysis
  {
    using namespace Kokkos::Experimental::Impl ;
    static_assert( rank_dynamic<>::value == 0 , "" );
    static_assert( rank_dynamic<1>::value == 0 , "" );
    static_assert( rank_dynamic<0>::value == 1 , "" );
    static_assert( rank_dynamic<0,1>::value == 1 , "" );
    static_assert( rank_dynamic<0,0,1>::value == 2 , "" );
  }

  {
    using namespace Kokkos::Experimental::Impl ;

    typedef ViewArrayAnalysis< int[] >                 a_int_r1 ;
    typedef ViewArrayAnalysis< int**[4][5][6] >        a_int_r5 ;
    typedef ViewArrayAnalysis< const int[] >           a_const_int_r1 ;
    typedef ViewArrayAnalysis< const int**[4][5][6] >  a_const_int_r5 ;

    static_assert( a_int_r1::dimension::rank == 1 , "" );
    static_assert( a_int_r1::dimension::rank_dynamic == 1 , "" );
    static_assert( std::is_same< typename a_int_r1::dimension , ViewDimension<0> >::value , "" );
    static_assert( std::is_same< typename a_int_r1::non_const_value_type , int >::value , "" );

    static_assert( a_const_int_r1::dimension::rank == 1 , "" );
    static_assert( a_const_int_r1::dimension::rank_dynamic == 1 , "" );
    static_assert( std::is_same< typename a_const_int_r1::dimension , ViewDimension<0> >::value , "" );
    static_assert( std::is_same< typename a_const_int_r1::non_const_value_type , int >::value , "" );

    static_assert( a_const_int_r5::dimension::rank == 5 , "" );
    static_assert( a_const_int_r5::dimension::rank_dynamic == 2 , "" );

    static_assert( std::is_same< typename a_const_int_r5::dimension , ViewDimension<0,0,4,5,6> >::value , "" );

    static_assert( std::is_same< typename a_const_int_r5::non_const_value_type , int >::value , "" );

    static_assert( a_int_r5::dimension::rank == 5 , "" );
    static_assert( a_int_r5::dimension::rank_dynamic == 2 , "" );
    static_assert( std::is_same< typename a_int_r5::dimension , ViewDimension<0,0,4,5,6> >::value , "" );
    static_assert( std::is_same< typename a_int_r5::non_const_value_type , int >::value , "" );
  }

  {
    using namespace Kokkos::Experimental::Impl ;

    typedef int t_i4[4] ;

    // Dimensions of t_i4 are appended to the multdimensional array.
    typedef ViewArrayAnalysis< t_i4 ***[3] > a_int_r5 ;

    static_assert( a_int_r5::dimension::rank == 5 , "" );
    static_assert( a_int_r5::dimension::rank_dynamic == 3 , "" );
    static_assert( a_int_r5::dimension::ArgN0 == 0 , "" );
    static_assert( a_int_r5::dimension::ArgN1 == 0 , "" );
    static_assert( a_int_r5::dimension::ArgN2 == 0 , "" );
    static_assert( a_int_r5::dimension::ArgN3 == 3 , "" );
    static_assert( a_int_r5::dimension::ArgN4 == 4 , "" );
    static_assert( std::is_same< typename a_int_r5::non_const_value_type , int >::value , "" );
  }

  {
    using namespace Kokkos::Experimental::Impl ;

    typedef ViewDataAnalysis< const int[] , void >  a_const_int_r1 ;

    static_assert( std::is_same< typename a_const_int_r1::specialize , void >::value , "" );
    static_assert( std::is_same< typename a_const_int_r1::dimension , Kokkos::Experimental::Impl::ViewDimension<0> >::value , "" );

    static_assert( std::is_same< typename a_const_int_r1::type , const int * >::value , "" );
    static_assert( std::is_same< typename a_const_int_r1::value_type , const int >::value , "" );

    static_assert( std::is_same< typename a_const_int_r1::array_scalar_type , const int * >::value , "" );
    static_assert( std::is_same< typename a_const_int_r1::const_type , const int * >::value , "" );
    static_assert( std::is_same< typename a_const_int_r1::const_value_type , const int >::value , "" );
    static_assert( std::is_same< typename a_const_int_r1::const_array_scalar_type , const int * >::value , "" );
    static_assert( std::is_same< typename a_const_int_r1::non_const_type , int * >::value , "" );
    static_assert( std::is_same< typename a_const_int_r1::non_const_value_type , int >::value , "" );

    typedef ViewDataAnalysis< const int**[4] , void >  a_const_int_r3 ;

    static_assert( std::is_same< typename a_const_int_r3::specialize , void >::value , "" );

    static_assert( std::is_same< typename a_const_int_r3::dimension , Kokkos::Experimental::Impl::ViewDimension<0,0,4> >::value , "" );

    static_assert( std::is_same< typename a_const_int_r3::type , const int**[4] >::value , "" );
    static_assert( std::is_same< typename a_const_int_r3::value_type , const int >::value , "" );
    static_assert( std::is_same< typename a_const_int_r3::array_scalar_type , const int**[4] >::value , "" );
    static_assert( std::is_same< typename a_const_int_r3::const_type , const int**[4] >::value , "" );
    static_assert( std::is_same< typename a_const_int_r3::const_value_type , const int >::value , "" );
    static_assert( std::is_same< typename a_const_int_r3::const_array_scalar_type , const int**[4] >::value , "" );
    static_assert( std::is_same< typename a_const_int_r3::non_const_type , int**[4] >::value , "" );
    static_assert( std::is_same< typename a_const_int_r3::non_const_value_type , int >::value , "" );
    static_assert( std::is_same< typename a_const_int_r3::non_const_array_scalar_type , int**[4] >::value , "" );


    // std::cout << "typeid(const int**[4]).name() = " << typeid(const int**[4]).name() << std::endl ;
  }

  //----------------------------------------

  {
    constexpr int N = 10 ;

    typedef Kokkos::Experimental::View<int*,ExecSpace>        T ;
    typedef Kokkos::Experimental::View<const int*,ExecSpace>  C ;

    int data[N] ;

    T vr1(data,N); // view of non-const
    C cr1(vr1);    // view of const from view of non-const
    C cr2( (const int *) data , N );

    // Generate static_assert error:
    // T tmp( cr1 );

    ASSERT_EQ( vr1.span() , N );
    ASSERT_EQ( cr1.span() , N );
    ASSERT_EQ( vr1.data() , & data[0] );
    ASSERT_EQ( cr1.data() , & data[0] );

    ASSERT_TRUE( ( std::is_same< typename T::data_type           , int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::const_data_type     , const int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::non_const_data_type , int* >::value ) );

    ASSERT_TRUE( ( std::is_same< typename T::array_scalar_type           , int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::const_array_scalar_type     , const int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::non_const_array_scalar_type , int* >::value ) );

    ASSERT_TRUE( ( std::is_same< typename T::value_type           , int >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::const_value_type     , const int >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::non_const_value_type , int >::value ) );

    ASSERT_TRUE( ( std::is_same< typename T::memory_space , typename ExecSpace::memory_space >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::reference_type , int & >::value ) );

    ASSERT_EQ( T::Rank , 1 );

    ASSERT_TRUE( ( std::is_same< typename C::data_type           , const int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename C::const_data_type     , const int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename C::non_const_data_type , int* >::value ) );

    ASSERT_TRUE( ( std::is_same< typename C::array_scalar_type           , const int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename C::const_array_scalar_type     , const int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename C::non_const_array_scalar_type , int* >::value ) );

    ASSERT_TRUE( ( std::is_same< typename C::value_type           , const int >::value ) );
    ASSERT_TRUE( ( std::is_same< typename C::const_value_type     , const int >::value ) );
    ASSERT_TRUE( ( std::is_same< typename C::non_const_value_type , int >::value ) );

    ASSERT_TRUE( ( std::is_same< typename C::memory_space , typename ExecSpace::memory_space >::value ) );
    ASSERT_TRUE( ( std::is_same< typename C::reference_type , const int & >::value ) );

    ASSERT_EQ( C::Rank , 1 );

    ASSERT_EQ( vr1.dimension_0() , N );

    if ( Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< typename ExecSpace::memory_space , Kokkos::HostSpace >::value ) {
      for ( int i = 0 ; i < N ; ++i ) data[i] = i + 1 ;
      for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( vr1[i] , i + 1 );
      for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( cr1[i] , i + 1 );

      {
        T tmp( vr1 );
        for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( tmp[i] , i + 1 );
        for ( int i = 0 ; i < N ; ++i ) vr1(i) = i + 2 ;
        for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( tmp[i] , i + 2 );
      }

      for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( vr1[i] , i + 2 );
    }
  }

  {
    constexpr int N = 10 ;
    typedef Kokkos::Experimental::View<int*,ExecSpace>        T ;
    typedef Kokkos::Experimental::View<const int*,ExecSpace>  C ;

    T vr1("vr1",N);
    C cr1(vr1);

    ASSERT_TRUE( ( std::is_same< typename T::data_type           , int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::const_data_type     , const int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::non_const_data_type , int* >::value ) );

    ASSERT_TRUE( ( std::is_same< typename T::array_scalar_type           , int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::const_array_scalar_type     , const int* >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::non_const_array_scalar_type , int* >::value ) );

    ASSERT_TRUE( ( std::is_same< typename T::value_type           , int >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::const_value_type     , const int >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::non_const_value_type , int >::value ) );

    ASSERT_TRUE( ( std::is_same< typename T::memory_space , typename ExecSpace::memory_space >::value ) );
    ASSERT_TRUE( ( std::is_same< typename T::reference_type , int & >::value ) );
    ASSERT_EQ( T::Rank , 1 );
 
    ASSERT_EQ( vr1.dimension_0() , N );

    if ( Kokkos::Impl::VerifyExecutionCanAccessMemorySpace< typename ExecSpace::memory_space , Kokkos::HostSpace >::value ) {
      for ( int i = 0 ; i < N ; ++i ) vr1(i) = i + 1 ;
      for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( vr1[i] , i + 1 );
      for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( cr1[i] , i + 1 );

      {
        T tmp( vr1 );
        for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( tmp[i] , i + 1 );
        for ( int i = 0 ; i < N ; ++i ) vr1(i) = i + 2 ;
        for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( tmp[i] , i + 2 );
      }

      for ( int i = 0 ; i < N ; ++i ) ASSERT_EQ( vr1[i] , i + 2 );
    }
  }

  // Testing using space instance for allocation.
  // The execution space of the memory space must be available for view data initialization

  if ( std::is_same< ExecSpace , typename ExecSpace::memory_space::execution_space >::value ) {

    using namespace Kokkos::Experimental ;

    typedef typename ExecSpace::memory_space  memory_space ;
    typedef View<int*,memory_space>           V ;

    constexpr int N = 10 ;

    memory_space mem_space ;

    V v( "v" , N );
    V va( view_alloc() , N );
    V vb( view_alloc( "vb" ) , N );
    V vc( view_alloc( "vc" , AllowPadding ) , N );
    V vd( view_alloc( "vd" , WithoutInitializing ) , N );
    V ve( view_alloc( "ve" , WithoutInitializing , AllowPadding ) , N );
    V vf( view_alloc( "vf" , mem_space , WithoutInitializing , AllowPadding ) , N );
    V vg( view_alloc( mem_space , "vg" , WithoutInitializing , AllowPadding ) , N );
    V vh( view_alloc( WithoutInitializing , AllowPadding ) , N );
    V vi( view_alloc( WithoutInitializing ) , N );
    V vj( view_alloc( std::string("vj") , AllowPadding ) , N );
    V vk( view_alloc( mem_space , std::string("vk") , AllowPadding ) , N );
  }

  {
    typedef Kokkos::Experimental::ViewTraits<int***,Kokkos::LayoutStride,ExecSpace>  traits_t ;
    typedef Kokkos::Experimental::Impl::ViewDimension<0,0,0>                         dims_t ;
    typedef Kokkos::Experimental::Impl::ViewOffset< dims_t , Kokkos::LayoutStride >  offset_t ;

    Kokkos::LayoutStride stride ;

    stride.dimension[0] = 3 ;
    stride.dimension[1] = 4 ;
    stride.dimension[2] = 5 ;
    stride.stride[0] = 4 ;
    stride.stride[1] = 1 ;
    stride.stride[2] = 12 ;

    const offset_t offset( stride );

    ASSERT_EQ( offset.dimension_0() , 3 );
    ASSERT_EQ( offset.dimension_1() , 4 );
    ASSERT_EQ( offset.dimension_2() , 5 );

    ASSERT_EQ( offset.stride_0() , 4 );
    ASSERT_EQ( offset.stride_1() , 1 );
    ASSERT_EQ( offset.stride_2() , 12 );

    ASSERT_EQ( offset.span() , 60 );
    ASSERT_TRUE( offset.span_is_contiguous() );

    Kokkos::Experimental::Impl::ViewMapping< traits_t , void >  v( (int*) 0 , std::false_type() , stride );
  }

  {
    typedef Kokkos::Experimental::View<int**,ExecSpace>  V ;
    typedef typename V::HostMirror  M ;

    constexpr int N0 = 10 ;
    constexpr int N1 = 11 ;

    V a("a",N0,N1);
    M b = Kokkos::Experimental::create_mirror(a);
    M c = Kokkos::Experimental::create_mirror_view(a);
    M d ;

    for ( int i0 = 0 ; i0 < N0 ; ++i0 )
    for ( int i1 = 0 ; i1 < N1 ; ++i1 )
      b(i0,i1) = 1 + i0 + i1 * N0 ;

    Kokkos::Experimental::deep_copy( a , b );
    Kokkos::Experimental::deep_copy( c , a );

    for ( int i0 = 0 ; i0 < N0 ; ++i0 )
    for ( int i1 = 0 ; i1 < N1 ; ++i1 )
      ASSERT_EQ( b(i0,i1) , c(i0,i1) );

    Kokkos::Experimental::resize( b , 5 , 6 );
    Kokkos::Experimental::realloc( c , 5 , 6 );
    Kokkos::Experimental::realloc( d , 5 , 6 );

    ASSERT_EQ( b.dimension_0() , 5 );
    ASSERT_EQ( b.dimension_1() , 6 );
    ASSERT_EQ( c.dimension_0() , 5 );
    ASSERT_EQ( c.dimension_1() , 6 );
    ASSERT_EQ( d.dimension_0() , 5 );
    ASSERT_EQ( d.dimension_1() , 6 );
  }

#if defined( KOKKOS_USING_EXPERIMENTAL_VIEW )
  /* Only works when experimental view is activated */
  {
    typedef Kokkos::Experimental::View<int*,ExecSpace> V ;
    typedef Kokkos::Experimental::View<int*,ExecSpace,Kokkos::MemoryUnmanaged> U ;


    V a("a",10);

    ASSERT_EQ( a.use_count() , 1 );

    V b = a ;

    ASSERT_EQ( a.use_count() , 2 );
    ASSERT_EQ( b.use_count() , 2 );

    {
      U c = b ; // 'c' is compile-time unmanaged

      ASSERT_EQ( a.use_count() , 2 );
      ASSERT_EQ( b.use_count() , 2 );
      ASSERT_EQ( c.use_count() , 2 );

      V d = c ; // 'd' is run-time unmanaged

      ASSERT_EQ( a.use_count() , 2 );
      ASSERT_EQ( b.use_count() , 2 );
      ASSERT_EQ( c.use_count() , 2 );
      ASSERT_EQ( d.use_count() , 2 );
    }

    ASSERT_EQ( a.use_count() , 2 );
    ASSERT_EQ( b.use_count() , 2 );

    b = V();

    ASSERT_EQ( a.use_count() , 1 );
    ASSERT_EQ( b.use_count() , 0 );

    Kokkos::parallel_for(
      Kokkos::RangePolicy< Kokkos::DefaultHostExecutionSpace >(0,10) ,
      [=]( int i ){
        // 'a' is captured by copy and the capture mechanism
        // converts 'a' to an unmanaged copy.
        // When the parallel dispatch accepts a move for the lambda
        // this count should become 1
        ASSERT_EQ( a.use_count() , 2 );
        V x = a ;
        ASSERT_EQ( a.use_count() , 2 );
        ASSERT_EQ( x.use_count() , 2 );
      });
  }
#endif /* #if defined( KOKKOS_USING_EXPERIMENTAL_VIEW ) */
}

template< class ExecSpace >
struct TestViewMappingSubview {

  typedef Kokkos::pair<int,int> range ;

  enum { AN = 10 };
  typedef Kokkos::Experimental::View<int*,ExecSpace>  AT ;
  typedef Kokkos::Experimental::View<const int*,ExecSpace>  ACT ;
  typedef Kokkos::Experimental::Subview< AT , range >  AS ;

  enum { BN0 = 10 , BN1 = 11 , BN2 = 12 };
  typedef Kokkos::Experimental::View<int***,ExecSpace>  BT ;
  typedef Kokkos::Experimental::Subview< BT , range , range , range >  BS ;

  enum { CN0 = 10 , CN1 = 11 , CN2 = 12 };
  typedef Kokkos::Experimental::View<int***[13][14],ExecSpace>  CT ;
  typedef Kokkos::Experimental::Subview< CT , range , range , range , int , int >  CS ;

  enum { DN0 = 10 , DN1 = 11 , DN2 = 12 , DN3 = 13 , DN4 = 14 };
  typedef Kokkos::Experimental::View<int***[DN3][DN4],ExecSpace>  DT ;
  typedef Kokkos::Experimental::Subview< DT , int , range , range , range , int >  DS ;


  typedef Kokkos::Experimental::View<int***[13][14],Kokkos::LayoutLeft,ExecSpace>  DLT ;
  typedef Kokkos::Experimental::Subview< DLT , range , int , int , int , int >  DLS1 ;

  static_assert( DLS1::rank == 1 && std::is_same< typename DLS1::array_layout , Kokkos::LayoutLeft >::value
               , "Subview layout error for rank 1 subview of left-most range of LayoutLeft" );

  typedef Kokkos::Experimental::View<int***[13][14],Kokkos::LayoutRight,ExecSpace>  DRT ;
  typedef Kokkos::Experimental::Subview< DRT , int , int , int , int , range >  DRS1 ;

  static_assert( DRS1::rank == 1 && std::is_same< typename DRS1::array_layout , Kokkos::LayoutRight >::value
               , "Subview layout error for rank 1 subview of right-most range of LayoutRight" );

  AT Aa ;
  AS Ab ;
  ACT Ac ;
  BT Ba ;
  BS Bb ;
  CT Ca ;
  CS Cb ;
  DT Da ;
  DS Db ;

  TestViewMappingSubview()
    : Aa("Aa",AN)
    , Ab( Kokkos::Experimental::subview( Aa , std::pair<int,int>(1,AN-1) ) )
    , Ac( Aa , std::pair<int,int>(1,AN-1) )
    , Ba("Ba",BN0,BN1,BN2)
    , Bb( Kokkos::Experimental::subview( Ba
                                        , std::pair<int,int>(1,BN0-1)
                                        , std::pair<int,int>(1,BN1-1)
                                        , std::pair<int,int>(1,BN2-1)
                                        ) )
    , Ca("Ca",CN0,CN1,CN2)
    , Cb( Kokkos::Experimental::subview( Ca
                                        , std::pair<int,int>(1,CN0-1)
                                        , std::pair<int,int>(1,CN1-1)
                                        , std::pair<int,int>(1,CN2-1)
                                        , 1
                                        , 2
                                        ) )
    , Da("Da",DN0,DN1,DN2)
    , Db( Kokkos::Experimental::subview( Da
                                        , 1
                                        , std::pair<int,int>(1,DN1-1)
                                        , std::pair<int,int>(1,DN2-1)
                                        , std::pair<int,int>(1,DN3-1)
                                        , 2
                                        ) )
    {
    }


  KOKKOS_INLINE_FUNCTION
  void operator()( const int , long & error_count ) const
    {
      auto Ad = Kokkos::Experimental::subview< Kokkos::MemoryUnmanaged >( Aa , Kokkos::pair<int,int>(1,AN-1) );

      for ( int i = 1 ; i < AN-1 ; ++i ) if( & Aa[i] != & Ab[i-1] ) ++error_count ;
      for ( int i = 1 ; i < AN-1 ; ++i ) if( & Aa[i] != & Ac[i-1] ) ++error_count ;
      for ( int i = 1 ; i < AN-1 ; ++i ) if( & Aa[i] != & Ad[i-1] ) ++error_count ;

      for ( int i2 = 1 ; i2 < BN2-1 ; ++i2 ) {
      for ( int i1 = 1 ; i1 < BN1-1 ; ++i1 ) {
      for ( int i0 = 1 ; i0 < BN0-1 ; ++i0 ) {
        if ( & Ba(i0,i1,i2) != & Bb(i0-1,i1-1,i2-1) ) ++error_count ;
      }}}

      for ( int i2 = 1 ; i2 < CN2-1 ; ++i2 ) {
      for ( int i1 = 1 ; i1 < CN1-1 ; ++i1 ) {
      for ( int i0 = 1 ; i0 < CN0-1 ; ++i0 ) {
        if ( & Ca(i0,i1,i2,1,2) != & Cb(i0-1,i1-1,i2-1) ) ++error_count ;
      }}}

      for ( int i2 = 1 ; i2 < DN3-1 ; ++i2 ) {
      for ( int i1 = 1 ; i1 < DN2-1 ; ++i1 ) {
      for ( int i0 = 1 ; i0 < DN1-1 ; ++i0 ) {
        if ( & Da(1,i0,i1,i2,2) != & Db(i0-1,i1-1,i2-1) ) ++error_count ;
      }}}
    }

  static void run()
  {
    TestViewMappingSubview self ;

    ASSERT_EQ( self.Aa.dimension_0() , AN );
    ASSERT_EQ( self.Ab.dimension_0() , AN - 2 );
    ASSERT_EQ( self.Ac.dimension_0() , AN - 2 );
    ASSERT_EQ( self.Ba.dimension_0() , BN0 );
    ASSERT_EQ( self.Ba.dimension_1() , BN1 );
    ASSERT_EQ( self.Ba.dimension_2() , BN2 );
    ASSERT_EQ( self.Bb.dimension_0() , BN0 - 2 );
    ASSERT_EQ( self.Bb.dimension_1() , BN1 - 2 );
    ASSERT_EQ( self.Bb.dimension_2() , BN2 - 2 );

    ASSERT_EQ( self.Ca.dimension_0() , CN0 );
    ASSERT_EQ( self.Ca.dimension_1() , CN1 );
    ASSERT_EQ( self.Ca.dimension_2() , CN2 );
    ASSERT_EQ( self.Ca.dimension_3() , 13 );
    ASSERT_EQ( self.Ca.dimension_4() , 14 );
    ASSERT_EQ( self.Cb.dimension_0() , CN0 - 2 );
    ASSERT_EQ( self.Cb.dimension_1() , CN1 - 2 );
    ASSERT_EQ( self.Cb.dimension_2() , CN2 - 2 );

    ASSERT_EQ( self.Da.dimension_0() , DN0 );
    ASSERT_EQ( self.Da.dimension_1() , DN1 );
    ASSERT_EQ( self.Da.dimension_2() , DN2 );
    ASSERT_EQ( self.Da.dimension_3() , DN3 );
    ASSERT_EQ( self.Da.dimension_4() , DN4 );

    ASSERT_EQ( self.Db.dimension_0() , DN1 - 2 );
    ASSERT_EQ( self.Db.dimension_1() , DN2 - 2 );
    ASSERT_EQ( self.Db.dimension_2() , DN3 - 2 );

    ASSERT_EQ( self.Da.stride_1() , self.Db.stride_0() );
    ASSERT_EQ( self.Da.stride_2() , self.Db.stride_1() );
    ASSERT_EQ( self.Da.stride_3() , self.Db.stride_2() );

    long error_count = -1 ;
    Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >(0,1) , self , error_count );
    ASSERT_EQ( error_count , 0 );
  }

};

template< class ExecSpace >
void test_view_mapping_subview()
{
  TestViewMappingSubview< ExecSpace >::run();
}

/*--------------------------------------------------------------------------*/

template< class ViewType >
struct TestViewMapOperator {

  static_assert( ViewType::reference_type_is_lvalue_reference
               , "Test only valid for lvalue reference type" );

  const ViewType v ;

  KOKKOS_INLINE_FUNCTION
  void test_left( size_t i0 , long & error_count ) const
    {
      typename ViewType::value_type * const base_ptr = & v(0,0,0,0,0,0,0,0);
      const size_t n1 = v.dimension_1();
      const size_t n2 = v.dimension_2();
      const size_t n3 = v.dimension_3();
      const size_t n4 = v.dimension_4();
      const size_t n5 = v.dimension_5();
      const size_t n6 = v.dimension_6();
      const size_t n7 = v.dimension_7();

      long offset = 0 ;

      for ( size_t i7 = 0 ; i7 < n7 ; ++i7 )
      for ( size_t i6 = 0 ; i6 < n6 ; ++i6 )
      for ( size_t i5 = 0 ; i5 < n5 ; ++i5 )
      for ( size_t i4 = 0 ; i4 < n4 ; ++i4 )
      for ( size_t i3 = 0 ; i3 < n3 ; ++i3 )
      for ( size_t i2 = 0 ; i2 < n2 ; ++i2 )
      for ( size_t i1 = 0 ; i1 < n1 ; ++i1 )
      {
        const long d = & v(i0,i1,i2,i3,i4,i5,i6,i7) - base_ptr ;
        if ( d < offset ) ++error_count ;
        offset = d ;
      }

      if ( v.span() <= size_t(offset) ) ++error_count ;
    }

  KOKKOS_INLINE_FUNCTION
  void test_right( size_t i0 , long & error_count ) const
    {
      typename ViewType::value_type * const base_ptr = & v(0,0,0,0,0,0,0,0);
      const size_t n1 = v.dimension_1();
      const size_t n2 = v.dimension_2();
      const size_t n3 = v.dimension_3();
      const size_t n4 = v.dimension_4();
      const size_t n5 = v.dimension_5();
      const size_t n6 = v.dimension_6();
      const size_t n7 = v.dimension_7();

      long offset = 0 ;

      for ( size_t i1 = 0 ; i1 < n1 ; ++i1 )
      for ( size_t i2 = 0 ; i2 < n2 ; ++i2 )
      for ( size_t i3 = 0 ; i3 < n3 ; ++i3 )
      for ( size_t i4 = 0 ; i4 < n4 ; ++i4 )
      for ( size_t i5 = 0 ; i5 < n5 ; ++i5 )
      for ( size_t i6 = 0 ; i6 < n6 ; ++i6 )
      for ( size_t i7 = 0 ; i7 < n7 ; ++i7 )
      {
        const long d = & v(i0,i1,i2,i3,i4,i5,i6,i7) - base_ptr ;
        if ( d < offset ) ++error_count ;
        offset = d ;
      }

      if ( v.span() <= size_t(offset) ) ++error_count ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_t i , long & error_count ) const
    {
      if ( std::is_same< typename ViewType::array_layout , Kokkos::LayoutLeft >::value )
        test_left(i,error_count);
      else if ( std::is_same< typename ViewType::array_layout , Kokkos::LayoutRight >::value )
        test_right(i,error_count);
    }

  constexpr static size_t N0 = 10 ;
  constexpr static size_t N1 =  9 ;
  constexpr static size_t N2 =  8 ;
  constexpr static size_t N3 =  7 ;
  constexpr static size_t N4 =  6 ;
  constexpr static size_t N5 =  5 ;
  constexpr static size_t N6 =  4 ;
  constexpr static size_t N7 =  3 ;

  TestViewMapOperator() : v( "Test" , N0, N1, N2, N3, N4, N5, N6, N7 ) {}

  static void run()
    {
      TestViewMapOperator self ;

      ASSERT_EQ( self.v.dimension_0() , ( 0 < ViewType::rank ? N0 : 1 ) );
      ASSERT_EQ( self.v.dimension_1() , ( 1 < ViewType::rank ? N1 : 1 ) );
      ASSERT_EQ( self.v.dimension_2() , ( 2 < ViewType::rank ? N2 : 1 ) );
      ASSERT_EQ( self.v.dimension_3() , ( 3 < ViewType::rank ? N3 : 1 ) );
      ASSERT_EQ( self.v.dimension_4() , ( 4 < ViewType::rank ? N4 : 1 ) );
      ASSERT_EQ( self.v.dimension_5() , ( 5 < ViewType::rank ? N5 : 1 ) );
      ASSERT_EQ( self.v.dimension_6() , ( 6 < ViewType::rank ? N6 : 1 ) );
      ASSERT_EQ( self.v.dimension_7() , ( 7 < ViewType::rank ? N7 : 1 ) );

      ASSERT_LE( self.v.dimension_0()*
                 self.v.dimension_1()*
                 self.v.dimension_2()*
                 self.v.dimension_3()*
                 self.v.dimension_4()*
                 self.v.dimension_5()*
                 self.v.dimension_6()*
                 self.v.dimension_7()
               , self.v.span() );

      long error_count ;
      Kokkos::RangePolicy< typename ViewType::execution_space > range(0,self.v.dimension_0());
      Kokkos::parallel_reduce( range , self , error_count );
      ASSERT_EQ( 0 , error_count );
    }
};


template< class ExecSpace >
void test_view_mapping_operator()
{
  TestViewMapOperator< Kokkos::Experimental::View<int,Kokkos::LayoutLeft,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int*,Kokkos::LayoutLeft,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int**,Kokkos::LayoutLeft,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int***,Kokkos::LayoutLeft,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int****,Kokkos::LayoutLeft,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int*****,Kokkos::LayoutLeft,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int******,Kokkos::LayoutLeft,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int*******,Kokkos::LayoutLeft,ExecSpace> >::run();

  TestViewMapOperator< Kokkos::Experimental::View<int,Kokkos::LayoutRight,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int*,Kokkos::LayoutRight,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int**,Kokkos::LayoutRight,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int***,Kokkos::LayoutRight,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int****,Kokkos::LayoutRight,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int*****,Kokkos::LayoutRight,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int******,Kokkos::LayoutRight,ExecSpace> >::run();
  TestViewMapOperator< Kokkos::Experimental::View<int*******,Kokkos::LayoutRight,ExecSpace> >::run();
}

/*--------------------------------------------------------------------------*/

template< class ExecSpace >
struct TestViewMappingAtomic {
  typedef Kokkos::MemoryTraits< Kokkos::Atomic >  mem_trait ;

  typedef Kokkos::Experimental::View< int * , ExecSpace > T ;
  typedef Kokkos::Experimental::View< int * , ExecSpace , mem_trait >  T_atom ;

  T      x ;
  T_atom x_atom ;

  constexpr static size_t N = 100000 ;

  struct TagInit {};
  struct TagUpdate {};
  struct TagVerify {};

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagInit & , const int i ) const
    { x(i) = i ; }

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagUpdate & , const int i ) const
    { x_atom(i%2) += 1 ; }

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagVerify & , const int i , long & error_count ) const
    {
       if ( i < 2 ) { if ( x(i) != int(i + N / 2) ) ++error_count ; }
       else         { if ( x(i) != int(i) ) ++error_count ; }
    }

  TestViewMappingAtomic()
    : x("x",N)
    , x_atom( x )
    {}

  static void run()
    {
      ASSERT_TRUE( T::reference_type_is_lvalue_reference );
      ASSERT_FALSE( T_atom::reference_type_is_lvalue_reference );

      TestViewMappingAtomic self ;
      Kokkos::parallel_for( Kokkos::RangePolicy< ExecSpace , TagInit >(0,N) , self );
      Kokkos::parallel_for( Kokkos::RangePolicy< ExecSpace , TagUpdate >(0,N) , self );
      long error_count = -1 ;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , TagVerify >(0,N) , self , error_count );
      ASSERT_EQ( 0 , error_count );
    }
};

/*--------------------------------------------------------------------------*/

template< class ExecSpace >
struct TestViewMappingClassValue {

  struct ValueType {
    KOKKOS_INLINE_FUNCTION
    ValueType()
    {
#if 0
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA )
      printf("TestViewMappingClassValue construct on Cuda\n");
#elif defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      printf("TestViewMappingClassValue construct on Host\n");
#else
      printf("TestViewMappingClassValue construct unknown\n");
#endif
#endif
    }
    KOKKOS_INLINE_FUNCTION
    ~ValueType()
    {
#if 0
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA )
      printf("TestViewMappingClassValue destruct on Cuda\n");
#elif defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      printf("TestViewMappingClassValue destruct on Host\n");
#else
      printf("TestViewMappingClassValue destruct unknown\n");
#endif
#endif
    }
  };

  static void run()
  {
    using namespace Kokkos::Experimental ;
    ExecSpace::fence();
    {
      View< ValueType , ExecSpace > a("a");
      ExecSpace::fence();
    }
    ExecSpace::fence();
  }
};

} /* namespace Test */

/*--------------------------------------------------------------------------*/

