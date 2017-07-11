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

namespace Test {

template< class Space >
struct TestViewMappingSubview
{
  typedef typename Space::execution_space ExecSpace;
  typedef typename Space::memory_space    MemSpace;

  typedef Kokkos::pair< int, int > range;

  enum { AN = 10 };
  typedef Kokkos::View< int*, ExecSpace >  AT;
  typedef Kokkos::View< const int*, ExecSpace >  ACT;
  typedef Kokkos::Subview< AT, range >  AS;

  enum { BN0 = 10, BN1 = 11, BN2 = 12 };
  typedef Kokkos::View< int***, ExecSpace >  BT;
  typedef Kokkos::Subview< BT, range, range, range >  BS;

  enum { CN0 = 10, CN1 = 11, CN2 = 12 };
  typedef Kokkos::View< int***[13][14], ExecSpace >  CT;
  typedef Kokkos::Subview< CT, range, range, range, int, int >  CS;

  enum { DN0 = 10, DN1 = 11, DN2 = 12, DN3 = 13, DN4 = 14 };
  typedef Kokkos::View< int***[DN3][DN4], ExecSpace >  DT;
  typedef Kokkos::Subview< DT, int, range, range, range, int >  DS;

  typedef Kokkos::View< int***[13][14], Kokkos::LayoutLeft, ExecSpace >  DLT;
  typedef Kokkos::Subview< DLT, range, int, int, int, int >  DLS1;

  static_assert( DLS1::rank == 1 && std::is_same< typename DLS1::array_layout, Kokkos::LayoutLeft >::value
               , "Subview layout error for rank 1 subview of left-most range of LayoutLeft" );

  typedef Kokkos::View< int***[13][14], Kokkos::LayoutRight, ExecSpace >  DRT;
  typedef Kokkos::Subview< DRT, int, int, int, int, range >  DRS1;

  static_assert( DRS1::rank == 1 && std::is_same< typename DRS1::array_layout, Kokkos::LayoutRight >::value
               , "Subview layout error for rank 1 subview of right-most range of LayoutRight" );

  AT Aa;
  AS Ab;
  ACT Ac;
  BT Ba;
  BS Bb;
  CT Ca;
  CS Cb;
  DT Da;
  DS Db;

  TestViewMappingSubview()
    : Aa( "Aa", AN )
    , Ab( Kokkos::Experimental::subview( Aa, std::pair< int, int >( 1, AN - 1 ) ) )
    , Ac( Aa, std::pair< int, int >( 1, AN - 1 ) )
    , Ba( "Ba", BN0, BN1, BN2 )
    , Bb( Kokkos::Experimental::subview( Ba
                                        , std::pair< int, int >( 1, BN0 - 1 )
                                        , std::pair< int, int >( 1, BN1 - 1 )
                                        , std::pair< int, int >( 1, BN2 - 1 )
                                        ) )
    , Ca( "Ca", CN0, CN1, CN2 )
    , Cb( Kokkos::Experimental::subview( Ca
                                        , std::pair< int, int >( 1, CN0 - 1 )
                                        , std::pair< int, int >( 1, CN1 - 1 )
                                        , std::pair< int, int >( 1, CN2 - 1 )
                                        , 1
                                        , 2
                                        ) )
    , Da( "Da", DN0, DN1, DN2 )
    , Db( Kokkos::Experimental::subview( Da
                                        , 1
                                        , std::pair< int, int >( 1, DN1 - 1 )
                                        , std::pair< int, int >( 1, DN2 - 1 )
                                        , std::pair< int, int >( 1, DN3 - 1 )
                                        , 2
                                        ) )
    {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int, long & error_count ) const
  {
    auto Ad = Kokkos::Experimental::subview< Kokkos::MemoryUnmanaged >( Aa, Kokkos::pair< int, int >( 1, AN - 1 ) );

    for ( int i = 1; i < AN - 1; ++i ) if( & Aa[i] != & Ab[i - 1] ) ++error_count;
    for ( int i = 1; i < AN - 1; ++i ) if( & Aa[i] != & Ac[i - 1] ) ++error_count;
    for ( int i = 1; i < AN - 1; ++i ) if( & Aa[i] != & Ad[i - 1] ) ++error_count;

    for ( int i2 = 1; i2 < BN2 - 1; ++i2 )
    for ( int i1 = 1; i1 < BN1 - 1; ++i1 )
    for ( int i0 = 1; i0 < BN0 - 1; ++i0 )
    {
      if ( & Ba( i0, i1, i2 ) != & Bb( i0 - 1, i1 - 1, i2 - 1 ) ) ++error_count;
    }

    for ( int i2 = 1; i2 < CN2 - 1; ++i2 )
    for ( int i1 = 1; i1 < CN1 - 1; ++i1 )
    for ( int i0 = 1; i0 < CN0 - 1; ++i0 )
    {
      if ( & Ca( i0, i1, i2, 1, 2 ) != & Cb( i0 - 1, i1 - 1, i2 - 1 ) ) ++error_count;
    }

    for ( int i2 = 1; i2 < DN3 - 1; ++i2 )
    for ( int i1 = 1; i1 < DN2 - 1; ++i1 )
    for ( int i0 = 1; i0 < DN1 - 1; ++i0 )
    {
      if ( & Da( 1, i0, i1, i2, 2 ) != & Db( i0 - 1, i1 - 1, i2 - 1 ) ) ++error_count;
    }
  }

  void run()
  {
    typedef typename Space::execution_space ExecSpace;

    TestViewMappingSubview< ExecSpace > self;

    ASSERT_EQ( Aa.dimension_0(), AN );
    ASSERT_EQ( Ab.dimension_0(), AN - 2 );
    ASSERT_EQ( Ac.dimension_0(), AN - 2 );
    ASSERT_EQ( Ba.dimension_0(), BN0 );
    ASSERT_EQ( Ba.dimension_1(), BN1 );
    ASSERT_EQ( Ba.dimension_2(), BN2 );
    ASSERT_EQ( Bb.dimension_0(), BN0 - 2 );
    ASSERT_EQ( Bb.dimension_1(), BN1 - 2 );
    ASSERT_EQ( Bb.dimension_2(), BN2 - 2 );

    ASSERT_EQ( Ca.dimension_0(), CN0 );
    ASSERT_EQ( Ca.dimension_1(), CN1 );
    ASSERT_EQ( Ca.dimension_2(), CN2 );
    ASSERT_EQ( Ca.dimension_3(), 13 ); 
    ASSERT_EQ( Ca.dimension_4(), 14 );
    ASSERT_EQ( Cb.dimension_0(), CN0 - 2 );
    ASSERT_EQ( Cb.dimension_1(), CN1 - 2 );
    ASSERT_EQ( Cb.dimension_2(), CN2 - 2 );

    ASSERT_EQ( Da.dimension_0(), DN0 );
    ASSERT_EQ( Da.dimension_1(), DN1 );
    ASSERT_EQ( Da.dimension_2(), DN2 );
    ASSERT_EQ( Da.dimension_3(), DN3 );
    ASSERT_EQ( Da.dimension_4(), DN4 );

    ASSERT_EQ( Db.dimension_0(), DN1 - 2 );
    ASSERT_EQ( Db.dimension_1(), DN2 - 2 );
    ASSERT_EQ( Db.dimension_2(), DN3 - 2 );

    ASSERT_EQ( Da.stride_1(), Db.stride_0() );
    ASSERT_EQ( Da.stride_2(), Db.stride_1() );
    ASSERT_EQ( Da.stride_3(), Db.stride_2() );

    long error_count = -1;
    Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace >( 0, 1 ), *this, error_count );
    ASSERT_EQ( error_count, 0 );
  }
};

TEST_F( TEST_CATEGORY , view_mapping_subview )
{
  TestViewMappingSubview< TEST_EXECSPACE > f;
  f.run();
}

}
