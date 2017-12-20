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

#include <vector>

#include <Kokkos_StaticCrsGraph.hpp>

/*--------------------------------------------------------------------------*/

namespace TestStaticCrsGraph {

template< class Space >
void run_test_graph()
{
  typedef Kokkos::StaticCrsGraph< unsigned , Space > dView ;
  typedef typename dView::HostMirror hView ;

  const unsigned LENGTH = 1000 ;
  dView dx ;
  hView hx ;

  std::vector< std::vector< int > > graph( LENGTH );

  for ( size_t i = 0 ; i < LENGTH ; ++i ) {
    graph[i].reserve(8);
    for ( size_t j = 0 ; j < 8 ; ++j ) {
      graph[i].push_back( i + j * 3 );
    }
  }

  dx = Kokkos::create_staticcrsgraph<dView>( "dx" , graph );
  hx = Kokkos::create_mirror( dx );

  ASSERT_EQ( hx.row_map.dimension_0() - 1 , LENGTH );

  for ( size_t i = 0 ; i < LENGTH ; ++i ) {
    const size_t begin = hx.row_map[i];
    const size_t n = hx.row_map[i+1] - begin ;
    ASSERT_EQ( n , graph[i].size() );
    for ( size_t j = 0 ; j < n ; ++j ) {
      ASSERT_EQ( (int) hx.entries( j + begin ) , graph[i][j] );
    }
  }

  // Test row view access
  for ( size_t i = 0 ; i < LENGTH ; ++i ) {
    auto rowView = hx.rowConst(i);
    ASSERT_EQ( rowView.length, graph[i].size() );
    for ( size_t j = 0 ; j < rowView.length ; ++j ) {
      ASSERT_EQ( rowView.colidx( j ) , graph[i][j] );
      ASSERT_EQ( rowView( j )        , graph[i][j] );
    }
  }
}

template< class Space >
void run_test_graph2()
{
  typedef Kokkos::StaticCrsGraph< unsigned[3] , Space > dView ;
  typedef typename dView::HostMirror hView ;

  const unsigned LENGTH = 10 ;

  std::vector< size_t > sizes( LENGTH );

  size_t total_length = 0 ;

  for ( size_t i = 0 ; i < LENGTH ; ++i ) {
    total_length += ( sizes[i] = 6 + i % 4 );
  }

  dView dx = Kokkos::create_staticcrsgraph<dView>( "test" , sizes );
  hView hx = Kokkos::create_mirror( dx );
  hView mx = Kokkos::create_mirror( dx );

  ASSERT_EQ( (size_t) dx.row_map.dimension_0() , (size_t) LENGTH + 1 );
  ASSERT_EQ( (size_t) hx.row_map.dimension_0() , (size_t) LENGTH + 1 );
  ASSERT_EQ( (size_t) mx.row_map.dimension_0() , (size_t) LENGTH + 1 );

  ASSERT_EQ( (size_t) dx.entries.dimension_0() , (size_t) total_length );
  ASSERT_EQ( (size_t) hx.entries.dimension_0() , (size_t) total_length );
  ASSERT_EQ( (size_t) mx.entries.dimension_0() , (size_t) total_length );

  ASSERT_EQ( (size_t) dx.entries.dimension_1() , (size_t) 3 );
  ASSERT_EQ( (size_t) hx.entries.dimension_1() , (size_t) 3 );
  ASSERT_EQ( (size_t) mx.entries.dimension_1() , (size_t) 3 );

  for ( size_t i = 0 ; i < LENGTH ; ++i ) {
    const size_t entry_begin = hx.row_map[i];
    const size_t entry_end   = hx.row_map[i+1];
    for ( size_t j = entry_begin ; j < entry_end ; ++j ) {
      hx.entries(j,0) = j + 1 ;
      hx.entries(j,1) = j + 2 ;
      hx.entries(j,2) = j + 3 ;
    }
  }

  Kokkos::deep_copy( dx.entries , hx.entries );
  Kokkos::deep_copy( mx.entries , dx.entries );

  ASSERT_EQ( mx.row_map.dimension_0() , (size_t) LENGTH + 1 );

  for ( size_t i = 0 ; i < LENGTH ; ++i ) {
    const size_t entry_begin = mx.row_map[i];
    const size_t entry_end   = mx.row_map[i+1];
    ASSERT_EQ( ( entry_end - entry_begin ) , sizes[i] );
    for ( size_t j = entry_begin ; j < entry_end ; ++j ) {
      ASSERT_EQ( (size_t) mx.entries( j , 0 ) , ( j + 1 ) );
      ASSERT_EQ( (size_t) mx.entries( j , 1 ) , ( j + 2 ) );
      ASSERT_EQ( (size_t) mx.entries( j , 2 ) , ( j + 3 ) );
    }
  }
}

template< class Space >
void run_test_graph3(size_t B, size_t N)
{
  srand(10310);

  typedef Kokkos::StaticCrsGraph< int , Space > dView ;
  typedef typename dView::HostMirror hView ;

  const unsigned LENGTH = 2000 ;

  std::vector< size_t > sizes( LENGTH );

  size_t total_length = 0 ;

  for ( size_t i = 0 ; i < LENGTH ; ++i ) {
    sizes[i] = rand()%1000;
  }

  sizes[1] = N;
  sizes[1998] = N;

  for ( size_t i = 0 ; i < LENGTH ; ++i ) {
    total_length += sizes[i];
  }

  int C = 0;
  dView dx = Kokkos::create_staticcrsgraph<dView>( "test" , sizes );
  dx.create_block_partitioning(B,C);
  hView hx = Kokkos::create_mirror( dx );

  for( size_t i = 0; i<B; i++) {
    size_t ne = 0;
    for(size_t j = hx.row_block_offsets(i); j<hx.row_block_offsets(i+1); j++)
      ne += hx.row_map(j+1)-hx.row_map(j)+C;

    ASSERT_FALSE((ne>2*((hx.row_map(hx.numRows())+C*hx.numRows())/B))&&(hx.row_block_offsets(i+1)>hx.row_block_offsets(i)+1));
  }
}

} /* namespace TestStaticCrsGraph */

