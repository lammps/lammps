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

#include <stdio.h>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace Test {

namespace {

template <typename ExecSpace >
struct TestMDRange_2D {
  using DataType     = int;
  using ViewType     = typename Kokkos::View< DataType**, ExecSpace >;
  using HostViewType = typename ViewType::HostMirror;

  ViewType input_view;

  TestMDRange_2D( const DataType N0, const DataType N1 ) : input_view( "input_view", N0, N1 ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j ) const
  {
    input_view( i, j ) = 1;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j, double &lsum ) const
  {
    lsum += input_view( i, j ) * 2;
  }

  // tagged operators
  struct InitTag {};
  KOKKOS_INLINE_FUNCTION
  void operator()( const InitTag &, const int i, const int j ) const
  {
    input_view( i, j ) = 3;
  }

  static void test_reduce2( const int N0, const int N1 )
  {
    using namespace Kokkos::Experimental;

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 3, 3 } } );

      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Default, Iterate::Default>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 2, 6 } } );

      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Left, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 2, 6 } } );

      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Left, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 2, 6 } } );

      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Right, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 2, 6 } } );

      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Right, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 2, 6 } } );

      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 );
    }
  } // end test_reduce2

  static void test_for2( const int N0, const int N1 )
  {
    using namespace Kokkos::Experimental;

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2>, Kokkos::IndexType<int>, InitTag > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 3, 3 } } );
      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      {
        if ( h_view( i, j ) != 3 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Default Layouts + InitTag op(): Errors in test_for2; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2>, InitTag > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 3, 3 } } );
      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      {
        if ( h_view( i, j ) != 3 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Default Layouts + InitTag op(): Errors in test_for2; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2>, InitTag > range_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } } );
      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      {
        if ( h_view( i, j ) != 3 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Default Layouts + InitTag op() + Default Tile: Errors in test_for2; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 3, 3 } } );
      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      {
        if ( h_view( i, j ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "No info: Errors in test_for2; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Default, Iterate::Default>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 4, 4 } } );
      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      {
        if ( h_view( i, j ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "D D: Errors in test_for2; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Left, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 3, 3 } } );
      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      {
        if ( h_view( i, j ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "L L: Errors in test_for2; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Left, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 7, 7 } } );
      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      {
        if ( h_view( i, j ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "L R: Errors in test_for2; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Right, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 16, 16 } } );
      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      {
        if ( h_view( i, j ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "R L: Errors in test_for2; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<2, Iterate::Right, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0 } }, point_type{ { N0, N1 } }, tile_type{ { 5, 16 } } );
      TestMDRange_2D functor( N0, N1 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      {
        if ( h_view( i, j ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "R R: Errors in test_for2; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }
  } // end test_for2
}; // MDRange_2D

template <typename ExecSpace >
struct TestMDRange_3D {
  using DataType     = int;
  using ViewType     = typename Kokkos::View< DataType***, ExecSpace >;
  using HostViewType = typename ViewType::HostMirror;

  ViewType input_view;

  TestMDRange_3D( const DataType N0, const DataType N1, const DataType N2 ) : input_view( "input_view", N0, N1, N2 ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j, const int k ) const
  {
    input_view( i, j, k ) = 1;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j, const int k, double &lsum ) const
  {
    lsum += input_view( i, j, k ) * 2;
  }

  // tagged operators
  struct InitTag {};
  KOKKOS_INLINE_FUNCTION
  void operator()( const InitTag &, const int i, const int j, const int k ) const
  {
    input_view( i, j, k ) = 3;
  }

  static void test_reduce3( const int N0, const int N1, const int N2 )
  {
    using namespace Kokkos::Experimental;

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 3, 3, 3 } } );

      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 * N2 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Default, Iterate::Default >, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 2, 4, 6 } } );

      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 * N2 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Left, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 2, 4, 6 } } );

      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 * N2 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Left, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 2, 4, 6 } } );

      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 * N2 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Right, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 2, 4, 6 } } );

      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 * N2 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Right, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 2, 4, 6 } } );

      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );
      double sum = 0.0;
      md_parallel_reduce( range, functor, sum );

      ASSERT_EQ( sum, 2 * N0 * N1 * N2 );
    }
  } // end test_reduce3

  static void test_for3( const int N0, const int N1, const int N2 )
  {
    using namespace Kokkos::Experimental;

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3> > range_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } } );
      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      {
        if ( h_view( i, j, k ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Defaults + No Tile: Errors in test_for3; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3>, Kokkos::IndexType<int>, InitTag > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 3, 3, 3 } } );
      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      {
        if ( h_view( i, j, k ) != 3 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Defaults + InitTag op(): Errors in test_for3; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 3, 3, 3 } } );

      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      {
        if ( h_view( i, j, k ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for3; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Default, Iterate::Default>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 3, 3, 3 } } );
      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      {
        if ( h_view( i, j, k ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for3; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Left, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 2, 4, 2 } } );
      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      {
        if ( h_view( i, j, k ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for3; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Left, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 3, 5, 7 } } );
      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      {
        if ( h_view( i, j, k ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for3; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Right, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 8, 8, 8 } } );
      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      {
        if ( h_view( i, j, k ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for3; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<3, Iterate::Right, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0 } }, point_type{ { N0, N1, N2 } }, tile_type{ { 2, 4, 2 } } );
      TestMDRange_3D functor( N0, N1, N2 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      {
        if ( h_view( i, j, k ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for3; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }
  } // end test_for3
};

template <typename ExecSpace >
struct TestMDRange_4D {
  using DataType     = int;
  using ViewType     = typename Kokkos::View< DataType****, ExecSpace >;
  using HostViewType = typename ViewType::HostMirror;

  ViewType input_view;

  TestMDRange_4D( const DataType N0, const DataType N1, const DataType N2, const DataType N3 ) : input_view( "input_view", N0, N1, N2, N3 ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j, const int k, const int l ) const
  {
    input_view( i, j, k, l ) = 1;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j, const int k, const int l, double &lsum ) const
  {
    lsum += input_view( i, j, k, l ) * 2;
  }

  // tagged operators
  struct InitTag {};
  KOKKOS_INLINE_FUNCTION
  void operator()( const InitTag &, const int i, const int j, const int k, const int l ) const
  {
    input_view( i, j, k, l ) = 3;
  }

  static void test_for4( const int N0, const int N1, const int N2, const int N3 )
  {
    using namespace Kokkos::Experimental;

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<4> > range_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3 } } );
      TestMDRange_4D functor( N0, N1, N2, N3 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      {
        if ( h_view( i, j, k, l ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Defaults + No Tile: Errors in test_for4; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<4>, Kokkos::IndexType<int>, InitTag > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3 } }, tile_type{ { 3, 11, 3, 3 } } );
      TestMDRange_4D functor( N0, N1, N2, N3 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      {
        if ( h_view( i, j, k, l ) != 3 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf("Defaults +m_tile > m_upper dim2 InitTag op(): Errors in test_for4; mismatches = %d\n\n",counter);
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<4>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3 } }, tile_type{ { 4, 4, 4, 4 } } );

      TestMDRange_4D functor( N0, N1, N2, N3 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      {
        if ( h_view( i, j, k, l ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for4; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<4, Iterate::Default, Iterate::Default>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3 } }, tile_type{ { 4, 4, 4, 4 } } );

      TestMDRange_4D functor( N0, N1, N2, N3 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      {
        if ( h_view( i, j, k, l ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for4; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<4, Iterate::Left, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3 } }, tile_type{ { 4, 4, 4, 4 } } );

      TestMDRange_4D functor( N0, N1, N2, N3 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      {
        if ( h_view( i, j, k, l ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for4; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<4, Iterate::Left, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3 } }, tile_type{ { 4, 4, 4, 4 } } );

      TestMDRange_4D functor( N0, N1, N2, N3 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      {
        if ( h_view( i, j, k, l ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for4; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<4, Iterate::Right, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3 } }, tile_type{ { 4, 4, 4, 4 } } );

      TestMDRange_4D functor( N0, N1, N2, N3 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      {
        if ( h_view( i, j, k, l ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for4; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<4, Iterate::Right, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3 } }, tile_type{ { 4, 4, 4, 4 } } );

      TestMDRange_4D functor( N0, N1, N2, N3 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      {
        if ( h_view( i, j, k, l ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for4; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }
  } // end test_for4
};

template <typename ExecSpace >
struct TestMDRange_5D {
  using DataType     = int;
  using ViewType     = typename Kokkos::View< DataType*****, ExecSpace >;
  using HostViewType = typename ViewType::HostMirror;

  ViewType input_view;

  TestMDRange_5D( const DataType N0, const DataType N1, const DataType N2, const DataType N3, const DataType N4 ) : input_view( "input_view", N0, N1, N2, N3, N4 ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j, const int k, const int l, const int m ) const
  {
    input_view( i, j, k, l, m ) = 1;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j, const int k, const int l, const int m, double &lsum ) const
  {
    lsum += input_view( i, j, k, l, m ) * 2;
  }

  // tagged operators
  struct InitTag {};
  KOKKOS_INLINE_FUNCTION
  void operator()( const InitTag &, const int i, const int j, const int k, const int l, const int m ) const
  {
    input_view( i, j, k, l, m ) = 3;
  }

  static void test_for5( const int N0, const int N1, const int N2, const int N3, const int N4 )
  {
    using namespace Kokkos::Experimental;

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<5> > range_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4 } } );
      TestMDRange_5D functor( N0, N1, N2, N3, N4 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      {
        if ( h_view( i, j, k, l, m ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Defaults + No Tile: Errors in test_for5; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<5>, Kokkos::IndexType<int>, InitTag > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4 } }, tile_type{ { 3, 3, 3, 3, 7 } } );
      TestMDRange_5D functor( N0, N1, N2, N3, N4 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      {
        if ( h_view( i, j, k, l, m ) != 3 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Defaults + InitTag op(): Errors in test_for5; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<5>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4 } }, tile_type{ { 4, 4, 4, 2, 2 } } );

      TestMDRange_5D functor( N0, N1, N2, N3, N4 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      {
        if ( h_view( i, j, k, l, m ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for5; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<5, Iterate::Default, Iterate::Default>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4 } }, tile_type{ { 4, 4, 4, 2, 2 } } );

      TestMDRange_5D functor( N0, N1, N2, N3, N4 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      {
        if ( h_view( i, j, k, l, m ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for5; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<5, Iterate::Left, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4 } }, tile_type{ { 4, 4, 4, 2, 2 } } );

      TestMDRange_5D functor( N0, N1, N2, N3, N4 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      {
        if ( h_view( i, j, k, l, m ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for5; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<5, Iterate::Left, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4 } }, tile_type{ { 4, 4, 4, 2, 2 } } );

      TestMDRange_5D functor( N0, N1, N2, N3, N4 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      {
        if ( h_view( i, j, k, l, m ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for5; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<5, Iterate::Right, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4 } }, tile_type{ { 4, 4, 4, 2, 2 } } );

      TestMDRange_5D functor( N0, N1, N2, N3, N4 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      {
        if ( h_view( i, j, k, l, m ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for5; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<5, Iterate::Right, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4 } }, tile_type{ { 4, 4, 4, 2, 2 } } );

      TestMDRange_5D functor( N0, N1, N2, N3, N4 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      {
        if ( h_view( i, j, k, l, m ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for5; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }
  }
};

template <typename ExecSpace >
struct TestMDRange_6D {
  using DataType     = int;
  using ViewType     = typename Kokkos::View< DataType******, ExecSpace >;
  using HostViewType = typename ViewType::HostMirror;

  ViewType input_view;

  TestMDRange_6D( const DataType N0, const DataType N1, const DataType N2, const DataType N3, const DataType N4, const DataType N5 ) : input_view( "input_view", N0, N1, N2, N3, N4, N5 ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j, const int k, const int l, const int m, const int n ) const
  {
    input_view( i, j, k, l, m, n ) = 1;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i, const int j, const int k, const int l, const int m, const int n, double &lsum ) const
  {
    lsum += input_view( i, j, k, l, m, n ) * 2;
  }

  // tagged operators
  struct InitTag {};
  KOKKOS_INLINE_FUNCTION
  void operator()( const InitTag &, const int i, const int j, const int k, const int l, const int m, const int n ) const
  {
    input_view( i, j, k, l, m, n ) = 3;
  }

  static void test_for6( const int N0, const int N1, const int N2, const int N3, const int N4, const int N5 )
  {
    using namespace Kokkos::Experimental;

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<6> > range_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4, N5 } } );
      TestMDRange_6D functor( N0, N1, N2, N3, N4, N5 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      for ( int n = 0; n < N5; ++n )
      {
        if ( h_view( i, j, k, l, m, n ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Defaults + No Tile: Errors in test_for6; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<6>, Kokkos::IndexType<int>, InitTag > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4, N5 } }, tile_type{ { 3, 3, 3, 3, 2, 3 } } ); //tile dims 3,3,3,3,3,3 more than cuda can handle with debugging
      TestMDRange_6D functor( N0, N1, N2, N3, N4, N5 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      for ( int n = 0; n < N5; ++n )
      {
        if ( h_view( i, j, k, l, m, n ) != 3 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( "Defaults + InitTag op(): Errors in test_for6; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<6>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4, N5 } }, tile_type{ { 4, 4, 4, 2, 2, 2 } } );

      TestMDRange_6D functor( N0, N1, N2, N3, N4, N5 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      for ( int n = 0; n < N5; ++n )
      {
        if ( h_view( i, j, k, l, m, n ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for6; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<6, Iterate::Default, Iterate::Default>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4, N5 } }, tile_type{ { 4, 4, 4, 2, 2, 2 } } );

      TestMDRange_6D functor( N0, N1, N2, N3, N4, N5 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      for ( int n = 0; n < N5; ++n )
      {
        if ( h_view( i, j, k, l, m, n ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for6; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<6, Iterate::Left, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4, N5 } }, tile_type{ { 4, 4, 4, 2, 2, 2 } } );

      TestMDRange_6D functor( N0, N1, N2, N3, N4, N5 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      for ( int n = 0; n < N5; ++n )
      {
        if ( h_view( i, j, k, l, m, n ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for6; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<6, Iterate::Left, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4, N5 } }, tile_type{ { 4, 4, 4, 2, 2, 2 } } );

      TestMDRange_6D functor( N0, N1, N2, N3, N4, N5 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      for ( int n = 0; n < N5; ++n )
      {
        if ( h_view( i, j, k, l, m, n ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for6; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<6, Iterate::Right, Iterate::Left>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4, N5 } }, tile_type{ { 4, 4, 4, 2, 2, 2 } } );

      TestMDRange_6D functor( N0, N1, N2, N3, N4, N5 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      for ( int n = 0; n < N5; ++n )
      {
        if ( h_view( i, j, k, l, m, n ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for6; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }

    {
      typedef typename Kokkos::Experimental::MDRangePolicy< ExecSpace, Rank<6, Iterate::Right, Iterate::Right>, Kokkos::IndexType<int> > range_type;
      typedef typename range_type::tile_type tile_type;
      typedef typename range_type::point_type point_type;

      range_type range( point_type{ { 0, 0, 0, 0, 0, 0 } }, point_type{ { N0, N1, N2, N3, N4, N5 } }, tile_type{ { 4, 4, 4, 2, 2, 2 } } );

      TestMDRange_6D functor( N0, N1, N2, N3, N4, N5 );

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view, functor.input_view );

      int counter = 0;
      for ( int i = 0; i < N0; ++i )
      for ( int j = 0; j < N1; ++j )
      for ( int k = 0; k < N2; ++k )
      for ( int l = 0; l < N3; ++l )
      for ( int m = 0; m < N4; ++m )
      for ( int n = 0; n < N5; ++n )
      {
        if ( h_view( i, j, k, l, m, n ) != 1 ) {
          ++counter;
        }
      }

      if ( counter != 0 ) {
        printf( " Errors in test_for6; mismatches = %d\n\n", counter );
      }

      ASSERT_EQ( counter, 0 );
    }
  }
};

} // namespace

} // namespace Test
