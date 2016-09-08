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

/*--------------------------------------------------------------------------*/

namespace Test {
namespace {

template <typename ExecSpace >
struct TestMDRange_2D {

  using DataType     = int ;
  using ViewType     = typename Kokkos::View< DataType** ,  ExecSpace > ;
  using HostViewType = typename ViewType::HostMirror ;

  ViewType input_view ;

  TestMDRange_2D( const DataType N0, const DataType N1 ) : input_view("input_view", N0, N1) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i , const int j ) const
  {
    input_view(i,j) = 1;
  }


  static void test_for2( const int64_t N0, const int64_t N1 )
  {

    using namespace Kokkos::Experimental;

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<2>, Kokkos::IndexType<int> >;
      range_type range( {0,0}, {N0,N1} );
      TestMDRange_2D functor(N0,N1);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          if ( h_view(i,j) != 1 ) {
            ++counter;
          }
        }}
      if ( counter != 0 )
        printf(" Errors in test_for2; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<2, Iterate::Default, Iterate::Default >, Kokkos::IndexType<int> >;

      range_type range( {0,0}, {N0,N1} );
      TestMDRange_2D functor(N0,N1);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          if ( h_view(i,j) != 1 ) {
            ++counter;
          }
        }}
      if ( counter != 0 )
        printf(" Errors in test_for2; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<2, Iterate::Default, Iterate::Flat >, Kokkos::IndexType<int> >;

      range_type range( {0,0}, {N0,N1} );
      TestMDRange_2D functor(N0,N1);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          if ( h_view(i,j) != 1 ) {
            ++counter;
          }
        }}
      if ( counter != 0 )
        printf(" Errors in test_for2; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<2, Iterate::Right, Iterate::Flat >, Kokkos::IndexType<int> >;

      range_type range( {0,0}, {N0,N1} );
      TestMDRange_2D functor(N0,N1);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          if ( h_view(i,j) != 1 ) {
            ++counter;
          }
        }}
      if ( counter != 0 )
        printf(" Errors in test_for2; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<2, Iterate::Left, Iterate::Flat >, Kokkos::IndexType<int> >;

      range_type range( {0,0}, {N0,N1} );
      TestMDRange_2D functor(N0,N1);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          if ( h_view(i,j) != 1 ) {
            ++counter;
          }
        }}
      if ( counter != 0 )
        printf(" Errors in test_for2; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<2, Iterate::Left , Iterate::Left >, Kokkos::IndexType<int> >;

      range_type range( {0,0}, {N0,N1}, {3,3} );
      TestMDRange_2D functor(N0,N1);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          if ( h_view(i,j) != 1 ) {
            ++counter;
          }
        }}
      if ( counter != 0 )
        printf(" Errors in test_for2; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<2, Iterate::Left , Iterate::Right >, Kokkos::IndexType<int> >;

      range_type range( {0,0}, {N0,N1}, {7,7} );
      TestMDRange_2D functor(N0,N1);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          if ( h_view(i,j) != 1 ) {
            ++counter;
          }
        }}
      if ( counter != 0 )
        printf(" Errors in test_for2; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<2, Iterate::Right, Iterate::Left >, Kokkos::IndexType<int> >;

      range_type range( {0,0}, {N0,N1}, {16,16} );
      TestMDRange_2D functor(N0,N1);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          if ( h_view(i,j) != 1 ) {
            ++counter;
          }
        }}
      if ( counter != 0 )
        printf(" Errors in test_for2; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<2, Iterate::Right, Iterate::Right >, Kokkos::IndexType<int> >;

      range_type range( {0,0}, {N0,N1}, {5,16} );
      TestMDRange_2D functor(N0,N1);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          if ( h_view(i,j) != 1 ) {
            ++counter;
          }
        }}
      if ( counter != 0 )
        printf(" Errors in test_for2; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

  } //end test_for2
}; //MDRange_2D

template <typename ExecSpace >
struct TestMDRange_3D {

  using DataType = int ;
  using ViewType     = typename Kokkos::View< DataType*** ,  ExecSpace > ;
  using HostViewType = typename ViewType::HostMirror ;

  ViewType input_view ;

  TestMDRange_3D( const DataType N0, const DataType N1, const DataType N2 ) : input_view("input_view", N0, N1, N2) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i , const int j , const int k ) const
  {
    input_view(i,j,k) = 1;
  }

  static void test_for3( const int64_t N0, const int64_t N1, const int64_t N2 )
  {
    using namespace Kokkos::Experimental;

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3>, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3, Iterate::Default, Iterate::Default >, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3, Iterate::Flat, Iterate::Default>, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3, Iterate::Flat, Iterate::Flat >, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3, Iterate::Left, Iterate::Flat >, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3, Iterate::Right, Iterate::Flat >, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3, Iterate::Left, Iterate::Left >, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2}, {2,4,2} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3, Iterate::Left, Iterate::Right >, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2}, {3,5,7} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3, Iterate::Right, Iterate::Left >, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2}, {8,8,8} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

    {
      using range_type = MDRangePolicy< ExecSpace, Rank<3, Iterate::Right, Iterate::Right >, Kokkos::IndexType<int> >;

      range_type range( {0,0,0}, {N0,N1,N2}, {2,4,2} );
      TestMDRange_3D functor(N0,N1,N2);

      md_parallel_for( range, functor );

      HostViewType h_view = Kokkos::create_mirror_view( functor.input_view );
      Kokkos::deep_copy( h_view , functor.input_view );

      int counter = 0;
      for ( int i=0; i<N0; ++i ) {
        for ( int j=0; j<N1; ++j ) {
          for ( int k=0; k<N2; ++k ) {
          if ( h_view(i,j,k) != 1 ) {
            ++counter;
          }
        }}}
      if ( counter != 0 )
        printf(" Errors in test_for3; mismatches = %d\n\n",counter);
      ASSERT_EQ( counter , 0 );
    }

  } //end test_for3
};

} /* namespace */
} /* namespace Test */

/*--------------------------------------------------------------------------*/

