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

#include <qthreads/TestQthreads.hpp>

namespace Test {

TEST_F( qthreads, impl_shared_alloc )
{
#if 0
  test_shared_alloc< Kokkos::HostSpace, Kokkos::Qthreads >();
#endif
}

TEST_F( qthreads, impl_view_mapping_b )
{
#if 0
  test_view_mapping_subview< Kokkos::Qthreads >();
  TestViewMappingAtomic< Kokkos::Qthreads >::run();
#endif
}

TEST_F( qthreads, view_api )
{
#if 0
  TestViewAPI< double, Kokkos::Qthreads >();
#endif
}

TEST_F( qthreads, view_nested_view )
{
#if 0
  ::Test::view_nested_view< Kokkos::Qthreads >();
#endif
}

TEST_F( qthreads, view_remap )
{
#if 0
  enum { N0 = 3, N1 = 2, N2 = 8, N3 = 9 };

  typedef Kokkos::View< double*[N1][N2][N3],
                        Kokkos::LayoutRight,
                        Kokkos::Qthreads > output_type;

  typedef Kokkos::View< int**[N2][N3],
                        Kokkos::LayoutLeft,
                        Kokkos::Qthreads > input_type;

  typedef Kokkos::View< int*[N0][N2][N3],
                        Kokkos::LayoutLeft,
                        Kokkos::Qthreads > diff_type;

  output_type output( "output", N0 );
  input_type  input ( "input", N0, N1 );
  diff_type   diff  ( "diff", N0 );

  int value = 0;

  for ( size_t i3 = 0; i3 < N3; ++i3 )
  for ( size_t i2 = 0; i2 < N2; ++i2 )
  for ( size_t i1 = 0; i1 < N1; ++i1 )
  for ( size_t i0 = 0; i0 < N0; ++i0 )
  {
    input( i0, i1, i2, i3 ) = ++value;
  }

  // Kokkos::deep_copy( diff, input ); // Throw with incompatible shape.
  Kokkos::deep_copy( output, input );

  value = 0;

  for ( size_t i3 = 0; i3 < N3; ++i3 )
  for ( size_t i2 = 0; i2 < N2; ++i2 )
  for ( size_t i1 = 0; i1 < N1; ++i1 )
  for ( size_t i0 = 0; i0 < N0; ++i0 )
  {
    ++value;
    ASSERT_EQ( value, ( (int) output( i0, i1, i2, i3 ) ) );
  }
#endif
}

TEST_F( qthreads, view_aggregate )
{
#if 0
  TestViewAggregate< Kokkos::Qthreads >();
#endif
}

TEST_F( qthreads, template_meta_functions )
{
#if 0
  TestTemplateMetaFunctions< int, Kokkos::Qthreads >();
#endif
}

} // namespace Test
