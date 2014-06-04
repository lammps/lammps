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


#include <sstream>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Shape.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void assert_counts_are_equal_throw(
  const unsigned x_count ,
  const unsigned y_count )
{
  std::ostringstream msg ;

  msg << "Kokkos::Impl::assert_counts_are_equal_throw( "
      << x_count << " != " << y_count << " )" ;

  throw_runtime_exception( msg.str() );
}

void assert_shapes_are_equal_throw(
  const unsigned x_scalar_size ,
  const unsigned x_rank ,
  const unsigned x_N0 , const unsigned x_N1 ,
  const unsigned x_N2 , const unsigned x_N3 ,
  const unsigned x_N4 , const unsigned x_N5 ,
  const unsigned x_N6 , const unsigned x_N7 ,

  const unsigned y_scalar_size ,
  const unsigned y_rank ,
  const unsigned y_N0 , const unsigned y_N1 ,
  const unsigned y_N2 , const unsigned y_N3 ,
  const unsigned y_N4 , const unsigned y_N5 ,
  const unsigned y_N6 , const unsigned y_N7 )
{
  std::ostringstream msg ;

  msg << "Kokkos::Impl::assert_shape_are_equal_throw( {"
      << " scalar_size(" << x_scalar_size
      << ") rank(" << x_rank
      << ") dimension(" ;
  if ( 0 < x_rank ) { msg << " " << x_N0 ; }
  if ( 1 < x_rank ) { msg << " " << x_N1 ; }
  if ( 2 < x_rank ) { msg << " " << x_N2 ; }
  if ( 3 < x_rank ) { msg << " " << x_N3 ; }
  if ( 4 < x_rank ) { msg << " " << x_N4 ; }
  if ( 5 < x_rank ) { msg << " " << x_N5 ; }
  if ( 6 < x_rank ) { msg << " " << x_N6 ; }
  if ( 7 < x_rank ) { msg << " " << x_N7 ; }
  msg << " ) } != { "
      << " scalar_size(" << y_scalar_size
      << ") rank(" << y_rank
      << ") dimension(" ;
  if ( 0 < y_rank ) { msg << " " << y_N0 ; }
  if ( 1 < y_rank ) { msg << " " << y_N1 ; }
  if ( 2 < y_rank ) { msg << " " << y_N2 ; }
  if ( 3 < y_rank ) { msg << " " << y_N3 ; }
  if ( 4 < y_rank ) { msg << " " << y_N4 ; }
  if ( 5 < y_rank ) { msg << " " << y_N5 ; }
  if ( 6 < y_rank ) { msg << " " << y_N6 ; }
  if ( 7 < y_rank ) { msg << " " << y_N7 ; }
  msg << " ) } )" ;

  throw_runtime_exception( msg.str() );
}

void AssertShapeBoundsAbort< Kokkos::HostSpace >::apply(
  const size_t rank ,
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
  std::ostringstream msg ;
  msg << "Kokkos::Impl::AssertShapeBoundsAbort( shape = {" ;
  if ( 0 < rank ) { msg << " " << n0 ; }
  if ( 1 < rank ) { msg << " " << n1 ; }
  if ( 2 < rank ) { msg << " " << n2 ; }
  if ( 3 < rank ) { msg << " " << n3 ; }
  if ( 4 < rank ) { msg << " " << n4 ; }
  if ( 5 < rank ) { msg << " " << n5 ; }
  if ( 6 < rank ) { msg << " " << n6 ; }
  if ( 7 < rank ) { msg << " " << n7 ; }
  msg << " } index = {" ;
  if ( 0 < arg_rank ) { msg << " " << i0 ; }
  if ( 1 < arg_rank ) { msg << " " << i1 ; }
  if ( 2 < arg_rank ) { msg << " " << i2 ; }
  if ( 3 < arg_rank ) { msg << " " << i3 ; }
  if ( 4 < arg_rank ) { msg << " " << i4 ; }
  if ( 5 < arg_rank ) { msg << " " << i5 ; }
  if ( 6 < arg_rank ) { msg << " " << i6 ; }
  if ( 7 < arg_rank ) { msg << " " << i7 ; }
  msg << " } )" ;

  throw_runtime_exception( msg.str() );
}

void assert_shape_effective_rank1_at_leastN_throw(
  const size_t x_rank , const size_t x_N0 ,
  const size_t x_N1 ,   const size_t x_N2 ,
  const size_t x_N3 ,   const size_t x_N4 ,
  const size_t x_N5 ,   const size_t x_N6 ,
  const size_t x_N7 ,
  const size_t N0 )
{
  std::ostringstream msg ;

  msg << "Kokkos::Impl::assert_shape_effective_rank1_at_leastN_throw( shape = {" ;
  if ( 0 < x_rank ) { msg << " " << x_N0 ; }
  if ( 1 < x_rank ) { msg << " " << x_N1 ; }
  if ( 2 < x_rank ) { msg << " " << x_N2 ; }
  if ( 3 < x_rank ) { msg << " " << x_N3 ; }
  if ( 4 < x_rank ) { msg << " " << x_N4 ; }
  if ( 5 < x_rank ) { msg << " " << x_N5 ; }
  if ( 6 < x_rank ) { msg << " " << x_N6 ; }
  if ( 7 < x_rank ) { msg << " " << x_N7 ; }
  msg << " } N = " << N0 << " )" ;

  throw_runtime_exception( msg.str() );
}



}
}

