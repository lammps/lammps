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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef TEST_AGGREGATE_HPP
#define TEST_AGGREGATE_HPP

#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <impl/Kokkos_ViewArray.hpp>

namespace Test {

template< class DeviceType >
void TestViewAggregate()
{
  typedef Kokkos::Array< double, 32 >  value_type;
  typedef Kokkos::Impl::ViewDataAnalysis< value_type *, Kokkos::LayoutLeft, value_type > analysis_1d;

  static_assert( std::is_same< typename analysis_1d::specialize, Kokkos::Array<> >::value, "" );

  typedef Kokkos::ViewTraits< value_type **, DeviceType > a32_traits;
  typedef Kokkos::ViewTraits< typename a32_traits::scalar_array_type, DeviceType > flat_traits;

  static_assert( std::is_same< typename a32_traits::specialize, Kokkos::Array<> >::value, "" );
  static_assert( std::is_same< typename a32_traits::value_type, value_type >::value, "" );
  static_assert( a32_traits::rank == 2, "" );
  static_assert( a32_traits::rank_dynamic == 2, "" );

  static_assert( std::is_same< typename flat_traits::specialize, void >::value, "" );
  static_assert( flat_traits::rank == 3, "" );
  static_assert( flat_traits::rank_dynamic == 2, "" );
  static_assert( flat_traits::dimension::N2 == 32, "" );

  typedef Kokkos::View< Kokkos::Array< double, 32 > **, DeviceType > a32_type;
  typedef typename a32_type::array_type  a32_flat_type;

  static_assert( std::is_same< typename a32_type::value_type, value_type >::value, "" );
  static_assert( std::is_same< typename a32_type::pointer_type, double * >::value, "" );
  static_assert( a32_type::Rank == 2, "" );
  static_assert( a32_flat_type::Rank == 3, "" );

  a32_type x( "test", 4, 5 );
  a32_flat_type y( x );

  ASSERT_EQ( x.extent( 0 ), 4 );
  ASSERT_EQ( x.extent( 1 ), 5 );
  ASSERT_EQ( y.extent( 0 ), 4 );
  ASSERT_EQ( y.extent( 1 ), 5 );
  ASSERT_EQ( y.extent( 2 ), 32 );

  // Initialize arrays from brace-init-list as for std::array.
  //
  // Comment: Clang will issue the following warning if we don't use double
  //          braces here (one for initializing the Kokkos::Array and one for
  //          initializing the sub-aggreagate C-array data member),
  //
  //            warning: suggest braces around initialization of subobject
  //
  //          but single brace syntax would be valid as well.
  Kokkos::Array< float, 2 > aggregate_initialization_syntax_1 = { { 1.41, 3.14 } };
  ASSERT_FLOAT_EQ( aggregate_initialization_syntax_1[0], 1.41 );
  ASSERT_FLOAT_EQ( aggregate_initialization_syntax_1[1], 3.14 );

  Kokkos::Array< int, 3 > aggregate_initialization_syntax_2{ { 0, 1, 2 } }; // since C++11
  for ( int i = 0; i < 3; ++i ) {
    ASSERT_EQ( aggregate_initialization_syntax_2[i], i );
  }

  // Note that this is a valid initialization.
  Kokkos::Array< double, 3 > initialized_with_one_argument_missing = { { 255, 255 } };
  for (int i = 0; i < 2; ++i) {
    ASSERT_DOUBLE_EQ( initialized_with_one_argument_missing[i], 255 );
  }
  // But the following line would not compile
//  Kokkos::Array< double, 3 > initialized_with_too_many{ { 1, 2, 3, 4 } };
}

TEST_F( TEST_CATEGORY, view_aggregate )
{
  TestViewAggregate< TEST_EXECSPACE >();
}

} // namespace Test

#endif /* #ifndef TEST_AGGREGATE_HPP */
