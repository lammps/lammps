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

inline
void test_utilities()
{
  using namespace Kokkos::Impl;

  {
    using i = integer_sequence< int >;
    using j = make_integer_sequence< int, 0 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 0u, "Error: integer_sequence.size()" );
  }

  {
    using i = integer_sequence< int, 0 >;
    using j = make_integer_sequence< int, 1 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 1u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = integer_sequence< int, 0, 1 >;
    using j = make_integer_sequence< int, 2 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 2u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 1, i >::value == 1, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 1, i{} ) == 1, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = integer_sequence< int, 0, 1, 2 >;
    using j = make_integer_sequence< int, 3 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 3u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 1, i >::value == 1, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 2, i >::value == 2, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 1, i{} ) == 1, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 2, i{} ) == 2, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = integer_sequence< int, 0, 1, 2, 3 >;
    using j = make_integer_sequence< int, 4 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 4u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 1, i >::value == 1, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 2, i >::value == 2, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 3, i >::value == 3, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 1, i{} ) == 1, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 2, i{} ) == 2, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 3, i{} ) == 3, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = integer_sequence< int, 0, 1, 2, 3, 4 >;
    using j = make_integer_sequence< int, 5 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 5u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 1, i >::value == 1, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 2, i >::value == 2, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 3, i >::value == 3, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 4, i >::value == 4, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 1, i{} ) == 1, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 2, i{} ) == 2, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 3, i{} ) == 3, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 4, i{} ) == 4, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = integer_sequence< int, 0, 1, 2, 3, 4, 5 >;
    using j = make_integer_sequence< int, 6 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 6u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 1, i >::value == 1, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 2, i >::value == 2, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 3, i >::value == 3, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 4, i >::value == 4, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 5, i >::value == 5, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 1, i{} ) == 1, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 2, i{} ) == 2, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 3, i{} ) == 3, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 4, i{} ) == 4, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 5, i{} ) == 5, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = integer_sequence< int, 0, 1, 2, 3, 4, 5, 6 >;
    using j = make_integer_sequence< int, 7 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 7u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 1, i >::value == 1, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 2, i >::value == 2, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 3, i >::value == 3, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 4, i >::value == 4, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 5, i >::value == 5, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 6, i >::value == 6, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 1, i{} ) == 1, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 2, i{} ) == 2, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 3, i{} ) == 3, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 4, i{} ) == 4, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 5, i{} ) == 5, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 6, i{} ) == 6, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = integer_sequence< int, 0, 1, 2, 3, 4, 5, 6, 7 >;
    using j = make_integer_sequence< int, 8 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 8u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 1, i >::value == 1, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 2, i >::value == 2, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 3, i >::value == 3, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 4, i >::value == 4, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 5, i >::value == 5, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 6, i >::value == 6, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 7, i >::value == 7, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 1, i{} ) == 1, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 2, i{} ) == 2, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 3, i{} ) == 3, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 4, i{} ) == 4, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 5, i{} ) == 5, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 6, i{} ) == 6, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 7, i{} ) == 7, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = integer_sequence< int, 0, 1, 2, 3, 4, 5, 6, 7, 8 >;
    using j = make_integer_sequence< int, 9 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 9u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 1, i >::value == 1, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 2, i >::value == 2, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 3, i >::value == 3, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 4, i >::value == 4, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 5, i >::value == 5, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 6, i >::value == 6, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 7, i >::value == 7, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 8, i >::value == 8, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 1, i{} ) == 1, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 2, i{} ) == 2, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 3, i{} ) == 3, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 4, i{} ) == 4, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 5, i{} ) == 5, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 6, i{} ) == 6, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 7, i{} ) == 7, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 8, i{} ) == 8, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = integer_sequence< int, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 >;
    using j = make_integer_sequence< int, 10 >;

    static_assert( std::is_same< i, j >::value, "Error: make_integer_sequence" );
    static_assert( i::size() == 10u, "Error: integer_sequence.size()" );

    static_assert( integer_sequence_at< 0, i >::value == 0, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 1, i >::value == 1, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 2, i >::value == 2, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 3, i >::value == 3, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 4, i >::value == 4, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 5, i >::value == 5, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 6, i >::value == 6, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 7, i >::value == 7, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 8, i >::value == 8, "Error: integer_sequence_at" );
    static_assert( integer_sequence_at< 9, i >::value == 9, "Error: integer_sequence_at" );

    static_assert( at( 0, i{} ) == 0, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 1, i{} ) == 1, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 2, i{} ) == 2, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 3, i{} ) == 3, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 4, i{} ) == 4, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 5, i{} ) == 5, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 6, i{} ) == 6, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 7, i{} ) == 7, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 8, i{} ) == 8, "Error: at(unsigned, integer_sequence)" );
    static_assert( at( 9, i{} ) == 9, "Error: at(unsigned, integer_sequence)" );
  }

  {
    using i = make_integer_sequence< int, 5 >;
    using r = reverse_integer_sequence< i >;
    using gr = integer_sequence< int, 4, 3, 2, 1, 0 >;

    static_assert( std::is_same< r, gr >::value, "Error: reverse_integer_sequence" );
  }

  {
    using s = make_integer_sequence< int, 10 >;
    using e = exclusive_scan_integer_sequence< s >;
    using i = inclusive_scan_integer_sequence< s >;

    using ge = integer_sequence< int, 0, 0, 1, 3, 6, 10, 15, 21, 28, 36 >;
    using gi = integer_sequence< int, 0, 1, 3, 6, 10, 15, 21, 28, 36, 45 >;

    static_assert( e::value == 45, "Error: scan value" );
    static_assert( i::value == 45, "Error: scan value" );

    static_assert( std::is_same< e::type, ge >::value, "Error: exclusive_scan" );
    static_assert( std::is_same< i::type, gi >::value, "Error: inclusive_scan" );
  }
}

} // namespace Test
