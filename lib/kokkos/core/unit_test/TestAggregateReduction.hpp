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

#ifndef TEST_AGGREGATE_REDUCTION_HPP
#define TEST_AGGREGATE_REDUCTION_HPP

#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <iostream>

namespace Test {

template< typename T , unsigned N >
struct StaticArray {
  T value[N] ;

  KOKKOS_INLINE_FUNCTION
  StaticArray()
    { for ( unsigned i = 0 ; i < N ; ++i ) value[i] = T(); }

  KOKKOS_INLINE_FUNCTION
  StaticArray( const StaticArray & rhs )
    { for ( unsigned i = 0 ; i < N ; ++i ) value[i] = rhs.value[i]; }

  KOKKOS_INLINE_FUNCTION
  operator T () { return value[0]; }

  KOKKOS_INLINE_FUNCTION
  StaticArray & operator = ( const T & rhs )
    {
      for ( unsigned i = 0 ; i < N ; ++i ) value[i] = rhs ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION
  StaticArray & operator = ( const StaticArray & rhs )
    {
      for ( unsigned i = 0 ; i < N ; ++i ) value[i] = rhs.value[i] ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION
  StaticArray operator * ( const StaticArray & rhs )
    {
      StaticArray tmp ;
      for ( unsigned i = 0 ; i < N ; ++i ) tmp.value[i] = value[i] * rhs.value[i] ;
      return tmp ;
    }

  KOKKOS_INLINE_FUNCTION
  StaticArray operator + ( const StaticArray & rhs )
    {
      StaticArray tmp ;
      for ( unsigned i = 0 ; i < N ; ++i ) tmp.value[i] = value[i] + rhs.value[i] ;
      return tmp ;
    }

  KOKKOS_INLINE_FUNCTION
  StaticArray & operator += ( const StaticArray & rhs )
    {
      for ( unsigned i = 0 ; i < N ; ++i ) value[i] += rhs.value[i] ;
      return *this ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator += ( const volatile StaticArray & rhs ) volatile
    {
      for ( unsigned i = 0 ; i < N ; ++i ) value[i] += rhs.value[i] ;
    }
};

template< typename T , class Space >
struct DOT {
  typedef T      value_type ;
  typedef Space  execution_space ;

  Kokkos::View< value_type * , Space > a ;
  Kokkos::View< value_type * , Space > b ;

  DOT( const Kokkos::View< value_type * , Space > arg_a
     , const Kokkos::View< value_type * , Space > arg_b
     )
    : a( arg_a ), b( arg_b ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i , value_type & update ) const
    {
      update += a(i) * b(i);
    }
};

template< typename T , class Space >
struct FILL {
  typedef T      value_type ;
  typedef Space  execution_space ;

  Kokkos::View< value_type * , Space > a ;
  Kokkos::View< value_type * , Space > b ;

  FILL( const Kokkos::View< value_type * , Space > & arg_a
      , const Kokkos::View< value_type * , Space > & arg_b
      )
    : a( arg_a ), b( arg_b ) {}

  KOKKOS_INLINE_FUNCTION
  void operator()( const int i ) const
    {
      a(i) = i % 2 ? i + 1 : 1 ;
      b(i) = i % 2 ? 1 : i + 1 ;
    }
};

template< class Space >
void TestViewAggregateReduction()
{
  const int count = 2 ;
  const long result = count % 2 ? ( count * ( ( count + 1 ) / 2 ) )
                                : ( ( count / 2 ) * ( count + 1 ) );

  Kokkos::View< long * , Space > a("a",count);
  Kokkos::View< long * , Space > b("b",count);
  Kokkos::View< StaticArray<long,4> * , Space > a4("a4",count);
  Kokkos::View< StaticArray<long,4> * , Space > b4("b4",count);
  Kokkos::View< StaticArray<long,10> * , Space > a10("a10",count);
  Kokkos::View< StaticArray<long,10> * , Space > b10("b10",count);

  Kokkos::parallel_for( count , FILL<long,Space>(a,b) );
  Kokkos::parallel_for( count , FILL< StaticArray<long,4> , Space >(a4,b4) );
  Kokkos::parallel_for( count , FILL< StaticArray<long,10> , Space >(a10,b10) );

  long r = 0;
  StaticArray<long,4> r4 ;
  StaticArray<long,10> r10 ;

  Kokkos::parallel_reduce( count , DOT<long,Space>(a,b) , r );
  Kokkos::parallel_reduce( count , DOT< StaticArray<long,4> , Space >(a4,b4) , r4 );
  Kokkos::parallel_reduce( count , DOT< StaticArray<long,10> , Space >(a10,b10) , r10 );

  ASSERT_EQ( result , r );
  for ( int i = 0 ; i < 10 ; ++i ) { ASSERT_EQ( result , r10.value[i] ); }
  for ( int i = 0 ; i < 4 ; ++i ) { ASSERT_EQ( result , r4.value[i] ); }
}

}

#endif /* #ifndef TEST_AGGREGATE_REDUCTION_HPP */

