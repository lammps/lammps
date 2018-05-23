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

#ifndef TEST_CONCURRENTBITSET_HPP
#define TEST_CONCURRENTBITSET_HPP

#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <impl/Kokkos_ConcurrentBitset.hpp>

namespace Test {

template< class DeviceType >
struct ConcurrentBitset {

  typedef Kokkos::View<uint32_t*,DeviceType> view_unsigned_type ;
  typedef Kokkos::View<int*,DeviceType>      view_int_type ;

  view_unsigned_type  bitset ;
  view_int_type       acquired ;
  uint32_t            bitset_count_lg2 ;
  uint32_t            bitset_count_mask ;

  ConcurrentBitset( const uint32_t arg_bitset_count_lg2
                  , const view_unsigned_type & arg_bitset
                  , const view_int_type & arg_acquired )
    : bitset( arg_bitset ), acquired( arg_acquired )
    , bitset_count_lg2( arg_bitset_count_lg2 )
    , bitset_count_mask( uint32_t( 1u << arg_bitset_count_lg2 ) - 1 )
    {}

  struct TagAcquire {};
  struct TagRelease {};
  struct TagReacquire {};

  KOKKOS_INLINE_FUNCTION
  void operator()( TagAcquire , int i , long & update ) const
    {
      unsigned hint = Kokkos::Impl::clock_tic() & bitset_count_mask ;

      Kokkos::pair<int,int> result =
        Kokkos::Impl::concurrent_bitset::acquire_bounded_lg2
          ( bitset.data() , bitset_count_lg2 , hint );

      acquired(i) = result.first ;

      if ( 0 <= result.first ) ++update ;
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( TagRelease , int i , long & update ) const
    {
      if ( 0 == ( i % 3 ) && 0 <= acquired(i) ) {
        Kokkos::Impl::concurrent_bitset::release( bitset.data() , acquired(i) );
        acquired(i) = -1 ;
        ++update ;
      }
    }

  KOKKOS_INLINE_FUNCTION
  void operator()( TagReacquire , int i , long & update ) const
    {
      if ( acquired(i) < 0 ) {

        unsigned hint = Kokkos::Impl::clock_tic() & bitset_count_mask ;

        Kokkos::pair<int,int> result  = Kokkos::Impl::concurrent_bitset::acquire_bounded_lg2
            ( bitset.data() , bitset_count_lg2 , hint );

        acquired(i) = result.first ;

        if ( 0 <= result.first ) ++update ;
      }
    }
};

template< class DeviceType >
void test_concurrent_bitset( int bit_count )
{
  typedef ConcurrentBitset< DeviceType > Functor ;
  typedef typename Functor::view_unsigned_type view_unsigned_type ;
  typedef typename Functor::view_int_type      view_int_type ;

  int bit_count_lg2 = 1 ;

  while ( ( 1 << bit_count_lg2 ) < bit_count ) ++bit_count_lg2 ;

  bit_count = 1 << bit_count_lg2 ;

  const int buffer_length =
    Kokkos::Impl::concurrent_bitset::buffer_bound_lg2(bit_count_lg2);

  view_unsigned_type bitset("bitset",buffer_length);

  // Try to acquire more than available:

  const size_t n = ( bit_count * 3 ) / 2 ;

  view_int_type acquired("acquired", n );

  typename view_unsigned_type::HostMirror bitset_host =
    Kokkos::create_mirror_view( bitset );

  Kokkos::deep_copy( bitset , 0u );

  long total = 0 ;
  long total_release = 0 ;
  long total_reacquire = 0 ;

  Kokkos::parallel_reduce
    ( Kokkos::RangePolicy< DeviceType , typename Functor::TagAcquire >(0,n)
    , Functor( bit_count_lg2 , bitset , acquired )
    , total );

  ASSERT_EQ( bit_count , total );

  Kokkos::parallel_reduce
    ( Kokkos::RangePolicy< DeviceType , typename Functor::TagRelease >(0,n)
    , Functor( bit_count_lg2 , bitset , acquired )
    , total_release );

  Kokkos::parallel_reduce
    ( Kokkos::RangePolicy< DeviceType , typename Functor::TagReacquire >(0,n)
    , Functor( bit_count_lg2 , bitset , acquired )
    , total_reacquire );

  ASSERT_EQ( total_release , total_reacquire );

}

} // namespace Test

#endif /* #ifndef TEST_CONCURRENTBITSET_HPP */
