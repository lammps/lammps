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

#ifndef KOKKOS_BITSET_IMPL_HPP
#define KOKKOS_BITSET_IMPL_HPP

#include <Kokkos_Macros.hpp>
#include <stdint.h>

#include <cstdio>
#include <climits>
#include <iostream>
#include <iomanip>

namespace Kokkos { namespace Impl {

KOKKOS_FORCEINLINE_FUNCTION
unsigned rotate_right(unsigned i, int r)
{
  enum { size = static_cast<int>(sizeof(unsigned)*CHAR_BIT) };
  return r ? ((i >> r) | (i << (size-r))) : i ;
}

KOKKOS_FORCEINLINE_FUNCTION
int bit_scan_forward(unsigned i)
{
#if defined( __CUDA_ARCH__ )
  return __ffs(i) - 1;
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return __builtin_ffs(i) - 1;
#elif defined( __INTEL_COMPILER )
  return _bit_scan_forward(i);
#else

  unsigned t = 1u;
  int r = 0;
  while (i && (i & t == 0))
  {
    t = t << 1;
    ++r;
  }
  return r;
#endif
}


KOKKOS_FORCEINLINE_FUNCTION
int bit_scan_reverse(unsigned i)
{
  enum { shift = static_cast<int>(sizeof(unsigned)*CHAR_BIT - 1) };
#if defined( __CUDA_ARCH__ )
  return shift - __clz(i);
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return shift - __builtin_clz(i);
#elif defined( __INTEL_COMPILER )
  return _bit_scan_reverse(i);
#else
  unsigned t = 1u << shift;
  int r = 0;
  while (i && (i & t == 0))
  {
    t = t >> 1;
    ++r;
  }
  return r;
#endif
}


// count the bits set
KOKKOS_FORCEINLINE_FUNCTION
int popcount(unsigned i)
{
#if defined( __CUDA_ARCH__ )
  return __popc(i);
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return __builtin_popcount(i);
#elif defined ( __INTEL_COMPILER )
  return _popcnt32(i);
#else
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive
  i = i - ((i >> 1) & ~0u/3u);                                         // temp
  i = (i & ~0u/15u*3u) + ((i >> 2) & ~0u/15u*3u);                      // temp
  i = (i + (i >> 4)) & ~0u/255u*15u;                                   // temp
  return (int)((i * (~0u/255u)) >> (sizeof(unsigned) - 1) * CHAR_BIT); // count
#endif
}


template <typename Bitset>
struct BitsetCount
{
  typedef Bitset bitset_type;
  typedef typename bitset_type::device_type device_type;
  typedef typename bitset_type::size_type size_type;
  typedef size_type value_type;

  bitset_type m_bitset;

  BitsetCount( bitset_type const& bitset)
    : m_bitset(bitset)
  {}

  size_type apply() const
  {
    size_type count = 0u;
    parallel_reduce(m_bitset.m_blocks.size(), *this, count);
    return count;
  }

  KOKKOS_INLINE_FUNCTION
  static void init( value_type & count)
  {
    count = 0u;
  }

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & count, const volatile size_type & incr )
  {
    count += incr;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( size_type i, value_type & count) const
  {
    count += popcount(m_bitset.m_blocks[i]);
  }
};

}} //Kokkos::Impl

#endif // KOKKOS_BITSET_IMPL_HPP

