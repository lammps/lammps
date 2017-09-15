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

#ifndef KOKKOS_FUNCTIONAL_IMPL_HPP
#define KOKKOS_FUNCTIONAL_IMPL_HPP

#include <Kokkos_Macros.hpp>
#include <cstdint>

namespace Kokkos { namespace Impl {

// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.
KOKKOS_FORCEINLINE_FUNCTION
uint32_t getblock32 ( const uint8_t * p, int i )
{
// used to avoid aliasing error which could cause errors with
// forced inlining
  return    ((uint32_t)p[i*4+0])
          | ((uint32_t)p[i*4+1] << 8)
          | ((uint32_t)p[i*4+2] << 16)
          | ((uint32_t)p[i*4+3] << 24);
}

KOKKOS_FORCEINLINE_FUNCTION
uint32_t rotl32 ( uint32_t x, int8_t r )
{ return (x << r) | (x >> (32 - r)); }

KOKKOS_FORCEINLINE_FUNCTION
uint32_t fmix32 ( uint32_t h )
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;

  return h;
}

KOKKOS_INLINE_FUNCTION
uint32_t MurmurHash3_x86_32 ( const void * key, int len, uint32_t seed )
{
  const uint8_t * data = (const uint8_t*)key;
  const int nblocks = len / 4;

  uint32_t h1 = seed;

  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;

  //----------
  // body

  for(int i=0; i<nblocks; ++i)
  {
    uint32_t k1 = getblock32(data,i);

    k1 *= c1;
    k1 = rotl32(k1,15);
    k1 *= c2;

    h1 ^= k1;
    h1 = rotl32(h1,13);
    h1 = h1*5+0xe6546b64;
  }

  //----------
  // tail

  const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

  uint32_t k1 = 0;

  switch(len & 3)
  {
  case 3: k1 ^= tail[2] << 16;
  case 2: k1 ^= tail[1] << 8;
  case 1: k1 ^= tail[0];
          k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
  };

  //----------
  // finalization

  h1 ^= len;

  h1 = fmix32(h1);

  return h1;
}


#if defined( __GNUC__ ) /* GNU C   */ || \
    defined( __GNUG__ ) /* GNU C++ */ || \
    defined( __clang__ )

#define KOKKOS_IMPL_MAY_ALIAS __attribute__((__may_alias__))

#else

#define KOKKOS_IMPL_MAY_ALIAS

#endif

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
bool bitwise_equal(T const * const a_ptr, T const * const b_ptr)
{
  typedef uint64_t KOKKOS_IMPL_MAY_ALIAS T64;
  typedef uint32_t KOKKOS_IMPL_MAY_ALIAS T32;
  typedef uint16_t KOKKOS_IMPL_MAY_ALIAS T16;
  typedef uint8_t  KOKKOS_IMPL_MAY_ALIAS T8;

  enum {
    NUM_8  = sizeof(T),
    NUM_16 = NUM_8 / 2,
    NUM_32 = NUM_8 / 4,
    NUM_64 = NUM_8 / 8
  };

  union {
    T   const * const ptr;
    T64 const * const ptr64;
    T32 const * const ptr32;
    T16 const * const ptr16;
    T8  const * const ptr8;
  } a = {a_ptr}, b = {b_ptr};

  bool result = true;

  for (int i=0; i < NUM_64; ++i) {
    result = result && a.ptr64[i] == b.ptr64[i];
  }

  if ( NUM_64*2 < NUM_32 ) {
    result = result && a.ptr32[NUM_64*2] == b.ptr32[NUM_64*2];
  }

  if ( NUM_32*2 < NUM_16 ) {
    result = result && a.ptr16[NUM_32*2] == b.ptr16[NUM_32*2];
  }

  if ( NUM_16*2 < NUM_8 ) {
    result = result && a.ptr8[NUM_16*2] == b.ptr8[NUM_16*2];
  }

  return result;
}



#undef KOKKOS_IMPL_MAY_ALIAS

}} // namespace Kokkos::Impl

#endif //KOKKOS_FUNCTIONAL_IMPL_HPP

