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

#ifndef KOKKOS_BITOPS_HPP
#define KOKKOS_BITOPS_HPP

#include <Kokkos_Macros.hpp>
#include <cstdint>
#include <climits>

#ifdef KOKKOS_COMPILER_INTEL
#include<immintrin.h>
#endif

#if defined( __HCC_ACCELERATOR__ )
#include <hc.hpp>
#endif

namespace Kokkos {

KOKKOS_FORCEINLINE_FUNCTION
int log2( unsigned i )
{
  enum : int { shift = sizeof(unsigned) * CHAR_BIT - 1 };
#if defined( __CUDA_ARCH__ )
  return shift - __clz(i);
#elif defined( __HCC_ACCELERATOR__ )
  return  (int)hc::__firstbit_u32_u32(i);
#elif defined( KOKKOS_COMPILER_INTEL )
  return _bit_scan_reverse(i);
#elif defined( KOKKOS_COMPILER_IBM )
  return shift - __cntlz4(i);
#elif defined( KOKKOS_COMPILER_CRAYC )
  return i ? shift - _leadz32(i) : 0 ;
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return shift - __builtin_clz(i);
#else
  int offset = 0;
  if ( i ) {
    for ( offset = shift ; (i & ( 1 << offset ) ) == 0 ; --offset );
  }
  return offset;
#endif
}

namespace Impl {

/**\brief  Find first zero bit.
 *
 *  If none then return -1 ;
 */
KOKKOS_FORCEINLINE_FUNCTION
int bit_first_zero( unsigned i ) noexcept
{
  enum : unsigned { full = ~0u };

#if defined( __CUDA_ARCH__ )
  return full != i ? __ffs( ~i ) - 1 : -1 ;
#elif defined( __HCC_ACCELERATOR__ )
  return full != i ? (int)hc::__firstbit_u32_u32(~i) : -1 ;
#elif defined( KOKKOS_COMPILER_INTEL )
  return full != i ? _bit_scan_forward( ~i ) : -1 ;
#elif defined( KOKKOS_COMPILER_IBM )
  return full != i ? __cnttz4( ~i ) : -1 ;
#elif defined( KOKKOS_COMPILER_CRAYC )
  return full != i ? _popcnt( i ^ (i+1) ) - 1 : -1 ;
#elif defined( KOKKOS_COMPILER_GNU ) || defined( __GNUC__ ) || defined( __GNUG__ )
  return full != i ? __builtin_ffs( ~i ) - 1 : -1 ;
#else
  int offset = -1 ;
  if ( full != i ) {
    for ( offset = 0 ; i & ( 1 << offset ) ; ++offset );
  }
  return offset ;
#endif
}

KOKKOS_FORCEINLINE_FUNCTION
int bit_scan_forward( unsigned i )
{
#if defined( __CUDA_ARCH__ )
  return __ffs(i) - 1;
#elif defined( __HCC_ACCELERATOR__ )
  return  (int)hc::__firstbit_u32_u32(i);
#elif defined( KOKKOS_COMPILER_INTEL )
  return _bit_scan_forward(i);
#elif defined( KOKKOS_COMPILER_IBM )
  return __cnttz4(i);
#elif defined( KOKKOS_COMPILER_CRAYC )
  return i ? _popcnt(~i & (i-1)) : -1;
#elif defined( KOKKOS_COMPILER_GNU ) || defined( __GNUC__ ) || defined( __GNUG__ )
  return __builtin_ffs(i) - 1;
#else
  int offset = -1;
  if ( i ) {
    for ( offset = 0 ; (i & ( 1 << offset ) ) == 0 ; ++offset );
  }
  return offset;
#endif
}

/// Count the number of bits set.
KOKKOS_FORCEINLINE_FUNCTION
int bit_count( unsigned i )
{
#if defined( __CUDA_ARCH__ )
  return __popc(i);
#elif defined( __HCC_ACCELERATOR__ )
  return  (int)hc::__popcount_u32_b32(i);
#elif defined ( __INTEL_COMPILER )
  return _popcnt32(i);
#elif defined( KOKKOS_COMPILER_IBM )
  return __popcnt4(i);
#elif defined( KOKKOS_COMPILER_CRAYC )
  return _popcnt(i);
#elif defined( __GNUC__ ) || defined( __GNUG__ )
  return __builtin_popcount(i);
#else
  // http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetNaive
  i = i - ( ( i >> 1 ) & ~0u / 3u );                             // temp
  i = ( i & ~0u / 15u * 3u ) + ( ( i >> 2 ) & ~0u / 15u * 3u );  // temp
  i = ( i + ( i >> 4 ) ) & ~0u / 255u * 15u;                     // temp

  // count
  return (int)( ( i * ( ~0u / 255u ) ) >> ( sizeof(unsigned) - 1 ) * CHAR_BIT );
#endif
}

KOKKOS_INLINE_FUNCTION
unsigned integral_power_of_two_that_contains( const unsigned N )
{
  const unsigned i = Kokkos::log2( N );
  return ( (1u << i) < N ) ? i + 1 : i ;
}


} // namespace Impl
} // namespace Kokkos

#endif // KOKKOS_BITOPS_HPP

