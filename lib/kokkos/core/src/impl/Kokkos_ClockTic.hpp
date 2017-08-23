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

#ifndef KOKKOS_CLOCKTIC_HPP
#define KOKKOS_CLOCKTIC_HPP

#include <Kokkos_Macros.hpp>
#include <stdint.h>
#include <chrono>

namespace Kokkos {
namespace Impl {

/**\brief  Quick query of clock register tics
 *
 *  Primary use case is to, with low overhead,
 *  obtain a integral value that consistently varies
 *  across concurrent threads of execution within
 *  a parallel algorithm.
 *  This value is often used to "randomly" seed an
 *  attempt to acquire an indexed resource (e.g., bit)
 *  from an array of resources (e.g., bitset) such that
 *  concurrent threads will have high likelihood of
 *  having different index-seed values.
 */
KOKKOS_FORCEINLINE_FUNCTION
uint64_t clock_tic(void) noexcept
{
#if defined( __CUDA_ARCH__ )

  // Return value of 64-bit hi-res clock register.

  return clock64();

#elif defined(__HCC_ACCELERATOR__)
    // Get clock register
    return hc::__clock_u64();

#elif defined( __i386__ ) || defined( __x86_64 )

  // Return value of 64-bit hi-res clock register.

  unsigned a = 0, d = 0;

  __asm__ volatile( "rdtsc" : "=a" (a), "=d" (d) );

  return ( (uint64_t) a ) | ( ( (uint64_t) d ) << 32 );

#elif defined( __powerpc )     || defined( __powerpc__ ) || \
      defined( __powerpc64__ ) || defined( __POWERPC__ ) || \
      defined( __ppc__ )       || defined( __ppc64__ )

  unsigned int cycles = 0;

  asm volatile( "mftb %0" : "=r" (cycles) );

  return (uint64_t) cycles;

#else

  return (uint64_t)
    std::chrono::high_resolution_clock::now().time_since_epoch().count();

#endif
}

} // namespace Impl
} // namespace Kokkos

#endif // KOKKOS_CLOCKTIC_HPP
