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

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ATOMIC_HPP ) && ! defined( KOKKOS_ATOMIC_ASSEMBLY_HPP )
#define KOKKOS_ATOMIC_ASSEMBLY_HPP
namespace Kokkos {

namespace Impl {
  struct cas128_t
  {
    uint64_t lower;
    uint64_t upper;

    KOKKOS_INLINE_FUNCTION
    cas128_t () {
      lower = 0;
      upper = 0;
    }

    KOKKOS_INLINE_FUNCTION
    cas128_t (const cas128_t& a) {
      lower = a.lower;
      upper = a.upper;
    }
    KOKKOS_INLINE_FUNCTION
    cas128_t (volatile cas128_t* a) {
      lower = a->lower;
      upper = a->upper;
    }

    KOKKOS_INLINE_FUNCTION
    bool operator != (const cas128_t& a) const {
      return (lower != a.lower) || upper!=a.upper;
    }

    KOKKOS_INLINE_FUNCTION
    void operator = (const cas128_t& a) {
      lower = a.lower;
      upper = a.upper;
    }
    KOKKOS_INLINE_FUNCTION
    void operator = (const cas128_t& a) volatile {
      lower = a.lower;
      upper = a.upper;
    }
  }
  __attribute__ (( __aligned__( 16 ) ));


  #if defined( KOKKOS_ENABLE_ASM ) && defined ( KOKKOS_ENABLE_ISA_X86_64 )
  inline cas128_t cas128( volatile cas128_t * ptr, cas128_t cmp,  cas128_t swap )
  {
      bool swapped = false;
      __asm__ __volatile__
      (
       "lock cmpxchg16b %1\n\t"
       "setz %0"
       : "=q" ( swapped )
       , "+m" ( *ptr )
       , "+d" ( cmp.upper )
       , "+a" ( cmp.lower )
       : "c" ( swap.upper )
       , "b" ( swap.lower )
       , "q" ( swapped )
     );
      return cmp;
  }
  #endif

}
}

#endif

