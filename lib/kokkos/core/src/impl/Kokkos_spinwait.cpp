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
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )

#include <impl/Kokkos_spinwait.hpp>

#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_BitOps.hpp>

/*--------------------------------------------------------------------------*/

#if !defined( _WIN32 )
  #if defined( KOKKOS_ENABLE_ASM )
    #if defined( __arm__ ) || defined( __aarch64__ )
      /* No-operation instruction to idle the thread. */
      #define KOKKOS_INTERNAL_PAUSE
    #else
      /* Pause instruction to prevent excess processor bus usage */
      #define KOKKOS_INTERNAL_PAUSE   asm volatile("pause\n":::"memory")
    #endif
    #define KOKKOS_INTERNAL_NOP2    asm volatile("nop\n" "nop\n")
    #define KOKKOS_INTERNAL_NOP4    KOKKOS_INTERNAL_NOP2;  KOKKOS_INTERNAL_NOP2
    #define KOKKOS_INTERNAL_NOP8    KOKKOS_INTERNAL_NOP4;  KOKKOS_INTERNAL_NOP4;
    #define KOKKOS_INTERNAL_NOP16   KOKKOS_INTERNAL_NOP8;  KOKKOS_INTERNAL_NOP8;
    #define KOKKOS_INTERNAL_NOP32   KOKKOS_INTERNAL_NOP16; KOKKOS_INTERNAL_NOP16;
    namespace {
    inline void kokkos_internal_yield( const unsigned i ) noexcept {
      switch (Kokkos::Impl::bit_scan_reverse((i >> 2)+1u)) {
      case 0u:  KOKKOS_INTERNAL_NOP2;  break;
      case 1u:  KOKKOS_INTERNAL_NOP4;  break;
      case 2u:  KOKKOS_INTERNAL_NOP8;  break;
      case 3u:  KOKKOS_INTERNAL_NOP16; break;
      default: KOKKOS_INTERNAL_NOP32;
      }
      KOKKOS_INTERNAL_PAUSE;
    }
    }
  #else
    #include <sched.h>
    namespace {
    inline void kokkos_internal_yield( const unsigned ) noexcept {
      sched_yield();
    }
    }
  #endif
#else // defined( _WIN32 )
  #if defined ( KOKKOS_ENABLE_WINTHREAD )
    #include <process.h>
    namespace {
    inline void kokkos_internal_yield( const unsigned ) noexcept {
      Sleep(0);
    }
    }
  #elif defined( _MSC_VER )
    #define NOMINMAX
    #include <winsock2.h>
    #include <windows.h>
    namespace {
    inline void kokkos_internal_yield( const unsigned ) noexcept {
      YieldProcessor();
    }
    }
  #else
    #define KOKKOS_INTERNAL_PAUSE   __asm__ __volatile__("pause\n":::"memory")
    #define KOKKOS_INTERNAL_NOP2    __asm__ __volatile__("nop\n" "nop")
    #define KOKKOS_INTERNAL_NOP4    KOKKOS_INTERNAL_NOP2;  KOKKOS_INTERNAL_NOP2
    #define KOKKOS_INTERNAL_NOP8    KOKKOS_INTERNAL_NOP4;  KOKKOS_INTERNAL_NOP4;
    #define KOKKOS_INTERNAL_NOP16   KOKKOS_INTERNAL_NOP8;  KOKKOS_INTERNAL_NOP8;
    #define KOKKOS_INTERNAL_NOP32   KOKKOS_INTERNAL_NOP16; KOKKOS_INTERNAL_NOP16;
    namespace {
    inline void kokkos_internal_yield( const unsigned i ) noexcept {
      switch (Kokkos::Impl::bit_scan_reverse((i >> 2)+1u)) {
      case 0:  KOKKOS_INTERNAL_NOP2;  break;
      case 1:  KOKKOS_INTERNAL_NOP4;  break;
      case 2:  KOKKOS_INTERNAL_NOP8;  break;
      case 3:  KOKKOS_INTERNAL_NOP16; break;
      default: KOKKOS_INTERNAL_NOP32;
      }
      KOKKOS_INTERNAL_PAUSE;
    }
    }
  #endif
#endif


/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

void spinwait_while_equal( volatile int32_t & flag , const int32_t value )
{
  Kokkos::store_fence();
  unsigned i = 0;
  while ( value == flag ) {
    kokkos_internal_yield(i);
    ++i;
  }
  Kokkos::load_fence();
}

void spinwait_until_equal( volatile int32_t & flag , const int32_t value )
{
  Kokkos::store_fence();
  unsigned i = 0;
  while ( value != flag ) {
    kokkos_internal_yield(i);
    ++i;
  }
  Kokkos::load_fence();
}

void spinwait_while_equal( volatile int64_t & flag , const int64_t value )
{
  Kokkos::store_fence();
  unsigned i = 0;
  while ( value == flag ) {
    kokkos_internal_yield(i);
    ++i;
  }
  Kokkos::load_fence();
}

void spinwait_until_equal( volatile int64_t & flag , const int64_t value )
{
  Kokkos::store_fence();
  unsigned i = 0;
  while ( value != flag ) {
    kokkos_internal_yield(i);
    ++i;
  }
  Kokkos::load_fence();
}

} /* namespace Impl */
} /* namespace Kokkos */

#else
void KOKKOS_CORE_SRC_IMPL_SPINWAIT_PREVENT_LINK_ERROR() {}
#endif

