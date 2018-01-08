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
#include <Kokkos_Atomic.hpp>

#include <impl/Kokkos_HostBarrier.hpp>
#include <impl/Kokkos_Spinwait.hpp>

namespace Kokkos { namespace Impl {

namespace {

enum : int { HEADER_SIZE = HostBarrier::HEADER / sizeof(uint64_t) };

inline constexpr int length64( const int nthreads ) noexcept
{
  return (nthreads-1 + sizeof(uint64_t)-1) / sizeof(uint64_t);
}

} // namespace

void rendezvous_initialize( volatile void * buffer
                          , const int size
                          , const int rank
                          ) noexcept
{
  Kokkos::store_fence();

  // ensure that the buffer has been zero'd out
  constexpr uint8_t  zero8  = static_cast<uint8_t>(0);
  constexpr uint64_t zero64 = static_cast<uint64_t>(0);

  volatile uint64_t * header = reinterpret_cast<volatile uint64_t *>(buffer);

  if (rank > 0) {
    volatile uint8_t * bytes = reinterpret_cast<volatile uint8_t *>(buffer) + RENDEZVOUS_HEADER;

    bytes[rank-1] = zero8;

    // last thread is responsible for zeroing out the final bytes of the last uint64_t
    if (rank == size-1) {
      const int tmp  = (size-1) % sizeof(uint64_t);
      const int rem = tmp ? sizeof(uint64_t) - tmp : 0;
      for (int i=0; i<rem; ++i) {
        bytes[rank+i] = zero8;
      }
    }

    spinwait_until_equal( *header, zero64 );
  }
  else {

    const int n = length64(size);
    volatile uint64_t * buff = reinterpret_cast<volatile uint64_t *>(buffer) + HEADER_SIZE;

    // wait for other threads to finish initializing
    for (int i=0; i<n; ++i) {
      spinwait_until_equal( buff[i], zero64 );
    }

    // release the waiting threads
    *header = zero64;
    Kokkos::store_fence();
  }
  Kokkos::load_fence();
}

bool rendezvous( volatile void * buffer
               , uint64_t &      step
               , const int       size
               , const int       rank
               , bool            active_wait
               ) noexcept
{
  // Force all outstanding stores from this thread to retire before continuing
  Kokkos::store_fence();

  // guarantees that will never spinwait on a spin_value of 0
  step = static_cast<uint8_t>(step + 1u)
         ? step + 1u
         : step + 2u
         ;

  // if size == 1, it is incorrect for rank 0 to check the tail value of the buffer
  // this optimization prevents a potential read of uninitialized memory
  if ( size == 1 ) { return true; }

  const uint8_t byte_value  = static_cast<uint8_t>(step);

  // byte that is set in the spin_value rotates every time
  // this prevents threads from overtaking the master thread
  const uint64_t spin_value = static_cast<uint64_t>(byte_value) << (byte_value&7);

  if ( rank > 0 ) {
    volatile uint64_t * header = reinterpret_cast<volatile uint64_t *>(buffer);
    volatile uint8_t *  bytes  = reinterpret_cast<volatile uint8_t *>(buffer) + RENDEZVOUS_HEADER;

    bytes[ rank-1 ] = byte_value;

    if ( active_wait ) {
      spinwait_until_equal( *header, spin_value );
    }
    else {
      yield_until_equal( *header, spin_value );
    }
  }
  else { // rank 0
    volatile uint64_t * buff = reinterpret_cast<volatile uint64_t *>(buffer) + HEADER_SIZE;
    const int n = length64(size);

    uint64_t comp = byte_value;
    comp = comp | (comp << 8);
    comp = comp | (comp << 16);
    comp = comp | (comp << 32);

    const int rem  = (size-1) % sizeof(uint64_t);

    union {
      volatile uint64_t value;
      volatile uint8_t  array[sizeof(uint64_t)];
    } tmp{};

    for (int i=0; i<rem; ++i) {
      tmp.array[i] = byte_value;
    }

    const uint64_t tail = rem ? tmp.value : comp;

    for (int i=0; i<n-1; ++i) {
      spinwait_until_equal( buff[i], comp );
    }
    spinwait_until_equal( buff[n-1], tail );

  }

  // Force all outstanding stores from other threads to retire before allowing
  // this thread to continue.  This forces correctness on systems with out-of-order
  // memory (Power and ARM)
  Kokkos::load_fence();

  return rank == 0;
}

void rendezvous_release( volatile void * buffer
                       , const uint64_t  step
                       ) noexcept
{
  const uint8_t       byte_value = static_cast<uint8_t>(step);
  const uint64_t      spin_value = static_cast<uint64_t>(byte_value) << (byte_value&7);
  volatile uint64_t * header     = reinterpret_cast<volatile uint64_t *>(buffer);

  // Force all outstanding stores from this thread to retire before releasing
  // the other threads.  This forces correctness on systems with out-of-order
  // memory (Power and ARM)
  Kokkos::store_fence();

  *header = spin_value;

  Kokkos::memory_fence();
}

}} // namespace Kokkos::Impl

