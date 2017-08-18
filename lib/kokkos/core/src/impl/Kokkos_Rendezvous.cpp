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
#include <impl/Kokkos_Rendezvous.hpp>
#include <impl/Kokkos_Spinwait.hpp>

namespace Kokkos { namespace Impl {

//----------------------------------------------------------------------------
/* pattern for rendezvous
 *
 *  if ( rendezvous() ) {
 *     ... all other threads are still in team_rendezvous() ...
 *     rendezvous_release();
 *     ... all other threads are released from team_rendezvous() ...
 *  }
 */

int rendezvous( volatile int64_t * const buffer
              , int const size
              , int const rank
              , int const slow
              ) noexcept
{
  enum : int { shift_byte = 3 };
  enum : int { size_byte  = ( 01 << shift_byte ) }; // == 8
  enum : int { mask_byte  = size_byte - 1 };

  enum : int { shift_mem_cycle = 2 };
  enum : int { size_mem_cycle  = ( 01 << shift_mem_cycle ) }; // == 4
  enum : int { mask_mem_cycle  = size_mem_cycle - 1 };

  // Cycle step values: 1 <= step <= size_val_cycle
  // An odd multiple of memory cycle so that when a memory location
  // is reused it has a different value.
  // Must be representable within a single byte: size_val_cycle < 16

  enum : int { size_val_cycle = 3 * size_mem_cycle };

  // Requires:
  //   Called by rank = [ 0 .. size )
  //   buffer aligned to int64_t[4]

  // A sequence of rendezvous uses four cycled locations in memory
  // and non-equal cycled synchronization values to
  // 1) prevent rendezvous from overtaking one another and
  // 2) give each spin wait location an int64_t[4] span
  //    so that it has its own cache line.

  const int64_t step = (buffer[0] % size_val_cycle ) + 1 ;

  // The leading int64_t[4] span is for thread 0 to write
  // and all other threads to read spin-wait.
  // sync_offset is the index into this array for this step.

  const int sync_offset = ( step & mask_mem_cycle ) + size_mem_cycle + size_mem_cycle ;

  if ( rank ) {

    const int group_begin = rank << shift_byte ; // == rank * size_byte

    if ( group_begin < size ) {

      //  This thread waits for threads
      //   [ group_begin .. group_begin + 8 )
      //   [ rank*8      .. rank*8 + 8      )
      // to write to their designated bytes.

      const int end = group_begin + size_byte < size
                    ? size_byte : size - group_begin ;

      int64_t value = 0;
      for ( int i = 0 ; i < end ; ++i ) {
        value |= step << (i * size_byte );
      }

      store_fence(); // This should not be needed but fixes #742

      if ( slow ) {
        yield_until_equal( buffer[ (rank << shift_mem_cycle) + sync_offset ]
                          , value );
      }
      else {
        spinwait_until_equal( buffer[ (rank << shift_mem_cycle) + sync_offset ]
                            , value );
      }
    }

    {
      // This thread sets its designated byte.
      //   ( rank % size_byte ) +
      //   ( ( rank / size_byte ) * size_byte * size_mem_cycle ) +
      //   ( sync_offset * size_byte )
      const int offset = ( rank & mask_byte )
                       + ( ( rank & ~mask_byte ) << shift_mem_cycle )
                       + ( sync_offset << shift_byte );

      // All of this thread's previous memory stores must be complete before
      // this thread stores the step value at this thread's designated byte
      // in the shared synchronization array.

      Kokkos::memory_fence();

      ((volatile int8_t*) buffer)[ offset ] = int8_t( step );

      // Memory fence to push the previous store out
      Kokkos::memory_fence();
    }

    // Wait for thread 0 to release all other threads

    if ( slow ) {
      yield_until_equal( buffer[ (step & mask_mem_cycle) + size_mem_cycle ] , int64_t(step) );
    }
    else {
      spinwait_until_equal( buffer[ (step & mask_mem_cycle) + size_mem_cycle ] , int64_t(step) );
    }
  }
  else {
    // Thread 0 waits for threads [1..7]
    // to write to their designated bytes.

    const int end = size_byte < size ? 8 : size ;

    int64_t value = 0;
    for ( int i = 1 ; i < end ; ++i ) {
      value |= step << (i * size_byte );
    }

    if ( slow ) {
      yield_until_equal( buffer[ sync_offset ], value );
    }
    else {
      spinwait_until_equal( buffer[ sync_offset ], value );
    }
  }

  return rank ? 0 : 1 ;
}

void rendezvous_release( volatile int64_t * const buffer ) noexcept
{
  enum : int { shift_mem_cycle = 2 };
  enum : int { size_mem_cycle  = ( 01 << shift_mem_cycle ) }; // == 4
  enum : int { mask_mem_cycle  = size_mem_cycle - 1 };
  enum : int { size_val_cycle = 3 * size_mem_cycle };

  // Requires:
  //   Called after team_rendezvous
  //   Called only by true == team_rendezvous(root)

  // update step
  const int64_t step = (buffer[0] % size_val_cycle ) + 1;
  buffer[0] = step;

  // Memory fence to be sure all previous writes are complete:
  Kokkos::memory_fence();

  buffer[ (step & mask_mem_cycle) + size_mem_cycle ] = step;

  // Memory fence to push the store out
  Kokkos::memory_fence();
}

}} // namespace Kokkos::Impl

