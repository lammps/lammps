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

#include <limits>
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Spinwait.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void HostThreadTeamData::organize_pool
  ( HostThreadTeamData * members[] , const int size )
{
  bool ok = true ;

  memory_fence();

  // Verify not already a member of a pool:
  for ( int rank = 0 ; rank < size && ok ; ++rank ) {
    ok = ( nullptr != members[rank] ) && ( 0 == members[rank]->m_pool_scratch );
  }

  if ( ok ) {

    int64_t * const root_scratch = members[0]->m_scratch ;

    for ( int i = m_pool_rendezvous ; i < m_pool_reduce ; ++i ) {
      root_scratch[i] = 0 ;
    }

    {
      HostThreadTeamData ** const pool =
        (HostThreadTeamData **) (root_scratch + m_pool_members);

      // team size == 1, league size == pool_size

      for ( int rank = 0 ; rank < size ; ++rank ) {
        HostThreadTeamData * const mem = members[ rank ] ;
        mem->m_pool_scratch = root_scratch ;
        mem->m_team_scratch = mem->m_scratch ;
        mem->m_pool_rank    = rank ;
        mem->m_pool_size    = size ;
        mem->m_team_base    = rank ;
        mem->m_team_rank    = 0 ;
        mem->m_team_size    = 1 ;
        mem->m_team_alloc   = 1 ;
        mem->m_league_rank  = rank ;
        mem->m_league_size  = size ;
        mem->m_team_rendezvous_step = 0 ;
        pool[ rank ] = mem ;
      }
    }

    Kokkos::memory_fence();
  }
  else {
    Kokkos::Impl::throw_runtime_exception("Kokkos::Impl::HostThreadTeamData::organize_pool ERROR pool already exists");
  }
}

void HostThreadTeamData::disband_pool()
{
   m_work_range.first  = -1 ;
   m_work_range.second = -1 ;
   m_pool_scratch = 0 ;
   m_team_scratch = 0 ;
   m_pool_rank    = 0 ;
   m_pool_size    = 1 ;
   m_team_base    = 0 ;
   m_team_rank    = 0 ;
   m_team_size    = 1 ;
   m_team_alloc   = 1 ;
   m_league_rank  = 0 ;
   m_league_size  = 1 ;
   m_team_rendezvous_step = 0 ;
}

int HostThreadTeamData::organize_team( const int team_size )
{
  // Pool is initialized
  const bool ok_pool = 0 != m_pool_scratch ;

  // Team is not set
  const bool ok_team =
    m_team_scratch == m_scratch &&
    m_team_base    == m_pool_rank &&
    m_team_rank    == 0 &&
    m_team_size    == 1 &&
    m_team_alloc   == 1 &&
    m_league_rank  == m_pool_rank &&
    m_league_size  == m_pool_size ;

  if ( ok_pool && ok_team ) {

    if ( team_size <= 0 ) return 0 ; // No teams to organize

    if ( team_size == 1 ) return 1 ; // Already organized in teams of one

    HostThreadTeamData * const * const pool =
      (HostThreadTeamData **) (m_pool_scratch + m_pool_members);

    // "league_size" in this context is the number of concurrent teams
    // that the pool can accommodate.  Excess threads are idle.
    const int league_size     = m_pool_size / team_size ;
    const int team_alloc_size = m_pool_size / league_size ;
    const int team_alloc_rank = m_pool_rank % team_alloc_size ;
    const int league_rank     = m_pool_rank / team_alloc_size ;
    const int team_base_rank  = league_rank * team_alloc_size ;

    m_team_scratch = pool[ team_base_rank ]->m_scratch ;
    m_team_base    = team_base_rank ;
    // This needs to check overflow, if m_pool_size % team_alloc_size !=0
    // there are two corner cases:
    // (i) if team_alloc_size == team_size there might be a non-full
    //     zombi team around (for example m_pool_size = 5 and team_size = 2
    // (ii) if team_alloc > team_size then the last team might have less
    //      threads than the others
    m_team_rank    = ( team_base_rank + team_size <= m_pool_size ) &&
                     ( team_alloc_rank < team_size ) ?
                     team_alloc_rank : -1;
    m_team_size    = team_size ;
    m_team_alloc   = team_alloc_size ;
    m_league_rank  = league_rank ;
    m_league_size  = league_size ;
    m_team_rendezvous_step = 0 ;

    if ( team_base_rank == m_pool_rank ) {
      // Initialize team's rendezvous memory
      for ( int i = m_team_rendezvous ; i < m_pool_reduce ; ++i ) {
        m_scratch[i] = 0 ;
      }
      // Make sure team's rendezvous memory initialized
      // is written before proceeding.
      Kokkos::memory_fence();
    }

    // Organizing threads into a team performs a barrier across the
    // entire pool to insure proper initialization of the team
    // rendezvous mechanism before a team rendezvous can be performed.

    if ( pool_rendezvous() ) {
      pool_rendezvous_release();
    }
  }
  else {
    Kokkos::Impl::throw_runtime_exception("Kokkos::Impl::HostThreadTeamData::organize_team ERROR");
  }

  return 0 <= m_team_rank ;
}

void HostThreadTeamData::disband_team()
{
  m_team_scratch = m_scratch ;
  m_team_base    = m_pool_rank ;
  m_team_rank    = 0 ;
  m_team_size    = 1 ;
  m_team_alloc   = 1 ;
  m_league_rank  = m_pool_rank ;
  m_league_size  = m_pool_size ;
  m_team_rendezvous_step = 0 ;
}

//----------------------------------------------------------------------------
/* pattern for rendezvous
 *
 *  if ( rendezvous() ) {
 *     ... all other threads are still in team_rendezvous() ...
 *     rendezvous_release();
 *     ... all other threads are released from team_rendezvous() ...
 *  }
 */

int HostThreadTeamData::rendezvous( int64_t * const buffer
                                  , int & rendezvous_step
                                  , int const size
                                  , int const rank ) noexcept
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

  const int step = ( rendezvous_step % size_val_cycle ) + 1 ;

  rendezvous_step = step ;

  // The leading int64_t[4] span is for thread 0 to write
  // and all other threads to read spin-wait.
  // sync_offset is the index into this array for this step.

  const int sync_offset = ( step & mask_mem_cycle ) + size_mem_cycle ;

  if ( rank ) {

    const int group_begin = rank << shift_byte ; // == rank * size_byte

    if ( group_begin < size ) {

      //  This thread waits for threads
      //   [ group_begin .. group_begin + 8 )
      //   [ rank*8      .. rank*8 + 8      )
      // to write to their designated bytes.

      const int end = group_begin + size_byte < size
                    ? size_byte : size - group_begin ;

      int64_t value = 0 ;

      for ( int i = 0 ; i < end ; ++i ) {
        ((int8_t*) & value )[i] = int8_t( step );
      }
      // Do not REMOVE this store fence!!!
      // Makes stuff hang on GCC with more than 8 threads
      store_fence();
      spinwait_until_equal( buffer[ (rank << shift_mem_cycle) + sync_offset ]
                          , value );
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

    spinwait_until_equal( buffer[ step & mask_mem_cycle ] , int64_t(step) );

  }
  else {
    // Thread 0 waits for threads [1..7]
    // to write to their designated bytes.

    const int end = size_byte < size ? 8 : size ;

    int64_t value = 0 ;
    for ( int i = 1 ; i < end ; ++i ) {
      ((int8_t *) & value)[i] = int8_t( step );
    }

    spinwait_until_equal( buffer[ sync_offset ], value );
  }

  return rank ? 0 : 1 ;
}

void HostThreadTeamData::
  rendezvous_release( int64_t * const buffer
                    , int const rendezvous_step ) noexcept
{
  enum : int { shift_mem_cycle = 2 };
  enum : int { size_mem_cycle  = ( 01 << shift_mem_cycle ) }; // == 4
  enum : int { mask_mem_cycle  = size_mem_cycle - 1 };

  // Requires:
  //   Called after team_rendezvous
  //   Called only by true == team_rendezvous(root)

  // Memory fence to be sure all previous writes are complete:
  Kokkos::memory_fence();

  ((volatile int64_t*) buffer)[ rendezvous_step & mask_mem_cycle ] =
     int64_t( rendezvous_step );

  // Memory fence to push the store out
  Kokkos::memory_fence();
}

//----------------------------------------------------------------------------

int HostThreadTeamData::get_work_stealing() noexcept
{
  pair_int_t w( -1 , -1 );

  if ( 1 == m_team_size || team_rendezvous() ) {

    // Attempt first from beginning of my work range
    for ( int attempt = m_work_range.first < m_work_range.second ; attempt ; ) {

      // Query and attempt to update m_work_range
      //   from: [ w.first     , w.second )
      //   to:   [ w.first + 1 , w.second ) = w_new
      //
      // If w is invalid then is just a query.

      const pair_int_t w_new( w.first + 1 , w.second );

      w = Kokkos::atomic_compare_exchange( & m_work_range, w, w_new );

      if ( w.first < w.second ) {
        // m_work_range is viable

        // If steal is successful then don't repeat attempt to steal
        attempt = ! ( w_new.first  == w.first + 1 &&
                      w_new.second == w.second );
      }
      else {
        // m_work_range is not viable
        w.first  = -1 ;
        w.second = -1 ;

        attempt = 0 ;
      }
    }

    if ( w.first == -1 && m_steal_rank != m_pool_rank ) {

      HostThreadTeamData * const * const pool =
        (HostThreadTeamData**)( m_pool_scratch + m_pool_members );

      // Attempt from begining failed, try to steal from end of neighbor

      pair_int_t volatile * steal_range =
        & ( pool[ m_steal_rank ]->m_work_range );

      for ( int attempt = true ; attempt ; ) {

        // Query and attempt to update steal_work_range
        //   from: [ w.first , w.second )
        //   to:   [ w.first , w.second - 1 ) = w_new
        //
        // If w is invalid then is just a query.

        const pair_int_t w_new( w.first , w.second - 1 );

        w = Kokkos::atomic_compare_exchange( steal_range, w, w_new );

        if ( w.first < w.second ) {
          // steal_work_range is viable

          // If steal is successful then don't repeat attempt to steal
          attempt = ! ( w_new.first  == w.first &&
                        w_new.second == w.second - 1 );
        }
        else {
          // steal_work_range is not viable, move to next member
          w.first  = -1 ;
          w.second = -1 ;

          // We need to figure out whether the next team is active
          // m_steal_rank + m_team_alloc could be the next base_rank to steal from
          // but only if there are another m_team_size threads available so that that
          // base rank has a full team.
          m_steal_rank = m_steal_rank + m_team_alloc + m_team_size <= m_pool_size ?
                         m_steal_rank + m_team_alloc : 0;

          steal_range = & ( pool[ m_steal_rank ]->m_work_range );

          // If tried all other members then don't repeat attempt to steal
          attempt = m_steal_rank != m_pool_rank ;
        }
      }

      if ( w.first != -1 ) w.first = w.second - 1 ;
    }

    if ( 1 < m_team_size ) {
      // Must share the work index
      *((int volatile *) team_reduce()) = w.first ;

      team_rendezvous_release();
    }
  }
  else if ( 1 < m_team_size ) {
    w.first = *((int volatile *) team_reduce());
  }

  // May exit because successfully stole work and w is good.
  // May exit because no work left to steal and w = (-1,-1).

#if 0
fprintf(stdout,"HostThreadTeamData::get_work_stealing() pool(%d of %d) %d\n"
       , m_pool_rank , m_pool_size , w.first );
fflush(stdout);
#endif

  return w.first ;
}

} // namespace Impl
} // namespace Kokkos

