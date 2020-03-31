/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#include <limits>
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Spinwait.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void HostThreadTeamData::organize_pool(HostThreadTeamData *members[],
                                       const int size) {
  bool ok = true;

  memory_fence();

  // Verify not already a member of a pool:
  for (int rank = 0; rank < size && ok; ++rank) {
    ok = (nullptr != members[rank]) && (0 == members[rank]->m_pool_scratch);
  }

  if (ok) {
    int64_t *const root_scratch = members[0]->m_scratch;

    for (int i = m_pool_rendezvous; i < m_pool_reduce; ++i) {
      root_scratch[i] = 0;
    }

    {
      HostThreadTeamData **const pool =
          (HostThreadTeamData **)(root_scratch + m_pool_members);

      // team size == 1, league size == pool_size

      for (int rank = 0; rank < size; ++rank) {
        HostThreadTeamData *const mem = members[rank];
        mem->m_pool_scratch           = root_scratch;
        mem->m_team_scratch           = mem->m_scratch;
        mem->m_pool_rank              = rank;
        mem->m_pool_size              = size;
        mem->m_team_base              = rank;
        mem->m_team_rank              = 0;
        mem->m_team_size              = 1;
        mem->m_team_alloc             = 1;
        mem->m_league_rank            = rank;
        mem->m_league_size            = size;
        mem->m_team_rendezvous_step   = 0;
        pool[rank]                    = mem;
      }
    }

    Kokkos::memory_fence();
  } else {
    Kokkos::Impl::throw_runtime_exception(
        "Kokkos::Impl::HostThreadTeamData::organize_pool ERROR pool already "
        "exists");
  }
}

void HostThreadTeamData::disband_pool() {
  m_work_range.first     = -1;
  m_work_range.second    = -1;
  m_pool_scratch         = 0;
  m_team_scratch         = 0;
  m_pool_rank            = 0;
  m_pool_size            = 1;
  m_team_base            = 0;
  m_team_rank            = 0;
  m_team_size            = 1;
  m_team_alloc           = 1;
  m_league_rank          = 0;
  m_league_size          = 1;
  m_team_rendezvous_step = 0;
}

int HostThreadTeamData::organize_team(const int team_size) {
  // Pool is initialized
  const bool ok_pool = 0 != m_pool_scratch;

  // Team is not set
  const bool ok_team =
      m_team_scratch == m_scratch && m_team_base == m_pool_rank &&
      m_team_rank == 0 && m_team_size == 1 && m_team_alloc == 1 &&
      m_league_rank == m_pool_rank && m_league_size == m_pool_size;

  if (ok_pool && ok_team) {
    if (team_size <= 0) return 0;  // No teams to organize

    if (team_size == 1) return 1;  // Already organized in teams of one

    HostThreadTeamData *const *const pool =
        (HostThreadTeamData **)(m_pool_scratch + m_pool_members);

    // "league_size" in this context is the number of concurrent teams
    // that the pool can accommodate.  Excess threads are idle.
    const int league_size     = m_pool_size / team_size;
    const int team_alloc_size = m_pool_size / league_size;
    const int team_alloc_rank = m_pool_rank % team_alloc_size;
    const int league_rank     = m_pool_rank / team_alloc_size;
    const int team_base_rank  = league_rank * team_alloc_size;

    m_team_scratch = pool[team_base_rank]->m_scratch;
    m_team_base    = team_base_rank;
    // This needs to check overflow, if m_pool_size % team_alloc_size !=0
    // there are two corner cases:
    // (i) if team_alloc_size == team_size there might be a non-full
    //     zombi team around (for example m_pool_size = 5 and team_size = 2
    // (ii) if team_alloc > team_size then the last team might have less
    //      threads than the others
    m_team_rank = (team_base_rank + team_size <= m_pool_size) &&
                          (team_alloc_rank < team_size)
                      ? team_alloc_rank
                      : -1;
    m_team_size            = team_size;
    m_team_alloc           = team_alloc_size;
    m_league_rank          = league_rank;
    m_league_size          = league_size;
    m_team_rendezvous_step = 0;

    if (team_base_rank == m_pool_rank) {
      // Initialize team's rendezvous memory
      for (int i = m_team_rendezvous; i < m_pool_reduce; ++i) {
        m_scratch[i] = 0;
      }
      // Make sure team's rendezvous memory initialized
      // is written before proceeding.
      Kokkos::memory_fence();
    }

    // Organizing threads into a team performs a barrier across the
    // entire pool to insure proper initialization of the team
    // rendezvous mechanism before a team rendezvous can be performed.

    if (pool_rendezvous()) {
      pool_rendezvous_release();
    }
  } else {
    Kokkos::Impl::throw_runtime_exception(
        "Kokkos::Impl::HostThreadTeamData::organize_team ERROR");
  }

  return 0 <= m_team_rank;
}

void HostThreadTeamData::disband_team() {
  m_team_scratch         = m_scratch;
  m_team_base            = m_pool_rank;
  m_team_rank            = 0;
  m_team_size            = 1;
  m_team_alloc           = 1;
  m_league_rank          = m_pool_rank;
  m_league_size          = m_pool_size;
  m_team_rendezvous_step = 0;
}

//----------------------------------------------------------------------------

int HostThreadTeamData::get_work_stealing() noexcept {
  pair_int_t w(-1, -1);

  // TODO DJS 3-17-2018:
  // Discover why the work stealing algorithm only works when called
  // by the master thread of the team.  If we can refactor this section to
  // remove that requirement we should be able to remove the split_master_wait
  // behavior in the team and pool rendezvous algorithms
  if (1 == m_team_size || team_rendezvous()) {
    // Attempt first from beginning of my work range
    for (int attempt = m_work_range.first < m_work_range.second; attempt;) {
      // Query and attempt to update m_work_range
      //   from: [ w.first     , w.second )
      //   to:   [ w.first + 1 , w.second ) = w_new
      //
      // If w is invalid then is just a query.

      const pair_int_t w_new(w.first + 1, w.second);

      w = Kokkos::atomic_compare_exchange(&m_work_range, w, w_new);

      if (w.first < w.second) {
        // m_work_range is viable

        // If steal is successful then don't repeat attempt to steal
        attempt = !(w_new.first == w.first + 1 && w_new.second == w.second);
      } else {
        // m_work_range is not viable
        w.first  = -1;
        w.second = -1;

        attempt = 0;
      }
    }

    if (w.first == -1 && m_steal_rank != m_pool_rank) {
      HostThreadTeamData *const *const pool =
          (HostThreadTeamData **)(m_pool_scratch + m_pool_members);

      // Attempt from beginning failed, try to steal from end of neighbor

      pair_int_t volatile *steal_range = &(pool[m_steal_rank]->m_work_range);

      for (int attempt = true; attempt;) {
        // Query and attempt to update steal_work_range
        //   from: [ w.first , w.second )
        //   to:   [ w.first , w.second - 1 ) = w_new
        //
        // If w is invalid then is just a query.

        const pair_int_t w_new(w.first, w.second - 1);

        w = Kokkos::atomic_compare_exchange(steal_range, w, w_new);

        if (w.first < w.second) {
          // steal_work_range is viable

          // If steal is successful then don't repeat attempt to steal
          attempt = !(w_new.first == w.first && w_new.second == w.second - 1);
        } else {
          // steal_work_range is not viable, move to next member
          w.first  = -1;
          w.second = -1;

          // We need to figure out whether the next team is active
          // m_steal_rank + m_team_alloc could be the next base_rank to steal
          // from but only if there are another m_team_size threads available so
          // that that base rank has a full team.
          m_steal_rank =
              m_steal_rank + m_team_alloc + m_team_size <= m_pool_size
                  ? m_steal_rank + m_team_alloc
                  : 0;

          steal_range = &(pool[m_steal_rank]->m_work_range);

          // If tried all other members then don't repeat attempt to steal
          attempt = m_steal_rank != m_pool_rank;
        }
      }

      if (w.first != -1) w.first = w.second - 1;
    }

    if (1 < m_team_size) {
      // Must share the work index
      *((int volatile *)team_reduce()) = w.first;

      team_rendezvous_release();
    }
  } else if (1 < m_team_size) {
    w.first = *((int volatile *)team_reduce());
  }

  // May exit because successfully stole work and w is good.
  // May exit because no work left to steal and w = (-1,-1).

#if 0
fprintf(stdout,"HostThreadTeamData::get_work_stealing() pool(%d of %d) %d\n"
       , m_pool_rank , m_pool_size , w.first );
fflush(stdout);
#endif

  return w.first;
}

}  // namespace Impl
}  // namespace Kokkos
