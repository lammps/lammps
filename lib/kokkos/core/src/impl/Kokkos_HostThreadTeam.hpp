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

#ifndef KOKKOS_IMPL_HOSTTHREADTEAM_HPP
#define KOKKOS_IMPL_HOSTTHREADTEAM_HPP

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_Atomic.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <impl/Kokkos_FunctorAdapter.hpp>
#include <impl/Kokkos_FunctorAnalysis.hpp>
#include <impl/Kokkos_HostBarrier.hpp>

#include <limits>     // std::numeric_limits
#include <algorithm>  // std::max

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class HostExecSpace>
class HostThreadTeamMember;

class HostThreadTeamData {
 public:
  template <class>
  friend class HostThreadTeamMember;

  // Assume upper bounds on number of threads:
  //   pool size       <= 1024 threads
  //   team size       <= 64 threads

  enum : int { max_pool_members = 1024 };
  enum : int { max_team_members = 64 };
  enum : int { max_pool_rendezvous = HostBarrier::required_buffer_size };
  enum : int { max_team_rendezvous = HostBarrier::required_buffer_size };

 private:
  // per-thread scratch memory buffer chunks:
  //
  //   [ pool_members ]     = [ m_pool_members    .. m_pool_rendezvous )
  //   [ pool_rendezvous ]  = [ m_pool_rendezvous .. m_team_rendezvous )
  //   [ team_rendezvous ]  = [ m_team_rendezvous .. m_pool_reduce )
  //   [ pool_reduce ]      = [ m_pool_reduce     .. m_team_reduce )
  //   [ team_reduce ]      = [ m_team_reduce     .. m_team_shared )
  //   [ team_shared ]      = [ m_team_shared     .. m_thread_local )
  //   [ thread_local ]     = [ m_thread_local    .. m_scratch_size )

  enum : int { m_pool_members = 0 };
  enum : int {
    m_pool_rendezvous =
        static_cast<int>(m_pool_members) + static_cast<int>(max_pool_members)
  };
  enum : int {
    m_team_rendezvous = static_cast<int>(m_pool_rendezvous) +
                        static_cast<int>(max_pool_rendezvous)
  };
  enum : int {
    m_pool_reduce = static_cast<int>(m_team_rendezvous) +
                    static_cast<int>(max_team_rendezvous)
  };

  using pair_int_t = Kokkos::pair<int64_t, int64_t>;

  pair_int_t m_work_range;
  int64_t m_work_end;
  int64_t* m_scratch;       // per-thread buffer
  int64_t* m_pool_scratch;  // == pool[0]->m_scratch
  int64_t* m_team_scratch;  // == pool[ 0 + m_team_base ]->m_scratch
  int m_pool_rank;
  int m_pool_size;
  int m_team_reduce;
  int m_team_shared;
  int m_thread_local;
  int m_scratch_size;
  int m_team_base;
  int m_team_rank;
  int m_team_size;
  int m_team_alloc;
  int m_league_rank;
  int m_league_size;
  int m_work_chunk;
  int m_steal_rank;  // work stealing rank
  int mutable m_pool_rendezvous_step;
  int mutable m_team_rendezvous_step;

  HostThreadTeamData* team_member(int r) const noexcept {
    return (reinterpret_cast<HostThreadTeamData**>(
        m_pool_scratch + m_pool_members))[m_team_base + r];
  }

 public:
  inline bool team_rendezvous() const noexcept {
    int* ptr = reinterpret_cast<int*>(m_team_scratch + m_team_rendezvous);
    HostBarrier::split_arrive(ptr, m_team_size, m_team_rendezvous_step);
    if (m_team_rank != 0) {
      HostBarrier::wait(ptr, m_team_size, m_team_rendezvous_step);
    } else {
      HostBarrier::split_master_wait(ptr, m_team_size, m_team_rendezvous_step);
    }

    return m_team_rank == 0;
  }

  inline bool team_rendezvous(const int source_team_rank) const noexcept {
    int* ptr = reinterpret_cast<int*>(m_team_scratch + m_team_rendezvous);
    HostBarrier::split_arrive(ptr, m_team_size, m_team_rendezvous_step);
    if (m_team_rank != source_team_rank) {
      HostBarrier::wait(ptr, m_team_size, m_team_rendezvous_step);
    } else {
      HostBarrier::split_master_wait(ptr, m_team_size, m_team_rendezvous_step);
    }

    return (m_team_rank == source_team_rank);
  }

  inline void team_rendezvous_release() const noexcept {
    HostBarrier::split_release(
        reinterpret_cast<int*>(m_team_scratch + m_team_rendezvous), m_team_size,
        m_team_rendezvous_step);
  }

  inline int pool_rendezvous() const noexcept {
    int* ptr = reinterpret_cast<int*>(m_pool_scratch + m_pool_rendezvous);
    HostBarrier::split_arrive(ptr, m_pool_size, m_pool_rendezvous_step);
    if (m_pool_rank != 0) {
      HostBarrier::wait(ptr, m_pool_size, m_pool_rendezvous_step);
    } else {
      HostBarrier::split_master_wait(ptr, m_pool_size, m_pool_rendezvous_step);
    }

    return m_pool_rank == 0;
  }

  inline void pool_rendezvous_release() const noexcept {
    HostBarrier::split_release(
        reinterpret_cast<int*>(m_pool_scratch + m_pool_rendezvous), m_pool_size,
        m_pool_rendezvous_step);
  }

  //----------------------------------------

  constexpr HostThreadTeamData() noexcept
      : m_work_range(-1, -1),
        m_work_end(0),
        m_scratch(nullptr),
        m_pool_scratch(nullptr),
        m_team_scratch(nullptr),
        m_pool_rank(0),
        m_pool_size(1),
        m_team_reduce(0),
        m_team_shared(0),
        m_thread_local(0),
        m_scratch_size(0),
        m_team_base(0),
        m_team_rank(0),
        m_team_size(1),
        m_team_alloc(1),
        m_league_rank(0),
        m_league_size(1),
        m_work_chunk(0),
        m_steal_rank(0),
        m_pool_rendezvous_step(0),
        m_team_rendezvous_step(0) {}

  //----------------------------------------
  // Organize array of members into a pool.
  // The 0th member is the root of the pool.
  // Requires: members are not already in a pool.
  // Requires: called by one thread.
  // Pool members are ordered as "close" - sorted by NUMA and then CORE
  // Each thread is its own team with team_size == 1.
  static void organize_pool(HostThreadTeamData* members[], const int size);

  // Called by each thread within the pool
  void disband_pool();

  //----------------------------------------
  // Each thread within a pool organizes itself into a team.
  // Must be called by all threads of the pool.
  // Organizing threads into a team performs a barrier across the
  // entire pool to insure proper initialization of the team
  // rendezvous mechanism before a team rendezvous can be performed.
  //
  // Return true  if a valid member of a team.
  // Return false if not a member and thread should be idled.
  int organize_team(const int team_size);

  // Each thread within a pool disbands itself from current team.
  // Each thread becomes its own team with team_size == 1.
  // Must be called by all threads of the pool.
  void disband_team();

  //----------------------------------------

  constexpr int pool_rank() const { return m_pool_rank; }
  constexpr int pool_size() const { return m_pool_size; }

  HostThreadTeamData* pool_member(int r) const noexcept {
    return (reinterpret_cast<HostThreadTeamData**>(m_pool_scratch +
                                                   m_pool_members))[r];
  }

  //----------------------------------------

 private:
  enum : int { mask_to_16 = 0x0f };  // align to 16 bytes
  enum : int { shift_to_8 = 3 };     // size to 8 bytes

 public:
  static constexpr int align_to_int64(int n) {
    return ((n + mask_to_16) & ~mask_to_16) >> shift_to_8;
  }

  constexpr int pool_reduce_bytes() const {
    return m_scratch_size ? sizeof(int64_t) * (m_team_reduce - m_pool_reduce)
                          : 0;
  }

  constexpr int team_reduce_bytes() const {
    return sizeof(int64_t) * (m_team_shared - m_team_reduce);
  }

  constexpr int team_shared_bytes() const {
    return sizeof(int64_t) * (m_thread_local - m_team_shared);
  }

  constexpr int thread_local_bytes() const {
    return sizeof(int64_t) * (m_scratch_size - m_thread_local);
  }

  constexpr int scratch_bytes() const {
    return sizeof(int64_t) * m_scratch_size;
  }

  // Memory chunks:

  int64_t* scratch_buffer() const noexcept { return m_scratch; }

  int64_t* pool_reduce() const noexcept {
    return m_pool_scratch + m_pool_reduce;
  }

  int64_t* pool_reduce_local() const noexcept {
    return m_scratch + m_pool_reduce;
  }

  int64_t* team_reduce() const noexcept {
    return m_team_scratch + m_team_reduce;
  }

  int64_t* team_reduce_local() const noexcept {
    return m_scratch + m_team_reduce;
  }

  int64_t* team_shared() const noexcept {
    return m_team_scratch + m_team_shared;
  }

  int64_t* local_scratch() const noexcept { return m_scratch + m_thread_local; }

  // Given:
  //   pool_reduce_size  = number bytes for pool reduce
  //   team_reduce_size  = number bytes for team reduce
  //   team_shared_size  = number bytes for team shared memory
  //   thread_local_size = number bytes for thread local memory
  // Return:
  //   total number of bytes that must be allocated
  static size_t scratch_size(int pool_reduce_size, int team_reduce_size,
                             int team_shared_size, int thread_local_size) {
    pool_reduce_size  = align_to_int64(pool_reduce_size);
    team_reduce_size  = align_to_int64(team_reduce_size);
    team_shared_size  = align_to_int64(team_shared_size);
    thread_local_size = align_to_int64(thread_local_size);

    const size_t total_bytes =
        (m_pool_reduce + pool_reduce_size + team_reduce_size +
         team_shared_size + thread_local_size) *
        sizeof(int64_t);

    return total_bytes;
  }

  // Given:
  //   alloc_ptr         = pointer to allocated memory
  //   alloc_size        = number bytes of allocated memory
  //   pool_reduce_size  = number bytes for pool reduce/scan operations
  //   team_reduce_size  = number bytes for team reduce/scan operations
  //   team_shared_size  = number bytes for team-shared memory
  //   thread_local_size = number bytes for thread-local memory
  // Return:
  //   total number of bytes that must be allocated
  void scratch_assign(void* const alloc_ptr, size_t const alloc_size,
                      int pool_reduce_size, int team_reduce_size,
                      int team_shared_size, int /* thread_local_size */) {
    pool_reduce_size = align_to_int64(pool_reduce_size);
    team_reduce_size = align_to_int64(team_reduce_size);
    team_shared_size = align_to_int64(team_shared_size);
    // thread_local_size = align_to_int64( thread_local_size );

    m_scratch      = static_cast<int64_t*>(alloc_ptr);
    m_team_reduce  = m_pool_reduce + pool_reduce_size;
    m_team_shared  = m_team_reduce + team_reduce_size;
    m_thread_local = m_team_shared + team_shared_size;
    m_scratch_size = align_to_int64(alloc_size);
  }

  //----------------------------------------
  // Get a work index within the range.
  // First try to steal from beginning of own teams's partition.
  // If that fails then try to steal from end of another teams' partition.
  int get_work_stealing() noexcept;

  //----------------------------------------
  // Set the initial work partitioning of [ 0 .. length ) among the teams
  // with granularity of chunk

  void set_work_partition(int64_t const length, int const chunk) noexcept {
    // Minimum chunk size to insure that
    //   m_work_end < std::numeric_limits<int>::max() * m_work_chunk

    int const chunk_min = (length + std::numeric_limits<int>::max()) /
                          std::numeric_limits<int>::max();

    m_work_end   = length;
    m_work_chunk = std::max(chunk, chunk_min);

    // Number of work chunks and partitioning of that number:
    int const num  = (m_work_end + m_work_chunk - 1) / m_work_chunk;
    int const part = (num + m_league_size - 1) / m_league_size;

    m_work_range.first  = part * m_league_rank;
    m_work_range.second = m_work_range.first + part;

    // Steal from next team, round robin
    // The next team is offset by m_team_alloc if it fits in the pool.

    m_steal_rank = m_team_base + m_team_alloc + m_team_size <= m_pool_size
                       ? m_team_base + m_team_alloc
                       : 0;
  }

  std::pair<int64_t, int64_t> get_work_partition() noexcept {
    int64_t first  = m_work_range.first;
    int64_t second = m_work_range.second;
    first *= m_work_chunk;
    second *= m_work_chunk;
    return std::pair<int64_t, int64_t>(
        first, second < m_work_end ? second : m_work_end);
  }

  std::pair<int64_t, int64_t> get_work_stealing_chunk() noexcept {
    std::pair<int64_t, int64_t> x(-1, -1);

    const int i = get_work_stealing();

    if (0 <= i) {
      x.first  = m_work_chunk * i;
      x.second = x.first + m_work_chunk < m_work_end ? x.first + m_work_chunk
                                                     : m_work_end;
    }

    return x;
  }
};

//----------------------------------------------------------------------------

template <class HostExecSpace>
class HostThreadTeamMember {
 public:
  using scratch_memory_space    = typename HostExecSpace::scratch_memory_space;
  using execution_space         = HostExecSpace;
  using thread_team_member      = HostThreadTeamMember;
  using host_thread_team_member = HostThreadTeamMember;

 private:
  scratch_memory_space m_scratch;
  HostThreadTeamData& m_data;
  int const m_league_rank;
  int const m_league_size;

 public:
  constexpr HostThreadTeamMember(HostThreadTeamData& arg_data) noexcept
      : m_scratch(arg_data.team_shared(), arg_data.team_shared_bytes()),
        m_data(arg_data),
        m_league_rank(arg_data.m_league_rank),
        m_league_size(arg_data.m_league_size) {}

  constexpr HostThreadTeamMember(HostThreadTeamData& arg_data,
                                 int const arg_league_rank,
                                 int const arg_league_size) noexcept
      : m_scratch(arg_data.team_shared(), arg_data.team_shared_bytes(),
                  arg_data.team_shared(), arg_data.team_shared_bytes()),
        m_data(arg_data),
        m_league_rank(arg_league_rank),
        m_league_size(arg_league_size) {}

  ~HostThreadTeamMember()                           = default;
  HostThreadTeamMember()                            = delete;
  HostThreadTeamMember(HostThreadTeamMember&&)      = default;
  HostThreadTeamMember(HostThreadTeamMember const&) = default;
  HostThreadTeamMember& operator=(HostThreadTeamMember&&) = default;
  HostThreadTeamMember& operator=(HostThreadTeamMember const&) = default;

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  int team_rank() const noexcept { return m_data.m_team_rank; }

  KOKKOS_INLINE_FUNCTION
  int team_size() const noexcept { return m_data.m_team_size; }

  KOKKOS_INLINE_FUNCTION
  int league_rank() const noexcept { return m_league_rank; }

  KOKKOS_INLINE_FUNCTION
  int league_size() const noexcept { return m_league_size; }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  const scratch_memory_space& team_shmem() const {
    return m_scratch.set_team_thread_mode(0, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const scratch_memory_space& team_scratch(int) const {
    return m_scratch.set_team_thread_mode(0, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const scratch_memory_space& thread_scratch(int) const {
    return m_scratch.set_team_thread_mode(0, m_data.m_team_size,
                                          m_data.m_team_rank);
  }

  //--------------------------------------------------------------------------
  // Team collectives
  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION void team_barrier() const noexcept {
    KOKKOS_IF_ON_HOST(
        (if (m_data.team_rendezvous()) { m_data.team_rendezvous_release(); }))
  }

  //--------------------------------------------------------------------------

  template <typename T>
  KOKKOS_INLINE_FUNCTION void team_broadcast(T& value,
                                             const int source_team_rank) const
      noexcept {
    KOKKOS_IF_ON_HOST((if (1 < m_data.m_team_size) {
      T volatile* const shared_value = (T*)m_data.team_reduce();

      // Don't overwrite shared memory until all threads arrive

      if (m_data.team_rendezvous(source_team_rank)) {
        // All threads have entered 'team_rendezvous'
        // only this thread returned from 'team_rendezvous'
        // with a return value of 'true'

        *shared_value = value;

        m_data.team_rendezvous_release();
        // This thread released all other threads from 'team_rendezvous'
        // with a return value of 'false'
      } else {
        value = *shared_value;
      }
    }))

    KOKKOS_IF_ON_DEVICE(((void)value; (void)source_team_rank; Kokkos::abort(
                             "HostThreadTeamMember team_broadcast\n");))
  }

  //--------------------------------------------------------------------------

  template <class Closure, typename T>
  KOKKOS_INLINE_FUNCTION void team_broadcast(Closure const& f, T& value,
                                             const int source_team_rank) const
      noexcept {
    KOKKOS_IF_ON_HOST((
        T volatile* const shared_value = (T*)m_data.team_reduce();

        // Don't overwrite shared memory until all threads arrive

        if (m_data.team_rendezvous(source_team_rank)) {
          // All threads have entered 'team_rendezvous'
          // only this thread returned from 'team_rendezvous'
          // with a return value of 'true'

          f(value);

          if (1 < m_data.m_team_size) {
            *shared_value = value;
          }

          m_data.team_rendezvous_release();
          // This thread released all other threads from 'team_rendezvous'
          // with a return value of 'false'
        } else { value = *shared_value; }))

    KOKKOS_IF_ON_DEVICE(
        ((void)f; (void)value; (void)source_team_rank;
         Kokkos::abort("HostThreadTeamMember team_broadcast\n");))
  }

  //--------------------------------------------------------------------------
  // team_reduce( Sum(result) );
  // team_reduce( Min(result) );
  // team_reduce( Max(result) );

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      team_reduce(ReducerType const& reducer) const noexcept {
    team_reduce(reducer, reducer.reference());
  }

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      team_reduce(ReducerType const& reducer,
                  typename ReducerType::value_type contribution) const
      noexcept {
    KOKKOS_IF_ON_HOST((
        if (1 < m_data.m_team_size) {
          using value_type = typename ReducerType::value_type;

          if (0 != m_data.m_team_rank) {
            // Non-root copies to their local buffer:
            /*reducer.copy( (value_type*) m_data.team_reduce_local()
                        , reducer.data() );*/
            *((value_type*)m_data.team_reduce_local()) = contribution;
          }

          // Root does not overwrite shared memory until all threads arrive
          // and copy to their local buffer.

          if (m_data.team_rendezvous()) {
            // All threads have entered 'team_rendezvous'
            // only this thread returned from 'team_rendezvous'
            // with a return value of 'true'
            //
            // This thread sums contributed values
            for (int i = 1; i < m_data.m_team_size; ++i) {
              value_type* const src =
                  (value_type*)m_data.team_member(i)->team_reduce_local();

              reducer.join(contribution, *src);
            }

            // Copy result to root member's buffer:
            // reducer.copy( (value_type*) m_data.team_reduce() , reducer.data()
            // );
            *((value_type*)m_data.team_reduce()) = contribution;
            reducer.reference()                  = contribution;
            m_data.team_rendezvous_release();
            // This thread released all other threads from 'team_rendezvous'
            // with a return value of 'false'
          } else {
            // Copy from root member's buffer:
            reducer.reference() = *((value_type*)m_data.team_reduce());
          }
        } else { reducer.reference() = contribution; }))

    KOKKOS_IF_ON_DEVICE(((void)reducer; (void)contribution;
                         Kokkos::abort("HostThreadTeamMember team_reduce\n");))
  }

  //--------------------------------------------------------------------------

  template <typename T>
  KOKKOS_INLINE_FUNCTION T team_scan(T const& value,
                                     T* const global = nullptr) const noexcept {
    KOKKOS_IF_ON_HOST((
        if (0 != m_data.m_team_rank) {
          // Non-root copies to their local buffer:
          ((T*)m_data.team_reduce_local())[1] = value;
        }

        // Root does not overwrite shared memory until all threads arrive
        // and copy to their local buffer.

        if (m_data.team_rendezvous()) {
          // All threads have entered 'team_rendezvous'
          // only this thread returned from 'team_rendezvous'
          // with a return value of 'true'
          //
          // This thread scans contributed values

          {
            T* prev = (T*)m_data.team_reduce_local();

            prev[0] = 0;
            prev[1] = value;

            for (int i = 1; i < m_data.m_team_size; ++i) {
              T* const ptr = (T*)m_data.team_member(i)->team_reduce_local();

              ptr[0] = prev[0] + prev[1];

              prev = ptr;
            }
          }

          // If adding to global value then atomic_fetch_add to that value
          // and sum previous value to every entry of the scan.
          if (global) {
            T* prev = (T*)m_data.team_reduce_local();

            {
              T* ptr = (T*)m_data.team_member(m_data.m_team_size - 1)
                           ->team_reduce_local();
              prev[0] = Kokkos::atomic_fetch_add(global, ptr[0] + ptr[1]);
            }

            for (int i = 1; i < m_data.m_team_size; ++i) {
              T* ptr = (T*)m_data.team_member(i)->team_reduce_local();
              ptr[0] += prev[0];
            }
          }

          m_data.team_rendezvous_release();
        }

        return ((T*)m_data.team_reduce_local())[0];))

    KOKKOS_IF_ON_DEVICE(((void)value; (void)global;
                         Kokkos::abort("HostThreadTeamMember team_scan\n");
                         return T();))
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template <typename iType, typename Member>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<iType, Member>
TeamThreadRange(
    Member const& member, iType count,
    typename std::enable_if<
        Impl::is_thread_team_member<Member>::value>::type const** = nullptr) {
  return Impl::TeamThreadRangeBoundariesStruct<iType, Member>(member, 0, count);
}

template <typename iType1, typename iType2, typename Member>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, Member>
TeamThreadRange(
    Member const& member, iType1 begin, iType2 end,
    typename std::enable_if<
        Impl::is_thread_team_member<Member>::value>::type const** = nullptr) {
  return Impl::TeamThreadRangeBoundariesStruct<
      typename std::common_type<iType1, iType2>::type, Member>(member, begin,
                                                               end);
}

template <typename iType, typename Member>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<iType, Member>
TeamVectorRange(
    Member const& member, iType count,
    typename std::enable_if<
        Impl::is_thread_team_member<Member>::value>::type const** = nullptr) {
  return Impl::TeamThreadRangeBoundariesStruct<iType, Member>(member, 0, count);
}

template <typename iType1, typename iType2, typename Member>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, Member>
TeamVectorRange(
    Member const& member, iType1 begin, iType2 end,
    typename std::enable_if<
        Impl::is_thread_team_member<Member>::value>::type const** = nullptr) {
  return Impl::TeamThreadRangeBoundariesStruct<
      typename std::common_type<iType1, iType2>::type, Member>(member, begin,
                                                               end);
}

template <typename iType, typename Member>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<iType, Member>
ThreadVectorRange(
    Member const& member, iType count,
    typename std::enable_if<
        Impl::is_thread_team_member<Member>::value>::type const** = nullptr) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Member>(member, count);
}

template <typename iType1, typename iType2, typename Member>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, Member>
ThreadVectorRange(
    Member const& member, iType1 arg_begin, iType2 arg_end,
    typename std::enable_if<
        Impl::is_thread_team_member<Member>::value>::type const** = nullptr) {
  using iType = typename std::common_type<iType1, iType2>::type;
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Member>(
      member, iType(arg_begin), iType(arg_end));
}

//----------------------------------------------------------------------------
/** \brief  Inter-thread parallel_for.
 *
 * Executes lambda(iType i) for each i=[0..N)
 *
 * The range [0..N) is mapped to all threads of the the calling thread team.
 */
template <typename iType, class Closure, class Member>
KOKKOS_INLINE_FUNCTION void parallel_for(
    Impl::TeamThreadRangeBoundariesStruct<iType, Member> const& loop_boundaries,
    Closure const& closure,
    typename std::enable_if<Impl::is_host_thread_team_member<Member>::value>::
        type const** = nullptr) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    closure(i);
  }
}

template <typename iType, class Closure, class Member>
KOKKOS_INLINE_FUNCTION void parallel_for(
    Impl::ThreadVectorRangeBoundariesStruct<iType, Member> const&
        loop_boundaries,
    Closure const& closure,
    typename std::enable_if<Impl::is_host_thread_team_member<Member>::value>::
        type const** = nullptr) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    closure(i);
  }
}

//----------------------------------------------------------------------------

template <typename iType, class Closure, class Reducer, class Member>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    Kokkos::is_reducer<Reducer>::value &&
    Impl::is_host_thread_team_member<Member>::value>::type
parallel_reduce(
    Impl::TeamThreadRangeBoundariesStruct<iType, Member> const& loop_boundaries,
    Closure const& closure, Reducer const& reducer) {
  typename Reducer::value_type value;
  reducer.init(value);

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    closure(i, value);
  }

  loop_boundaries.thread.team_reduce(reducer, value);
}

template <typename iType, typename Closure, typename ValueType, typename Member>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    !Kokkos::is_reducer<ValueType>::value &&
    Impl::is_host_thread_team_member<Member>::value>::type
parallel_reduce(
    Impl::TeamThreadRangeBoundariesStruct<iType, Member> const& loop_boundaries,
    Closure const& closure, ValueType& result) {
  ValueType val;
  Sum<ValueType> reducer(val);
  reducer.init(val);

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    closure(i, reducer.reference());
  }

  loop_boundaries.thread.team_reduce(reducer);
  result = reducer.reference();
}

/*template< typename iType, class Space
         , class Closure, class Joiner , typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce
  (
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::HostThreadTeamMember<Space> >
             const & loop_boundaries
  , Closure  const & closure
  , Joiner   const & joiner
  , ValueType      & result
  )
{
  Impl::Reducer< ValueType , Joiner > reducer( joiner , & result );

  reducer.init( reducer.data() );

  for( iType i = loop_boundaries.start
     ; i <  loop_boundaries.end
     ; i += loop_boundaries.increment ) {
    closure( i , reducer.reference() );
  }

  loop_boundaries.thread.team_reduce( reducer );
}*/

//----------------------------------------------------------------------------
/** \brief  Inter-thread vector parallel_reduce.
 *
 *  Executes lambda(iType i, ValueType & val) for each i=[0..N)
 *
 *  The range [0..N) is mapped to all threads of the
 *  calling thread team and a summation of  val is
 *  performed and put into result.
 */
template <typename iType, class Lambda, typename ValueType, typename Member>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    !Kokkos::is_reducer<ValueType>::value &&
    Impl::is_host_thread_team_member<Member>::value>::type
parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType, Member>&
                    loop_boundaries,
                const Lambda& lambda, ValueType& result) {
  result = ValueType();
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }
}

template <typename iType, class Lambda, typename ReducerType, typename Member>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    Kokkos::is_reducer<ReducerType>::value &&
    Impl::is_host_thread_team_member<Member>::value>::type
parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<iType, Member>&
                    loop_boundaries,
                const Lambda& lambda, const ReducerType& reducer) {
  reducer.init(reducer.reference());
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, reducer.reference());
  }
}

//----------------------------------------------------------------------------

template <typename iType, class Closure, class Member>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    Impl::is_host_thread_team_member<Member>::value>::type
parallel_scan(
    Impl::TeamThreadRangeBoundariesStruct<iType, Member> const& loop_boundaries,
    Closure const& closure) {
  // Extract ValueType from the closure

  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, Closure>::value_type;

  value_type accum = 0;

  // Intra-member scan
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    closure(i, accum, false);
  }

  // 'accum' output is the exclusive prefix sum
  accum = loop_boundaries.thread.team_scan(accum);

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    closure(i, accum, true);
  }
}

template <typename iType, class ClosureType, class Member>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    Impl::is_host_thread_team_member<Member>::value>::type
parallel_scan(Impl::ThreadVectorRangeBoundariesStruct<iType, Member> const&
                  loop_boundaries,
              ClosureType const& closure) {
  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Impl::FunctorPatternInterface::SCAN, void, ClosureType>::value_type;

  value_type scan_val = value_type();

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    closure(i, scan_val, true);
  }
}

template <typename iType, class Lambda, typename ReducerType, typename Member>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    Kokkos::is_reducer<ReducerType>::value &&
    Impl::is_host_thread_team_member<Member>::value>::type
parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<iType, Member>&
                  loop_boundaries,
              const Lambda& lambda, const ReducerType& reducer) {
  typename ReducerType::value_type scan_val;
  reducer.init(scan_val);

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, scan_val, true);
  }
}

//----------------------------------------------------------------------------

template <class Member>
KOKKOS_INLINE_FUNCTION Impl::ThreadSingleStruct<Member> PerTeam(
    Member const& member,
    typename std::enable_if<
        Impl::is_thread_team_member<Member>::value>::type const** = nullptr) {
  return Impl::ThreadSingleStruct<Member>(member);
}

template <class Member>
KOKKOS_INLINE_FUNCTION Impl::VectorSingleStruct<Member> PerThread(
    Member const& member,
    typename std::enable_if<
        Impl::is_thread_team_member<Member>::value>::type const** = nullptr) {
  return Impl::VectorSingleStruct<Member>(member);
}

template <class Member, class FunctorType>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    Impl::is_host_thread_team_member<Member>::value>::type
single(const Impl::ThreadSingleStruct<Member>& single,
       const FunctorType& functor) {
  // 'single' does not perform a barrier.
  if (single.team_member.team_rank() == 0) functor();
}

template <class Member, class FunctorType, typename ValueType>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    Impl::is_host_thread_team_member<Member>::value>::type
single(const Impl::ThreadSingleStruct<Member>& single,
       const FunctorType& functor, ValueType& val) {
  single.team_member.team_broadcast(functor, val, 0);
}

template <class Member, class FunctorType>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    Impl::is_host_thread_team_member<Member>::value>::type
single(const Impl::VectorSingleStruct<Member>&, const FunctorType& functor) {
  functor();
}

template <class Member, class FunctorType, typename ValueType>
KOKKOS_INLINE_FUNCTION typename std::enable_if<
    Impl::is_host_thread_team_member<Member>::value>::type
single(const Impl::VectorSingleStruct<Member>&, const FunctorType& functor,
       ValueType& val) {
  functor(val);
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_HOSTTHREADTEAM_HPP */
