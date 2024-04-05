//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_THREADSTEAM_HPP
#define KOKKOS_THREADSTEAM_HPP

#include <Kokkos_Macros.hpp>

#include <cstdio>

#include <utility>
#include <impl/Kokkos_HostThreadTeam.hpp>

#include <Kokkos_Atomic.hpp>
#include <Threads/Kokkos_Threads_Spinwait.hpp>
#include <Threads/Kokkos_Threads_State.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template <class>
struct ThreadsExecAdapter;

//----------------------------------------------------------------------------

class ThreadsExecTeamMember {
 private:
  enum { TEAM_REDUCE_SIZE = 512 };

 public:
  using execution_space      = Kokkos::Threads;
  using scratch_memory_space = execution_space::scratch_memory_space;
  using team_handle          = ThreadsExecTeamMember;

 private:
  using space = execution_space::scratch_memory_space;
  ThreadsInternal* const m_instance;
  ThreadsInternal* const* m_team_base;  ///< Base for team fan-in
  space m_team_shared;
  size_t m_team_shared_size;
  int m_team_size;
  int m_team_rank;
  int m_team_rank_rev;
  int m_league_size;
  int m_league_end;
  int m_league_rank;

  int m_chunk_size;
  int m_league_chunk_end;

  int m_invalid_thread;
  int m_team_alloc;

  inline void set_team_shared() {
    new (&m_team_shared) space(
        static_cast<char*>((*m_team_base)->scratch_memory()) + TEAM_REDUCE_SIZE,
        m_team_shared_size);
  }

 public:
  // Fan-in and wait until the matching fan-out is called.
  // The root thread which does not wait will return true.
  // All other threads will return false during the fan-out.
  KOKKOS_INLINE_FUNCTION bool team_fan_in() const {
    int n, j;

    // Wait for fan-in threads
    for (n = 1;
         (!(m_team_rank_rev & n)) && ((j = m_team_rank_rev + n) < m_team_size);
         n <<= 1) {
      spinwait_while_equal(m_team_base[j]->state(), ThreadState::Active);
    }

    // If not root then wait for release
    if (m_team_rank_rev) {
      m_instance->state() = ThreadState::Rendezvous;
      spinwait_while_equal(m_instance->state(), ThreadState::Rendezvous);
    }

    return !m_team_rank_rev;
  }

  KOKKOS_INLINE_FUNCTION void team_fan_out() const {
    int n, j;
    for (n = 1;
         (!(m_team_rank_rev & n)) && ((j = m_team_rank_rev + n) < m_team_size);
         n <<= 1) {
      m_team_base[j]->state() = ThreadState::Active;
    }
  }

 public:
  KOKKOS_INLINE_FUNCTION static int team_reduce_size() {
    return TEAM_REDUCE_SIZE;
  }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_shmem() const {
    return m_team_shared.set_team_thread_mode(0, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_scratch(int) const {
    return m_team_shared.set_team_thread_mode(0, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& thread_scratch(int) const {
    return m_team_shared.set_team_thread_mode(0, team_size(), team_rank());
  }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return m_team_rank; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size; }

  KOKKOS_INLINE_FUNCTION void team_barrier() const {
    team_fan_in();
    team_fan_out();
  }

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(ValueType& value,
                                             const int& thread_id) const {
    KOKKOS_IF_ON_DEVICE(((void)value; (void)thread_id;))

    KOKKOS_IF_ON_HOST((
        // Make sure there is enough scratch space:
        using type = typename if_c<sizeof(ValueType) < TEAM_REDUCE_SIZE,
                                   ValueType, void>::type;

        if (m_team_base) {
          type* const local_value = ((type*)m_team_base[0]->scratch_memory());
          memory_fence();
          team_barrier();
          if (team_rank() == thread_id) *local_value = value;
          memory_fence();
          team_barrier();
          value = *local_value;
        }))
  }

  template <class Closure, class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(Closure const& f, ValueType& value,
                                             const int& thread_id) const {
    KOKKOS_IF_ON_DEVICE(((void)f; (void)value; (void)thread_id;))

    KOKKOS_IF_ON_HOST((
        // Make sure there is enough scratch space:
        using type = typename if_c<sizeof(ValueType) < TEAM_REDUCE_SIZE,
                                   ValueType, void>::type;
        f(value); if (m_team_base) {
          type* const local_value = ((type*)m_team_base[0]->scratch_memory());
          memory_fence();
          team_barrier();
          if (team_rank() == thread_id) *local_value = value;
          memory_fence();
          team_barrier();
          value = *local_value;
        }))
  }

  template <typename Type>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!Kokkos::is_reducer<Type>::value, Type>
      team_reduce(const Type& value) const {
    KOKKOS_IF_ON_DEVICE((return value;))

    KOKKOS_IF_ON_HOST((
        // Make sure there is enough scratch space:
        using type =
            typename if_c<sizeof(Type) < TEAM_REDUCE_SIZE, Type, void>::type;

        if (m_instance == nullptr) return value;

        if (team_rank() != team_size() - 1) *
            ((volatile type*)m_instance->scratch_memory()) = value;

        memory_fence();

        type& accum = *((type*)m_team_base[0]->scratch_memory());

        if (team_fan_in()) {
          accum = value;
          for (int i = 1; i < m_team_size; ++i) {
            accum += *((type*)m_team_base[i]->scratch_memory());
          }
          memory_fence();
        }

        team_fan_out();

        return accum;))
  }

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<is_reducer<ReducerType>::value>
  team_reduce(ReducerType const& reducer) const noexcept {
    team_reduce(reducer, reducer.reference());
  }

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
      team_reduce(const ReducerType& reducer,
                  const typename ReducerType::value_type contribution) const {
    KOKKOS_IF_ON_DEVICE(((void)reducer; (void)contribution;))

    KOKKOS_IF_ON_HOST((
        using value_type = typename ReducerType::value_type;
        // Make sure there is enough scratch space:
        using type = typename if_c<sizeof(value_type) < TEAM_REDUCE_SIZE,
                                   value_type, void>::type;

        if (m_instance == nullptr) return;

        type* const local_value = ((type*)m_instance->scratch_memory());

        // Set this thread's contribution
        if (team_rank() != team_size() - 1) { *local_value = contribution; }

        // Fence to make sure the base team member has access:
        memory_fence();

        if (team_fan_in()) {
          // The last thread to synchronize returns true, all other threads
          // wait for team_fan_out()
          type* const team_value = ((type*)m_team_base[0]->scratch_memory());

          *team_value = contribution;
          // Join to the team value:
          for (int i = 1; i < m_team_size; ++i) {
            reducer.join(*team_value,
                         *((type*)m_team_base[i]->scratch_memory()));
          }

          // Team base thread may "lap" member threads so copy out to their
          // local value.
          for (int i = 1; i < m_team_size; ++i) {
            *((type*)m_team_base[i]->scratch_memory()) = *team_value;
          }

          // Fence to make sure all team members have access
          memory_fence();
        }

        team_fan_out();

        // Value was changed by the team base
        reducer.reference() = *local_value;))
  }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template <typename ArgType>
  KOKKOS_INLINE_FUNCTION ArgType team_scan(const ArgType& value,
                                           ArgType* const global_accum) const {
    KOKKOS_IF_ON_DEVICE(((void)global_accum; return value;))

    KOKKOS_IF_ON_HOST((  // Make sure there is enough scratch space:
        using type = typename if_c<sizeof(ArgType) < TEAM_REDUCE_SIZE, ArgType,
                                   void>::type;

        if (m_instance == nullptr) return type(0);

        volatile type* const work_value = ((type*)m_instance->scratch_memory());

        *work_value = value;

        memory_fence();

        if (team_fan_in()) {
          // The last thread to synchronize returns true, all other threads wait
          // for team_fan_out() m_team_base[0]                 == highest
          // ranking team member m_team_base[ m_team_size - 1 ] == lowest
          // ranking team member
          //
          // 1) copy from lower to higher rank, initialize lowest rank to zero
          // 2) prefix sum from lowest to highest rank, skipping lowest rank

          type accum = 0;

          if (global_accum) {
            for (int i = m_team_size; i--;) {
              type& val = *((type*)m_team_base[i]->scratch_memory());
              accum += val;
            }
            accum = atomic_fetch_add(global_accum, accum);
          }

          for (int i = m_team_size; i--;) {
            type& val         = *((type*)m_team_base[i]->scratch_memory());
            const type offset = accum;
            accum += val;
            val = offset;
          }

          memory_fence();
        }

        team_fan_out();

        return *work_value;))
  }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template <typename ArgType>
  KOKKOS_INLINE_FUNCTION ArgType team_scan(const ArgType& value) const {
    return this->template team_scan<ArgType>(value, nullptr);
  }

  //----------------------------------------
  // Private for the driver

  template <class... Properties>
  ThreadsExecTeamMember(
      Impl::ThreadsInternal* instance,
      const TeamPolicyInternal<Kokkos::Threads, Properties...>& team,
      const size_t shared_size)
      : m_instance(instance),
        m_team_base(nullptr),
        m_team_shared(nullptr, 0),
        m_team_shared_size(shared_size),
        m_team_size(team.team_size()),
        m_team_rank(0),
        m_team_rank_rev(0),
        m_league_size(0),
        m_league_end(0),
        m_league_rank(0),
        m_chunk_size(team.chunk_size()),
        m_league_chunk_end(0),
        m_team_alloc(team.team_alloc()) {
    if (team.league_size()) {
      // Execution is using device-team interface:

      const int pool_rank_rev =
          m_instance->pool_size() - (m_instance->pool_rank() + 1);
      const int team_rank_rev = pool_rank_rev % team.team_alloc();
      const size_t pool_league_size =
          m_instance->pool_size() / team.team_alloc();
      const size_t pool_league_rank_rev = pool_rank_rev / team.team_alloc();
      if (pool_league_rank_rev >= pool_league_size) {
        m_invalid_thread = 1;
        return;
      }
      const size_t pool_league_rank =
          pool_league_size - (pool_league_rank_rev + 1);

      const int pool_num_teams = m_instance->pool_size() / team.team_alloc();
      const int chunk_size =
          team.chunk_size() > 0 ? team.chunk_size() : team.team_iter();
      const int chunks_per_team =
          (team.league_size() + chunk_size * pool_num_teams - 1) /
          (chunk_size * pool_num_teams);
      int league_iter_end = team.league_size() -
                            pool_league_rank_rev * chunks_per_team * chunk_size;
      int league_iter_begin = league_iter_end - chunks_per_team * chunk_size;
      if (league_iter_begin < 0) league_iter_begin = 0;
      if (league_iter_end > team.league_size())
        league_iter_end = team.league_size();

      if ((team.team_alloc() > size_t(m_team_size))
              ? (team_rank_rev >= m_team_size)
              : (m_instance->pool_size() - pool_num_teams * m_team_size >
                 m_instance->pool_rank()))
        m_invalid_thread = 1;
      else
        m_invalid_thread = 0;

      // May be using fewer threads per team than a multiple of threads per
      // core, some threads will idle.

      if (team_rank_rev < team.team_size() && !m_invalid_thread) {
        m_team_base =
            m_instance->pool_base() + team.team_alloc() * pool_league_rank_rev;
        m_team_size     = team.team_size();
        m_team_rank     = team.team_size() - (team_rank_rev + 1);
        m_team_rank_rev = team_rank_rev;
        m_league_size   = team.league_size();

        m_league_rank =
            (team.league_size() * pool_league_rank) / pool_league_size;
        m_league_end =
            (team.league_size() * (pool_league_rank + 1)) / pool_league_size;

        set_team_shared();
      }

      if ((m_team_rank_rev == 0) && (m_invalid_thread == 0)) {
        m_instance->set_work_range(m_league_rank, m_league_end, m_chunk_size);
        m_instance->reset_steal_target(m_team_size);
      }
      if (std::is_same<typename TeamPolicyInternal<
                           Kokkos::Threads, Properties...>::schedule_type::type,
                       Kokkos::Dynamic>::value) {
        m_instance->barrier();
      }
    } else {
      m_invalid_thread = 1;
    }
  }

  ThreadsExecTeamMember()
      : m_instance(nullptr),
        m_team_base(nullptr),
        m_team_shared(nullptr, 0),
        m_team_shared_size(0),
        m_team_size(1),
        m_team_rank(0),
        m_team_rank_rev(0),
        m_league_size(1),
        m_league_end(0),
        m_league_rank(0),
        m_chunk_size(0),
        m_league_chunk_end(0),
        m_invalid_thread(0),
        m_team_alloc(0) {}

  inline ThreadsInternal& threads_exec_team_base() const {
    return m_team_base ? **m_team_base : *m_instance;
  }

  bool valid_static() const { return m_league_rank < m_league_end; }

  void next_static() {
    if (m_league_rank < m_league_end) {
      // Make sure all stores are complete before entering the barrier
      memory_fence();
      team_barrier();
      set_team_shared();
    }
    m_league_rank++;
  }

  bool valid_dynamic() {
    if (m_invalid_thread) return false;
    if ((m_league_rank < m_league_chunk_end) &&
        (m_league_rank < m_league_size)) {
      return true;
    }

    if (m_team_rank_rev == 0) {
      m_team_base[0]->get_work_index(m_team_alloc);
    }
    team_barrier();

    long work_index = m_team_base[0]->team_work_index();

    m_league_rank      = work_index * m_chunk_size;
    m_league_chunk_end = (work_index + 1) * m_chunk_size;

    if (m_league_chunk_end > m_league_size) m_league_chunk_end = m_league_size;

    if ((m_league_rank >= 0) && (m_league_rank < m_league_chunk_end))
      return true;
    return false;
  }

  void next_dynamic() {
    if (m_invalid_thread) return;

    if (m_league_rank < m_league_chunk_end) {
      // Make sure all stores are complete before entering the barrier
      memory_fence();
      team_barrier();
      set_team_shared();
    }
    m_league_rank++;
  }

  void set_league_shmem(const int arg_league_rank, const int arg_league_size,
                        const size_t arg_shmem_size) {
    m_league_rank      = arg_league_rank;
    m_league_size      = arg_league_size;
    m_team_shared_size = arg_shmem_size;
    set_team_shared();
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
template <class... Properties>
class TeamPolicyInternal<Kokkos::Threads, Properties...>
    : public PolicyTraits<Properties...> {
 private:
  int m_league_size;
  int m_team_size;
  int m_team_alloc;
  int m_team_iter;

  size_t m_team_scratch_size[2];
  size_t m_thread_scratch_size[2];

  int m_chunk_size;

  bool m_tune_team_size;
  bool m_tune_vector_length;

  inline void init(const int league_size_request, const int team_size_request) {
    const int pool_size = traits::execution_space::impl_thread_pool_size(0);
    const int max_host_team_size = Impl::HostThreadTeamData::max_team_members;
    const int team_max =
        pool_size < max_host_team_size ? pool_size : max_host_team_size;
    const int team_grain = traits::execution_space::impl_thread_pool_size(2);

    m_league_size = league_size_request;

    if (team_size_request > team_max)
      Kokkos::abort("Kokkos::abort: Requested Team Size is too large!");

    m_team_size = team_size_request < team_max ? team_size_request : team_max;

    // Round team size up to a multiple of 'team_gain'
    const int team_size_grain =
        (m_team_size + team_grain - 1 <= 0)
            ? 1
            : team_grain * ((m_team_size + team_grain - 1) / team_grain);
    const int team_count = pool_size / team_size_grain;

    // Constraint : pool_size = m_team_alloc * team_count
    m_team_alloc = pool_size / team_count;

    // Maxumum number of iterations each team will take:
    m_team_iter = (m_league_size + team_count - 1) / team_count;

    set_auto_chunk_size();
  }

 public:
  //! Tag this class as a kokkos execution policy
  //! Tag this class as a kokkos execution policy
  using execution_policy = TeamPolicyInternal;

  using traits = PolicyTraits<Properties...>;

  const typename traits::execution_space& space() const {
    static typename traits::execution_space m_space;
    return m_space;
  }

  template <class ExecSpace, class... OtherProperties>
  friend class TeamPolicyInternal;

  template <class... OtherProperties>
  TeamPolicyInternal(
      const TeamPolicyInternal<Kokkos::Threads, OtherProperties...>& p) {
    m_league_size            = p.m_league_size;
    m_team_size              = p.m_team_size;
    m_team_alloc             = p.m_team_alloc;
    m_team_iter              = p.m_team_iter;
    m_team_scratch_size[0]   = p.m_team_scratch_size[0];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_team_scratch_size[1]   = p.m_team_scratch_size[1];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size             = p.m_chunk_size;
    m_tune_team_size         = p.m_tune_team_size;
    m_tune_vector_length     = p.m_tune_vector_length;
  }

  //----------------------------------------

  template <class FunctorType>
  int team_size_max(const FunctorType&, const ParallelForTag&) const {
    int pool_size          = traits::execution_space::impl_thread_pool_size(1);
    int max_host_team_size = Impl::HostThreadTeamData::max_team_members;
    return pool_size < max_host_team_size ? pool_size : max_host_team_size;
  }
  template <class FunctorType>
  int team_size_max(const FunctorType&, const ParallelReduceTag&) const {
    int pool_size          = traits::execution_space::impl_thread_pool_size(1);
    int max_host_team_size = Impl::HostThreadTeamData::max_team_members;
    return pool_size < max_host_team_size ? pool_size : max_host_team_size;
  }
  template <class FunctorType, class ReducerType>
  inline int team_size_max(const FunctorType& f, const ReducerType&,
                           const ParallelReduceTag& t) const {
    return team_size_max(f, t);
  }
  template <class FunctorType>
  int team_size_recommended(const FunctorType&, const ParallelForTag&) const {
    return traits::execution_space::impl_thread_pool_size(2);
  }
  template <class FunctorType>
  int team_size_recommended(const FunctorType&,
                            const ParallelReduceTag&) const {
    return traits::execution_space::impl_thread_pool_size(2);
  }
  template <class FunctorType, class ReducerType>
  inline int team_size_recommended(const FunctorType& f, const ReducerType&,
                                   const ParallelReduceTag& t) const {
    return team_size_recommended(f, t);
  }

  inline static int vector_length_max() {
    return 1024;
  }  // Use arbitrary large number, is meant as a vectorizable length

  inline static int scratch_size_max(int level) {
    return (level == 0 ? 1024 * 32 :  // Roughly L1 size
                20 * 1024 * 1024);    // Limit to keep compatibility with CUDA
  }

  //----------------------------------------

  inline int team_size() const { return m_team_size; }
  inline int impl_vector_length() const { return 1; }
  inline size_t team_alloc() const { return m_team_alloc; }
  inline int league_size() const { return m_league_size; }

  inline bool impl_auto_team_size() const { return m_tune_team_size; }
  inline bool impl_auto_vector_length() const { return m_tune_vector_length; }
  inline void impl_set_team_size(size_t size) { init(m_league_size, size); }
  inline void impl_set_vector_length(size_t /**size*/) {}
  inline size_t scratch_size(const int& level, int team_size_ = -1) const {
    if (team_size_ < 0) team_size_ = m_team_size;
    return m_team_scratch_size[level] +
           team_size_ * m_thread_scratch_size[level];
  }

  inline int team_iter() const { return m_team_iter; }

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request, int team_size_request,
                     int vector_length_request = 1)
      : m_league_size(0),
        m_team_size(0),
        m_team_alloc(0),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(0),
        m_tune_team_size(false),
        m_tune_vector_length(false) {
    init(league_size_request, team_size_request);
    (void)vector_length_request;
  }

  /** \brief  Specify league size, request team size and vector length*/
  TeamPolicyInternal(const typename traits::execution_space& space,
                     int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : TeamPolicyInternal(space, league_size_request, -1, -1) {}

  /** \brief  Specify league size, request team size*/
  TeamPolicyInternal(const typename traits::execution_space& space,
                     int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     int vector_length_request)
      : TeamPolicyInternal(space, league_size_request, -1,
                           vector_length_request) {}

  /** \brief  Specify league size and team size, request vector length*/
  TeamPolicyInternal(const typename traits::execution_space& space,
                     int league_size_request, int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : TeamPolicyInternal(space, league_size_request, team_size_request, -1) {}

  TeamPolicyInternal(int league_size_request, int team_size_request,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(),
                           league_size_request, team_size_request,
                           vector_length_request) {}

  TeamPolicyInternal(int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(),
                           league_size_request, -1, vector_length_request) {}

  /** \brief  Specify league size, request team size and vector length*/
  TeamPolicyInternal(int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : TeamPolicyInternal(typename traits::execution_space(),
                           league_size_request, -1, -1) {}

  /** \brief  Specify league size and team size, request vector length*/
  TeamPolicyInternal(int league_size_request, int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : TeamPolicyInternal(typename traits::execution_space(),
                           league_size_request, team_size_request, -1) {}

  inline int chunk_size() const { return m_chunk_size; }

  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal& set_chunk_size(
      typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal& set_scratch_size(const int& level,
                                              const PerTeamValue& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal& set_scratch_size(
      const int& level, const PerThreadValue& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of
   * the scratch hierarchy */
  inline TeamPolicyInternal& set_scratch_size(
      const int& level, const PerTeamValue& per_team,
      const PerThreadValue& per_thread) {
    m_team_scratch_size[level]   = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

 private:
  /** \brief finalize chunk_size if it was set to AUTO*/
  inline void set_auto_chunk_size() {
    int64_t concurrency = space().concurrency() / m_team_alloc;
    if (concurrency == 0) concurrency = 1;

    if (m_chunk_size > 0) {
      if (!Impl::is_integral_power_of_two(m_chunk_size))
        Kokkos::abort("TeamPolicy blocking granularity must be power of two");
    }

    int new_chunk_size = 1;
    while (new_chunk_size * 100 * concurrency < m_league_size)
      new_chunk_size *= 2;
    if (new_chunk_size < 128) {
      new_chunk_size = 1;
      while ((new_chunk_size * 40 * concurrency < m_league_size) &&
             (new_chunk_size < 128))
        new_chunk_size *= 2;
    }
    m_chunk_size = new_chunk_size;
  }

 public:
  using member_type = Impl::ThreadsExecTeamMember;

  friend class Impl::ThreadsExecTeamMember;
};

} /*namespace Impl */
} /* namespace Kokkos */

namespace Kokkos {

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamThreadRangeBoundariesStruct<iType, Impl::ThreadsExecTeamMember>
    TeamThreadRange(const Impl::ThreadsExecTeamMember& thread,
                    const iType& count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType,
                                               Impl::ThreadsExecTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::ThreadsExecTeamMember>
TeamThreadRange(const Impl::ThreadsExecTeamMember& thread, const iType1& begin,
                const iType2& end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamThreadRangeBoundariesStruct<iType,
                                               Impl::ThreadsExecTeamMember>(
      thread, iType(begin), iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamThreadRangeBoundariesStruct<iType, Impl::ThreadsExecTeamMember>
    TeamVectorRange(const Impl::ThreadsExecTeamMember& thread,
                    const iType& count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType,
                                               Impl::ThreadsExecTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::ThreadsExecTeamMember>
TeamVectorRange(const Impl::ThreadsExecTeamMember& thread, const iType1& begin,
                const iType2& end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamThreadRangeBoundariesStruct<iType,
                                               Impl::ThreadsExecTeamMember>(
      thread, iType(begin), iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::ThreadsExecTeamMember>
    ThreadVectorRange(const Impl::ThreadsExecTeamMember& thread,
                      const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,
                                                 Impl::ThreadsExecTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::ThreadsExecTeamMember>
ThreadVectorRange(const Impl::ThreadsExecTeamMember& thread,
                  const iType1& arg_begin, const iType2& arg_end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::ThreadVectorRangeBoundariesStruct<iType,
                                                 Impl::ThreadsExecTeamMember>(
      thread, iType(arg_begin), iType(arg_end));
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::ThreadsExecTeamMember> PerTeam(
    const Impl::ThreadsExecTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::ThreadsExecTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::ThreadsExecTeamMember> PerThread(
    const Impl::ThreadsExecTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::ThreadsExecTeamMember>(thread);
}
}  // namespace Kokkos

namespace Kokkos {

/** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team.
 */
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::ThreadsExecTeamMember>& loop_boundaries,
    const Lambda& lambda) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team
 * and a summation of val is performed and put into result.
 */
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer<ValueType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::ThreadsExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ValueType& result) {
  ValueType intermediate;
  Sum<ValueType> sum(intermediate);
  sum.init(intermediate);

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i, tmp);
    intermediate += tmp;
  }

  loop_boundaries.thread.team_reduce(sum, intermediate);
  result = sum.reference();
}

template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::ThreadsExecTeamMember>& loop_boundaries,
                const Lambda& lambda, const ReducerType& reducer) {
  typename ReducerType::value_type value;
  reducer.init(value);

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, value);
  }

  loop_boundaries.thread.team_reduce(reducer, value);
}

}  // namespace Kokkos

namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 */
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::ThreadsExecTeamMember>& loop_boundaries,
    const Lambda& lambda) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment)
    lambda(i);
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a summation of val is performed and put into result.
 */
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer<ValueType>::value>
parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<
                    iType, Impl::ThreadsExecTeamMember>& loop_boundaries,
                const Lambda& lambda, ValueType& result) {
  result = ValueType();
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }
}

template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<
                    iType, Impl::ThreadsExecTeamMember>& loop_boundaries,
                const Lambda& lambda, const ReducerType& reducer) {
  reducer.init(reducer.reference());
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, reducer.reference());
  }
}

/** \brief  Inter-thread parallel exclusive prefix sum. Executes
 * lambda(iType i, ValueType & val, bool final) for each i=0..N-1.
 *
 */
template <typename iType, class FunctorType, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::ThreadsExecTeamMember>& loop_bounds,
    const FunctorType& lambda, ValueType& return_val) {
  // Extract ValueType from the Closure
  using closure_value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, FunctorType,
      void>::value_type;
  static_assert(std::is_same_v<closure_value_type, ValueType>,
                "Non-matching value types of closure and return type");

  auto scan_val = ValueType{};

  // Intra-member scan
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_bounds.start; i < loop_bounds.end;
       i += loop_bounds.increment) {
    lambda(i, scan_val, false);
  }

  auto& team_member = loop_bounds.thread;

  // 'scan_val' output is the exclusive prefix sum
  scan_val = team_member.team_scan(scan_val);

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_bounds.start; i < loop_bounds.end;
       i += loop_bounds.increment) {
    lambda(i, scan_val, true);
  }

  team_member.team_broadcast(scan_val, team_member.team_size() - 1);

  return_val = scan_val;
}

template <typename iType, class FunctorType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::ThreadsExecTeamMember>& loop_bounds,
    const FunctorType& lambda) {
  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, FunctorType,
      void>::value_type;

  value_type scan_val;
  parallel_scan(loop_bounds, lambda, scan_val);
}

/** \brief  Intra-thread vector parallel exclusive prefix sum. Executes
 * lambda(iType i, ValueType & val, bool final) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes in the thread and a scan
 * operation is performed. Depending on the target execution space the operator
 * might be called twice: once with final=false and once with final=true. When
 * final==true val contains the prefix sum value. The contribution of this "i"
 * needs to be added to val no matter whether final==true or not. In a serial
 * execution (i.e. team_size==1) the operator is only called once with
 * final==true. Scan_val will be set to the final sum value over all vector
 * lanes.
 */
template <typename iType, class FunctorType, typename ValueType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::ThreadsExecTeamMember>& loop_boundaries,
    const FunctorType& lambda, ValueType& return_val) {
  // Extract ValueType from the Closure
  using closure_value_type =
      typename Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                     TeamPolicy<Threads>, FunctorType,
                                     void>::value_type;
  static_assert(std::is_same<closure_value_type, ValueType>::value,
                "Non-matching value types of closure and return type");

  ValueType scan_val = ValueType();

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, scan_val, true);
  }

  return_val = scan_val;
}

template <typename iType, class FunctorType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::ThreadsExecTeamMember>& loop_boundaries,
    const FunctorType& lambda) {
  using value_type =
      typename Impl::FunctorAnalysis<Impl::FunctorPatternInterface::SCAN,
                                     TeamPolicy<Threads>, FunctorType,
                                     void>::value_type;

  value_type scan_val;
  parallel_scan(loop_boundaries, lambda, scan_val);
}

/** \brief  Intra-thread vector parallel scan with reducer
 *
 */
template <typename iType, class FunctorType, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<
                  iType, Impl::ThreadsExecTeamMember>& loop_boundaries,
              const FunctorType& lambda, const ReducerType& reducer) {
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

}  // namespace Kokkos

namespace Kokkos {

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<
        Impl::ThreadsExecTeamMember>& /*single_struct*/,
    const FunctorType& lambda) {
  lambda();
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::ThreadsExecTeamMember>& single_struct,
    const FunctorType& lambda) {
  if (single_struct.team_member.team_rank() == 0) lambda();
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<
        Impl::ThreadsExecTeamMember>& /*single_struct*/,
    const FunctorType& lambda, ValueType& val) {
  lambda(val);
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::ThreadsExecTeamMember>& single_struct,
    const FunctorType& lambda, ValueType& val) {
  if (single_struct.team_member.team_rank() == 0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val, 0);
}
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#endif /* #define KOKKOS_THREADSTEAM_HPP */
