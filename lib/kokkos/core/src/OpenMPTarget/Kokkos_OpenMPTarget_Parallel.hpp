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

#ifndef KOKKOS_OPENMPTARGET_PARALLEL_HPP
#define KOKKOS_OPENMPTARGET_PARALLEL_HPP

#include <omp.h>
#include <sstream>
#include <Kokkos_Parallel.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Spinwait.hpp>

#include <Kokkos_Atomic.hpp>
#include "Kokkos_OpenMPTarget_Abort.hpp"

// Intel architectures prefer the classical hierarchical parallelism that relies
// on OpenMP.
#if defined(KOKKOS_ARCH_INTEL_GPU)
#define KOKKOS_IMPL_OPENMPTARGET_HIERARCHICAL_INTEL_GPU
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

class OpenMPTargetExecTeamMember {
 public:
  static constexpr int TEAM_REDUCE_SIZE = 512;

  using execution_space      = Kokkos::Experimental::OpenMPTarget;
  using scratch_memory_space = execution_space::scratch_memory_space;
  using team_handle          = OpenMPTargetExecTeamMember;

  scratch_memory_space m_team_shared;
  size_t m_team_scratch_size[2];
  int m_team_rank;
  int m_team_size;
  int m_league_rank;
  int m_league_size;
  int m_vector_length;
  int m_vector_lane;
  int m_shmem_block_index;
  void* m_glb_scratch;
  void* m_reduce_scratch;

 public:
  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_shmem() const {
    return m_team_shared.set_team_thread_mode(0, 1, 0);
  }

  // set_team_thread_mode routine parameters for future understanding:
  // first parameter - scratch level.
  // second parameter - size multiplier for advancing scratch ptr after a
  // request was serviced. third parameter - offset size multiplier from current
  // scratch ptr when returning a ptr for a request.
  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_scratch(int level) const {
    return m_team_shared.set_team_thread_mode(level, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& thread_scratch(int level) const {
    return m_team_shared.set_team_thread_mode(level, team_size(), team_rank());
  }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size; }
  KOKKOS_INLINE_FUNCTION int team_rank() const { return m_team_rank; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size; }
  KOKKOS_INLINE_FUNCTION void* impl_reduce_scratch() const {
    return m_reduce_scratch;
  }

  KOKKOS_INLINE_FUNCTION void team_barrier() const {
#pragma omp barrier
  }

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(ValueType& value,
                                             int thread_id) const {
    // Make sure there is enough scratch space:
    using type = std::conditional_t<(sizeof(ValueType) < TEAM_REDUCE_SIZE),
                                    ValueType, void>;
    type* team_scratch =
        reinterpret_cast<type*>(static_cast<char*>(m_glb_scratch) +
                                TEAM_REDUCE_SIZE * omp_get_team_num());
#pragma omp barrier
    if (team_rank() == thread_id) *team_scratch = value;
#pragma omp barrier
    value = *team_scratch;
  }

  template <class Closure, class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(const Closure& f, ValueType& value,
                                             const int& thread_id) const {
    f(value);
    team_broadcast(value, thread_id);
  }

  // FIXME_OPENMPTARGET this function has the wrong interface and currently
  // ignores the reducer passed.
  template <class ValueType, class JoinOp>
  KOKKOS_INLINE_FUNCTION ValueType team_reduce(const ValueType& value,
                                               const JoinOp&) const {
#pragma omp barrier

    using value_type = ValueType;
    //    const JoinLambdaAdapter<value_type, JoinOp> op(op_in);

    // Make sure there is enough scratch space:
    using type = std::conditional_t<(sizeof(value_type) < TEAM_REDUCE_SIZE),
                                    value_type, void>;

    const int n_values = TEAM_REDUCE_SIZE / sizeof(value_type);
    type* team_scratch =
        reinterpret_cast<type*>(static_cast<char*>(m_glb_scratch) +
                                TEAM_REDUCE_SIZE * omp_get_team_num());
    for (int i = m_team_rank; i < n_values; i += m_team_size) {
      team_scratch[i] = value_type();
    }

#pragma omp barrier

    for (int k = 0; k < m_team_size; k += n_values) {
      if ((k <= m_team_rank) && (k + n_values > m_team_rank))
        team_scratch[m_team_rank % n_values] += value;
#pragma omp barrier
    }

    for (int d = 1; d < n_values; d *= 2) {
      if ((m_team_rank + d < n_values) && (m_team_rank % (2 * d) == 0)) {
        team_scratch[m_team_rank] += team_scratch[m_team_rank + d];
      }
#pragma omp barrier
    }
    return team_scratch[0];
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
  KOKKOS_INLINE_FUNCTION ArgType
  team_scan(const ArgType& /*value*/, ArgType* const /*global_accum*/) const {
    // FIXME_OPENMPTARGET
    /*  // Make sure there is enough scratch space:
      using type =
        std::conditional_t<(sizeof(ArgType) < TEAM_REDUCE_SIZE), ArgType, void>;

      volatile type * const work_value  = ((type*) m_exec.scratch_thread());

      *work_value = value ;

      memory_fence();

      if ( team_fan_in() ) {
        // The last thread to synchronize returns true, all other threads wait
      for team_fan_out()
        // m_team_base[0]                 == highest ranking team member
        // m_team_base[ m_team_size - 1 ] == lowest ranking team member
        //
        // 1) copy from lower to higher rank, initialize lowest rank to zero
        // 2) prefix sum from lowest to highest rank, skipping lowest rank

        type accum = 0 ;

        if ( global_accum ) {
          for ( int i = m_team_size ; i-- ; ) {
            type & val = *((type*) m_exec.pool_rev( m_team_base_rev + i
      )->scratch_thread()); accum += val ;
          }
          accum = atomic_fetch_add( global_accum , accum );
        }

        for ( int i = m_team_size ; i-- ; ) {
          type & val = *((type*) m_exec.pool_rev( m_team_base_rev + i
      )->scratch_thread()); const type offset = accum ; accum += val ; val =
      offset ;
        }

        memory_fence();
      }

      team_fan_out();

      return *work_value ;*/
    return ArgType();
  }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template <typename Type>
  KOKKOS_INLINE_FUNCTION Type team_scan(const Type& value) const {
    return this->template team_scan<Type>(value, 0);
  }

  //----------------------------------------
  // Private for the driver

 private:
  using space = execution_space::scratch_memory_space;

 public:
  // FIXME_OPENMPTARGET - 512(16*32) bytes at the begining of the scratch space
  // for each league is saved for reduction. It should actually be based on the
  // ValueType of the reduction variable.
  inline OpenMPTargetExecTeamMember(
      const int league_rank, const int league_size, const int team_size,
      const int vector_length  // const TeamPolicyInternal< OpenMPTarget,
                               // Properties ...> & team
      ,
      void* const glb_scratch, const int shmem_block_index,
      const size_t shmem_size_L0, const size_t shmem_size_L1)
      : m_team_scratch_size{shmem_size_L0, shmem_size_L1},
        m_team_rank(0),
        m_team_size(team_size),
        m_league_rank(league_rank),
        m_league_size(league_size),
        m_vector_length(vector_length),
        m_shmem_block_index(shmem_block_index),
        m_glb_scratch(glb_scratch) {
    const int omp_tid = omp_get_thread_num();

    // The scratch memory allocated is a sum of TEAM_REDUCE_SIZE, L0 shmem size
    // and L1 shmem size. TEAM_REDUCE_SIZE = 512 bytes saved per team for
    // hierarchical reduction. There is an additional 10% of the requested
    // scratch memory allocated per team as padding. Hence the product with 0.1.
    const int reduce_offset =
        m_shmem_block_index *
        (shmem_size_L0 + shmem_size_L1 +
         ((shmem_size_L0 + shmem_size_L1) * 0.1) + TEAM_REDUCE_SIZE);
    const int l0_offset = reduce_offset + TEAM_REDUCE_SIZE;
    const int l1_offset = l0_offset + shmem_size_L0;
    m_team_shared       = scratch_memory_space(
        (static_cast<char*>(glb_scratch) + l0_offset), shmem_size_L0,
        static_cast<char*>(glb_scratch) + l1_offset, shmem_size_L1);
    m_reduce_scratch = static_cast<char*>(glb_scratch) + reduce_offset;
    m_league_rank    = league_rank;
    m_team_rank      = omp_tid;
    m_vector_lane    = 0;
  }

  static inline int team_reduce_size() { return TEAM_REDUCE_SIZE; }
};

template <class... Properties>
class TeamPolicyInternal<Kokkos::Experimental::OpenMPTarget, Properties...>
    : public PolicyTraits<Properties...> {
 public:
  //! Tag this class as a kokkos execution policy
  using execution_policy = TeamPolicyInternal;

  using traits = PolicyTraits<Properties...>;

  //----------------------------------------

  template <class FunctorType>
  inline static int team_size_max(const FunctorType&, const ParallelForTag&) {
    return 256;
  }

  template <class FunctorType>
  inline static int team_size_max(const FunctorType&,
                                  const ParallelReduceTag&) {
    return 256;
  }

  template <class FunctorType, class ReducerType>
  inline static int team_size_max(const FunctorType&, const ReducerType&,
                                  const ParallelReduceTag&) {
    return 256;
  }

  template <class FunctorType>
  inline static int team_size_recommended(const FunctorType&,
                                          const ParallelForTag&) {
    return 128;
  }

  template <class FunctorType>
  inline static int team_size_recommended(const FunctorType&,
                                          const ParallelReduceTag&) {
    return 128;
  }

  template <class FunctorType, class ReducerType>
  inline static int team_size_recommended(const FunctorType&,
                                          const ReducerType&,
                                          const ParallelReduceTag&) {
    return 128;
  }

  //----------------------------------------

 private:
  int m_league_size;
  int m_team_size;
  int m_vector_length;
  int m_team_alloc;
  int m_team_iter;
  std::array<size_t, 2> m_team_scratch_size;
  std::array<size_t, 2> m_thread_scratch_size;
  bool m_tune_team_size;
  bool m_tune_vector_length;
  constexpr const static size_t default_team_size = 256;
  int m_chunk_size;

  inline void init(const int league_size_request, const int team_size_request,
                   const int vector_length_request) {
    m_league_size = league_size_request;

    // Minimum team size should be 32 for OpenMPTarget backend.
    if (team_size_request < 32) {
      Kokkos::Impl::OpenMPTarget_abort(
          "OpenMPTarget backend requires a minimum of 32 threads per team.\n");
    } else
      m_team_size = team_size_request;

    m_vector_length = vector_length_request;
    set_auto_chunk_size();
  }

  template <typename ExecSpace, typename... OtherProperties>
  friend class TeamPolicyInternal;

 public:
  // FIXME_OPENMPTARGET : Currently this routine is a copy of the Cuda
  // implementation, but this has to be tailored to be architecture specific.
  inline static int scratch_size_max(int level) {
    return (
        level == 0 ? 1024 * 40 :  // 48kB is the max for CUDA, but we need some
                                  // for team_member.reduce etc.
            20 * 1024 *
                1024);  // arbitrarily setting this to 20MB, for a Volta V100
                        // that would give us about 3.2GB for 2 teams per SM
  }
  inline bool impl_auto_team_size() const { return m_tune_team_size; }
  inline bool impl_auto_vector_length() const { return m_tune_vector_length; }
  inline void impl_set_team_size(const size_t size) { m_team_size = size; }
  inline void impl_set_vector_length(const size_t length) {
    m_tune_vector_length = length;
  }
  inline int impl_vector_length() const { return m_vector_length; }
  inline int team_size() const { return m_team_size; }
  inline int league_size() const { return m_league_size; }
  inline size_t scratch_size(const int& level, int team_size_ = -1) const {
    if (team_size_ < 0) team_size_ = m_team_size;
    return m_team_scratch_size[level] +
           team_size_ * m_thread_scratch_size[level];
  }

  inline Kokkos::Experimental::OpenMPTarget space() const {
    return Kokkos::Experimental::OpenMPTarget();
  }

  template <class... OtherProperties>
  TeamPolicyInternal(const TeamPolicyInternal<OtherProperties...>& p)
      : m_league_size(p.m_league_size),
        m_team_size(p.m_team_size),
        m_vector_length(p.m_vector_length),
        m_team_alloc(p.m_team_alloc),
        m_team_iter(p.m_team_iter),
        m_team_scratch_size(p.m_team_scratch_size),
        m_thread_scratch_size(p.m_thread_scratch_size),
        m_tune_team_size(p.m_tune_team_size),
        m_tune_vector_length(p.m_tune_vector_length),
        m_chunk_size(p.m_chunk_size) {}

  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request, int team_size_request,
                     int vector_length_request = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(false),
        m_tune_vector_length(false),
        m_chunk_size(0) {
    init(league_size_request, team_size_request, vector_length_request);
  }

  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     int vector_length_request = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(true),
        m_tune_vector_length(false),
        m_chunk_size(0) {
    init(league_size_request, default_team_size / vector_length_request,
         vector_length_request);
  }

  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(true),
        m_tune_vector_length(true),
        m_chunk_size(0) {
    init(league_size_request, default_team_size, 1);
  }
  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request, int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(false),
        m_tune_vector_length(true),
        m_chunk_size(0) {
    init(league_size_request, team_size_request, 1);
  }

  TeamPolicyInternal(int league_size_request, int team_size_request,
                     int vector_length_request = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(false),
        m_tune_vector_length(false),
        m_chunk_size(0) {
    init(league_size_request, team_size_request, vector_length_request);
  }

  TeamPolicyInternal(int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     int vector_length_request = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(true),
        m_tune_vector_length(false),
        m_chunk_size(0) {
    init(league_size_request, default_team_size / vector_length_request,
         vector_length_request);
  }

  TeamPolicyInternal(int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(true),
        m_tune_vector_length(true),
        m_chunk_size(0) {
    init(league_size_request, default_team_size, 1);
  }
  TeamPolicyInternal(int league_size_request, int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_tune_team_size(false),
        m_tune_vector_length(true),
        m_chunk_size(0) {
    init(league_size_request, team_size_request, 1);
  }
  inline static size_t vector_length_max() {
    return 32; /* TODO: this is bad. Need logic that is compiler and backend
                  aware */
  }
  inline int team_alloc() const { return m_team_alloc; }
  inline int team_iter() const { return m_team_iter; }

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
    int concurrency = 2048 * 128;

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
  using member_type = Impl::OpenMPTargetExecTeamMember;
};

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {

template <typename iType>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    iType, Impl::OpenMPTargetExecTeamMember>
TeamThreadRange(const Impl::OpenMPTargetExecTeamMember& thread,
                const iType& count) {
  return Impl::TeamThreadRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::OpenMPTargetExecTeamMember>
TeamThreadRange(const Impl::OpenMPTargetExecTeamMember& thread,
                const iType1& begin, const iType2& end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamThreadRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, iType(begin),
                                               iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<
    iType, Impl::OpenMPTargetExecTeamMember>
ThreadVectorRange(const Impl::OpenMPTargetExecTeamMember& thread,
                  const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::OpenMPTargetExecTeamMember>
ThreadVectorRange(const Impl::OpenMPTargetExecTeamMember& thread,
                  const iType1& arg_begin, const iType2& arg_end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::ThreadVectorRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, iType(arg_begin),
                                               iType(arg_end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION Impl::TeamVectorRangeBoundariesStruct<
    iType, Impl::OpenMPTargetExecTeamMember>
TeamVectorRange(const Impl::OpenMPTargetExecTeamMember& thread,
                const iType& count) {
  return Impl::TeamVectorRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamVectorRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::OpenMPTargetExecTeamMember>
TeamVectorRange(const Impl::OpenMPTargetExecTeamMember& thread,
                const iType1& arg_begin, const iType2& arg_end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamVectorRangeBoundariesStruct<
      iType, Impl::OpenMPTargetExecTeamMember>(thread, iType(arg_begin),
                                               iType(arg_end));
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember> PerTeam(
    const Impl::OpenMPTargetExecTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember> PerThread(
    const Impl::OpenMPTargetExecTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember>(thread);
}
}  // namespace Kokkos

namespace Kokkos {

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember>&
    /*single_struct*/,
    const FunctorType& lambda) {
  lambda();
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember>&
        single_struct,
    const FunctorType& lambda) {
  if (single_struct.team_member.team_rank() == 0) lambda();
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::OpenMPTargetExecTeamMember>&
    /*single_struct*/,
    const FunctorType& lambda, ValueType& val) {
  lambda(val);
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::OpenMPTargetExecTeamMember>&
        single_struct,
    const FunctorType& lambda, ValueType& val) {
  if (single_struct.team_member.team_rank() == 0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val, 0);
}
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <typename iType>
struct TeamThreadRangeBoundariesStruct<iType, OpenMPTargetExecTeamMember> {
  using index_type = iType;
  const iType start;
  const iType end;
  const OpenMPTargetExecTeamMember& team;

  TeamThreadRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                  iType count)
      : start(0), end(count), team(thread_) {}
  TeamThreadRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                  iType begin_, iType end_)
      : start(begin_), end(end_), team(thread_) {}
};

template <typename iType>
struct ThreadVectorRangeBoundariesStruct<iType, OpenMPTargetExecTeamMember> {
  using index_type = iType;
  const index_type start;
  const index_type end;
  const OpenMPTargetExecTeamMember& team;

  ThreadVectorRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                    index_type count)
      : start(0), end(count), team(thread_) {}
  ThreadVectorRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                    index_type begin_, index_type end_)
      : start(begin_), end(end_), team(thread_) {}
};

template <typename iType>
struct TeamVectorRangeBoundariesStruct<iType, OpenMPTargetExecTeamMember> {
  using index_type = iType;
  const index_type start;
  const index_type end;
  const OpenMPTargetExecTeamMember& team;

  TeamVectorRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                  index_type count)
      : start(0), end(count), team(thread_) {}
  TeamVectorRangeBoundariesStruct(const OpenMPTargetExecTeamMember& thread_,
                                  index_type begin_, index_type end_)
      : start(begin_), end(end_), team(thread_) {}
};

}  // namespace Impl

}  // namespace Kokkos
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------
/** \brief  Data for OpenMPTarget thread execution */

class OpenMPTargetExec {
 public:
  // FIXME_OPENMPTARGET - Currently the maximum number of
  // teams possible is calculated based on NVIDIA's Volta GPU. In
  // future this value should be based on the chosen architecture for the
  // OpenMPTarget backend.
  static int MAX_ACTIVE_THREADS;

 private:
  static void* scratch_ptr;

 public:
  static void verify_is_process(const char* const);
  static void verify_initialized(const char* const);

  static int* get_lock_array(int num_teams);
  static void* get_scratch_ptr();
  static void clear_scratch();
  static void clear_lock_array();
  static void resize_scratch(int64_t team_reduce_bytes,
                             int64_t team_shared_bytes,
                             int64_t thread_local_bytes, int64_t league_size);

  static void* m_scratch_ptr;
  static int64_t m_scratch_size;
  static int* m_lock_array;
  static uint64_t m_lock_size;
  static uint32_t* m_uniquetoken_ptr;
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* KOKKOS_OPENMPTARGET_PARALLEL_HPP */
