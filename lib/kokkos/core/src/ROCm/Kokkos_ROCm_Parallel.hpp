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

#include <algorithm>
#include <typeinfo>
#include <ROCm/Kokkos_ROCm_Reduce.hpp>
#include <ROCm/Kokkos_ROCm_Scan.hpp>
#include <ROCm/Kokkos_ROCm_Exec.hpp>
#include <ROCm/Kokkos_ROCm_Vectorization.hpp>
#include <ROCm/KokkosExp_ROCm_IterateTile_Refactor.hpp>

#include <KokkosExp_MDRangePolicy.hpp>

namespace Kokkos {
namespace Impl {

struct ROCmTeamMember;

template <class... Properties>
class TeamPolicyInternal<Kokkos::Experimental::ROCm, Properties...>
    : public PolicyTraits<Properties...> {
 private:
  int m_league_size;
  int m_team_size;
  int m_vector_length;
  int m_team_scratch_size[2];
  int m_thread_scratch_size[2];
  int m_chunk_size;

 public:
  using execution_policy = TeamPolicyInternal;
  using execution_space  = Kokkos::Experimental::ROCm;
  typedef PolicyTraits<Properties...> traits;

  TeamPolicyInternal& operator=(const TeamPolicyInternal& p) {
    m_league_size            = p.m_league_size;
    m_team_size              = p.m_team_size;
    m_vector_length          = p.m_vector_length;
    m_team_scratch_size[0]   = p.m_team_scratch_size[0];
    m_team_scratch_size[1]   = p.m_team_scratch_size[1];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size             = p.m_chunk_size;
    return *this;
  }

  template <class ExecSpace, class... OtherProperties>
  friend class TeamPolicyInternal;

  template <class... OtherProperties>
  TeamPolicyInternal(const TeamPolicyInternal<Kokkos::Experimental::ROCm,
                                              OtherProperties...>& p) {
    m_league_size            = p.m_league_size;
    m_team_size              = p.m_team_size;
    m_vector_length          = p.m_vector_length;
    m_team_scratch_size[0]   = p.m_team_scratch_size[0];
    m_team_scratch_size[1]   = p.m_team_scratch_size[1];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size             = p.m_chunk_size;
  }

  TeamPolicyInternal()
      : m_league_size(0),
        m_team_size(0),
        m_vector_length(0),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(64) {}

  TeamPolicyInternal(const int arg_league_size, const int arg_team_size)
      : m_league_size(arg_league_size),
        m_team_size(arg_team_size),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(64) {}

  TeamPolicyInternal(const int arg_league_size, const int arg_team_size,
                     const int vector_length_request = 1)
      : m_league_size(arg_league_size),
        m_team_size(arg_team_size),
        m_vector_length(vector_length_request),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(64) {}

  TeamPolicyInternal(const int arg_league_size, const Kokkos::AUTO_t)
      : m_league_size(arg_league_size),
        m_team_size(-1),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(64) {}

  TeamPolicyInternal(const int arg_league_size, const Kokkos::AUTO_t,
                     const int vector_length_request)
      : m_league_size(arg_league_size),
        m_team_size(-1),
        m_vector_length(vector_length_request),
        m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_chunk_size(64) {}

  inline int chunk_size() const { return m_chunk_size; }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  /** \brief set chunk_size to a discrete value*/
  KOKKOS_INLINE_FUNCTION TeamPolicyInternal
  set_chunk_size(typename traits::index_type chunk_size_) const {
    TeamPolicyInternal p = *this;
    p.m_chunk_size       = chunk_size_;
    return p;
  }

  /** \brief set per team scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal set_scratch_size(
      const int& level, const PerTeamValue& per_team) const {
    TeamPolicyInternal p         = *this;
    p.m_team_scratch_size[level] = per_team.value;
    return p;
  }

  /** \brief set per thread scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal set_scratch_size(
      const int& level, const PerThreadValue& per_thread) const {
    TeamPolicyInternal p           = *this;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  }

  /** \brief set per thread and per team scratch size for a specific level of
   * the scratch hierarchy */
  inline TeamPolicyInternal set_scratch_size(
      const int& level, const PerTeamValue& per_team,
      const PerThreadValue& per_thread) const {
    TeamPolicyInternal p           = *this;
    p.m_team_scratch_size[level]   = per_team.value;
    p.m_thread_scratch_size[level] = per_thread.value;
    return p;
  }
#else
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
#endif

 protected:
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  /** \brief set chunk_size to a discrete value*/
  inline TeamPolicyInternal internal_set_chunk_size(
      typename traits::index_type chunk_size_) {
    m_chunk_size = chunk_size_;
    return *this;
  }

  /** \brief set per team scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal internal_set_scratch_size(
      const int& level, const PerTeamValue& per_team) {
    m_team_scratch_size[level] = per_team.value;
    return *this;
  }

  /** \brief set per thread scratch size for a specific level of the scratch
   * hierarchy */
  inline TeamPolicyInternal internal_set_scratch_size(
      const int& level, const PerThreadValue& per_thread) {
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }

  /** \brief set per thread and per team scratch size for a specific level of
   * the scratch hierarchy */
  inline TeamPolicyInternal internal_set_scratch_size(
      const int& level, const PerTeamValue& per_team,
      const PerThreadValue& per_thread) {
    m_team_scratch_size[level]   = per_team.value;
    m_thread_scratch_size[level] = per_thread.value;
    return *this;
  }
#endif

 public:
  // TODO:  evaluate proper team_size_max requirements
  template <class Functor_Type>
  KOKKOS_INLINE_FUNCTION static int team_size_max(const Functor_Type& functor) {
    typedef typename Kokkos::Impl::FunctorValueTraits<
        Functor_Type, typename traits::work_tag>::value_type value_type;
    return team_size_recommended(functor);
    // return std::min(Kokkos::Impl::get_max_tile_size() / sizeof(value_type),
    // Kokkos::Impl::get_max_tile_thread());
  }

  template <class Functor_Type>
  KOKKOS_INLINE_FUNCTION static int team_size_recommended(
      const Functor_Type& functor) {
    return Kokkos::Impl::get_tile_size<
        typename Kokkos::Impl::FunctorValueTraits<
            Functor_Type, typename traits::work_tag>::value_type>();
  }

  template <class Functor_Type>
  KOKKOS_INLINE_FUNCTION static int team_size_recommended(
      const Functor_Type& functor, const int vector_length) {
    int max = team_size_recommended(functor) / vector_length;
    if (max < 1) max = 1;
    return (max);
  }

  template <class FunctorType, class PatternTypeTag>
  int team_size_max(const FunctorType& functor, PatternTypeTag) {
    return 256 / vector_length();
  }
  template <class FunctorType, class PatternTypeTag>
  int team_size_recommended(const FunctorType& functor, PatternTypeTag) {
    return 128 / vector_length();
  }

  template <class F>
  KOKKOS_INLINE_FUNCTION int team_size(const F& f) const {
    return (m_team_size > 0) ? m_team_size : team_size_recommended(f);
  }
  KOKKOS_INLINE_FUNCTION int team_size() const {
    return (m_team_size > 0) ? m_team_size : Impl::get_max_tile_thread();
    ;
  }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size; }

  inline int vector_length() const { return m_vector_length; }
  inline int scratch_size(int level, int team_size_ = -1) const {
    if (team_size_ < 0) team_size_ = m_team_size;
    return m_team_scratch_size[level] +
           team_size_ * m_thread_scratch_size[level];
  }
  inline size_t team_scratch_size(int level) const {
    return m_team_scratch_size[level];
  }
  inline size_t thread_scratch_size(int level) const {
    return m_thread_scratch_size[level];
  }

  static int scratch_size_max(int level) {
    return level == 0 ? 1024 * 40 : 1024 * 1204 * 20;
  }

  typedef Impl::ROCmTeamMember member_type;
};

struct ROCmTeamMember {
  typedef Kokkos::Experimental::ROCm execution_space;
  typedef Kokkos::ScratchMemorySpace<Kokkos::Experimental::ROCm>
      scratch_memory_space;

  KOKKOS_INLINE_FUNCTION
  const scratch_memory_space& team_shmem() const {
    return m_team_shared.set_team_thread_mode(0, 1, 0);
  }
  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_scratch(
      const int& level) const {
    return m_team_shared.set_team_thread_mode(level, 1, 0);
  }
  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& thread_scratch(
      const int& level) const {
    return m_team_shared.set_team_thread_mode(level, team_size(), team_rank());
  }

  /* Rank of this team within the league of teams */
  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_idx.tile[0]; }
  /* Number of teams in the league */
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size; }
  /* Rank of this thread within this team */
  KOKKOS_INLINE_FUNCTION int team_rank() const {
    return m_idx.local[0] / m_vector_length;
  }
  /* Rank of this thread within this thread */
  KOKKOS_INLINE_FUNCTION int vector_rank() const {
    return m_idx.local[0] % m_vector_length;
  }
  KOKKOS_INLINE_FUNCTION int lindex() const { return m_idx.local[0]; }
  KOKKOS_INLINE_FUNCTION int gindex() const { return m_idx.global[0]; }
  KOKKOS_INLINE_FUNCTION int tindex() const { return m_idx.tile[0]; }
  KOKKOS_INLINE_FUNCTION int tile_dim() const { return m_idx.tile_dim[0]; }
  KOKKOS_INLINE_FUNCTION int team_size() const { return m_team_size; }
  KOKKOS_INLINE_FUNCTION int vector_length() const { return m_vector_length; }

  KOKKOS_INLINE_FUNCTION
  ROCmTeamMember(const hc::tiled_index<1>& arg_idx, int league_size_,
                 int team_size_)
      : m_league_size(league_size_),
        m_team_size(team_size_),
        m_team_shared(nullptr, 0),
        m_vector_length(1),
        m_idx(arg_idx) {}

  KOKKOS_INLINE_FUNCTION
  ROCmTeamMember(const hc::tiled_index<1>& arg_idx, int league_size_,
                 int team_size_, char* shared, std::size_t shsize,
                 std::size_t scratch_size0, char* scratch_ptr,
                 std::size_t scratch_size1, std::size_t vector_length)
      : m_league_size(league_size_),
        m_team_size(team_size_),
        m_team_shared(shared + arg_idx.tile[0] * (shsize + scratch_size0),
                      (shsize + scratch_size0) * league_size_,
                      scratch_ptr + arg_idx.tile[0] * scratch_size1,
                      scratch_size1 * league_size_),
        m_vector_length(vector_length),
        m_idx(arg_idx) {}

  KOKKOS_INLINE_FUNCTION
  void team_barrier() const { m_idx.barrier.wait(); }

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(const ValueType& value,
                                             const int& thread_id) const {
    static_assert(std::is_trivially_default_constructible<ValueType>(),
                  "Only trivial constructible types can be broadcasted");
    tile_static ValueType local_value;
    zero_init(local_value);
    if (this->team_rank() == thread_id) {
      local_value = value;
    }
    this->team_barrier();
    value = local_value;
  }
  // Reduce across a team of threads.
  //
  // Each thread has vector_length elements.
  // This reduction is for TeamThreadRange operations, where the range
  // is spread across threads.  Effectively, there are vector_length
  // independent reduction operations.
  // This is different from a reduction across the elements of a thread,
  // which reduces every vector element.

  template <class ValueType, class JoinOp>
  KOKKOS_INLINE_FUNCTION ValueType team_reduce(const ValueType& value,
                                               const JoinOp& op_in) const {
    typedef JoinLambdaAdapter<ValueType, JoinOp> JoinOpFunctor;
    const JoinOpFunctor op(op_in);

    tile_static ValueType buffer[512];
    const auto local = lindex();
    const auto team  = team_rank();
    auto vector_rank = local % m_vector_length;
    auto thread_base = team * m_vector_length;

    const std::size_t size = next_pow_2(m_team_size + 1) / 2;
#if defined(ROCM15)
    buffer[local] = value;
#else
    // ROCM 1.5 handles address spaces better, previous version didn't
    lds_for(buffer[local], [&](ValueType& x) { x = value; });
#endif
    m_idx.barrier.wait();

    for (std::size_t s = 1; s < size; s *= 2) {
      const std::size_t index = 2 * s * team;
      if (index < size) {
#if defined(ROCM15)
        op.join(buffer[vector_rank + index * m_vector_length],
                buffer[vector_rank + (index + s) * m_vector_length]);
#else
        lds_for(buffer[vector_rank + index * m_vector_length],
                [&](ValueType& x) {
                  lds_for(buffer[vector_rank + (index + s) * m_vector_length],
                          [&](ValueType& y) { op.join(x, y); });
                });
#endif
      }
      m_idx.barrier.wait();
    }

    if (local == 0) {
      for (int i = size * m_vector_length; i < m_team_size * m_vector_length;
           i += m_vector_length)
#if defined(ROCM15)
        op.join(buffer[vector_rank], buffer[vector_rank + i]);
#else
        lds_for(buffer[vector_rank], [&](ValueType& x) {
          lds_for(buffer[vector_rank + i],
                  [&](ValueType& y) { op.join(x, y); });
        });
#endif
    }
    m_idx.barrier.wait();

    return buffer[0];
  }

  // Reduce across a team of threads, with a reducer data type
  //
  // Each thread has vector_length elements.
  // This reduction is for TeamThreadRange operations, where the range
  // is spread across threads.  Effectively, there are vector_length
  // independent reduction operations.
  // This is different from a reduction across the elements of a thread,
  // which reduces every vector element.

  template <class ReducerType>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      team_reduce(const ReducerType& reducer) const {
    typedef typename ReducerType::value_type value_type;

    tile_static value_type buffer[512];
    const auto local = lindex();
    const auto team  = team_rank();
    auto vector_rank = local % m_vector_length;
    auto thread_base = team * m_vector_length;

    const std::size_t size = next_pow_2(m_team_size + 1) / 2;
#if defined(ROCM15)
    buffer[local] = reducer.reference();
#else
    // ROCM 1.5 handles address spaces better, previous version didn't
    lds_for(buffer[local], [&](ValueType& x) { x = value; });
#endif
    m_idx.barrier.wait();

    for (std::size_t s = 1; s < size; s *= 2) {
      const std::size_t index = 2 * s * team;
      if (index < size) {
#if defined(ROCM15)
        reducer.join(buffer[vector_rank + index * m_vector_length],
                     buffer[vector_rank + (index + s) * m_vector_length]);
#else
        lds_for(buffer[vector_rank + index * m_vector_length],
                [&](ValueType& x) {
                  lds_for(buffer[vector_rank + (index + s) * m_vector_length],
                          [&](ValueType& y) { reducer.join(x, y); });
                });
#endif
      }
      m_idx.barrier.wait();
    }

    if (local == 0) {
      for (int i = size * m_vector_length; i < m_team_size * m_vector_length;
           i += m_vector_length)
#if defined(ROCM15)
        reducer.join(buffer[vector_rank], buffer[vector_rank + i]);
#else
        lds_for(buffer[vector_rank], [&](ValueType& x) {
          lds_for(buffer[vector_rank + i],
                  [&](ValueType& y) { reducer.join(x, y); });
        });
#endif
    }
    m_idx.barrier.wait();
    reducer.reference() = buffer[0];
  }

  /** \brief  Intra-team vector reduce
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The intra-team accumulation value will, at the end of the
   *  league's parallel execution, be the reduction's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's vector reduce operation is
   *  similarly non-deterministic.
   */
  template <class ValueType, class JoinOp>
  KOKKOS_INLINE_FUNCTION ValueType thread_reduce(const ValueType& value,
                                                 const JoinOp& op_in) const {
    typedef JoinLambdaAdapter<ValueType, JoinOp> JoinOpFunctor;
    const JoinOpFunctor op(op_in);

    const auto local = m_idx.local[0];
    tile_static ValueType buffer[512];
    const std::size_t size =
        m_vector_length;  // vector length must be power of 2
    auto vector_rank = local % m_vector_length;
    auto thread_base = team_rank() * m_vector_length;
    lds_for(buffer[local], [&](ValueType& x) { x = value; });
    m_idx.barrier.wait();
    for (std::size_t s = 1; s < size; s *= 2) {
      const std::size_t index = 2 * s * vector_rank;
      if (index < size) {
#if defined(ROCM15)
        op.join(buffer[thread_base + index], buffer[thread_base + index + s]);
#else

        lds_for(buffer[thread_base + index], [&](ValueType& x) {
          lds_for(buffer[thread_base + index + s],
                  [&](ValueType& y) { op.join(x, y); });
        });
#endif
      }
      m_idx.barrier.wait();
    }

    m_idx.barrier.wait();
    return buffer[thread_base];
  }

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      vector_reduce(ReducerType const& reducer) const {
#ifdef __HCC_ACCELERATOR__
    if (m_vector_length == 1) return;

    // Intra vector lane shuffle reduction:
    typename ReducerType::value_type tmp(reducer.reference());

    for (int i = m_vector_length; (i >>= 1);) {
      reducer.reference() = shfl_down(tmp, i, m_vector_length);
      if ((int)vector_rank() < i) {
        reducer.join(tmp, reducer.reference());
      }
    }

    // Broadcast from root lane to all other lanes.
    // Cannot use "butterfly" algorithm to avoid the broadcast
    // because floating point summation is not associative
    // and thus different threads could have different results.

    reducer.reference() = shfl(tmp, 0, m_vector_length);
#endif
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
  template <typename Type>
  KOKKOS_INLINE_FUNCTION Type
  team_scan(const Type& value, Type* const global_accum = nullptr) const {
#if 0
      const auto local = m_idx.local[0];
      const auto last = m_team_size - 1;
      const auto init = 0;
      tile_static Type buffer[256];

      if (local == last) buffer[0] = init;
      else buffer[local] = value;

      m_idx.barrier.wait();

      for(std::size_t s = 1; s < m_team_size; s *= 2)
      {
          if (local >= s) buffer[local] += buffer[local - s];
          m_idx.barrier.wait();
      }

      if ( global_accum )
      { 
         if(local == last)
         {
            atomic_fetch_add(global_accum, buffer[local] + value);
         }
         m_idx.barrier.wait();
         buffer[local] += *global_accum;
      }
      m_idx.barrier.wait();
      return buffer[local];
#else
    tile_static Type sarray[2][256 + 1];
    int lid = m_idx.local[0];
    int lp1 = lid + 1;

    int toggle  = 1;
    int _toggle = 0;
    m_idx.barrier.wait();

    if (lid == 0) {
      sarray[1][0] = 0;
      sarray[0][0] = 0;
    }
    sarray[1][lp1] = value;

    m_idx.barrier.wait();
    for (int stride = 1; stride < m_team_size; stride *= 2) {
      if (lid >= stride) {
        sarray[_toggle][lp1] =
            sarray[toggle][lp1] + sarray[toggle][lp1 - stride];
      } else {
        sarray[_toggle][lp1] = sarray[toggle][lp1];
      }
      toggle  = _toggle;
      _toggle = 1 - toggle;
      m_idx.barrier.wait();
    }

    if (global_accum) {
      if (m_team_size == lp1) {
        sarray[toggle][m_team_size] =
            atomic_fetch_add(global_accum, sarray[toggle][m_team_size]);
      }
      m_idx.barrier.wait();
      sarray[toggle][lid] += sarray[toggle][m_team_size];
    }
    m_idx.barrier.wait();
    return sarray[toggle][lid];
#endif
  }

 private:
  int m_league_size;
  int m_team_size;
  const scratch_memory_space m_team_shared;

 public:
  int m_vector_length;
  hc::tiled_index<1> m_idx;
};
}  // namespace Impl
}  // namespace Kokkos
#include <ROCm/Kokkos_ROCm_ReduceScan.hpp>

namespace Kokkos {
namespace Impl {

//----------------------------------------------------------------------------

template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::RangePolicy<Traits...>,
                  Kokkos::Experimental::ROCm> {
 private:
  typedef Kokkos::RangePolicy<Traits...> Policy;

 public:
  inline ParallelFor(const FunctorType& f, const Policy& policy) {
    const auto len    = policy.end() - policy.begin();
    const auto offset = policy.begin();
    if (len == 0) return;
    // define a lambda to work around a compiler issue.  The compiler does not
    // properly dereference f inside the pfe.
    auto foo = [=](size_t i) { rocm_invoke<typename Policy::work_tag>(f, i); };

#if __hcc_workweek__ > 16600
    hc::parallel_for_each(
        hc::extent<1>(len),
        [=](const hc::index<1>& idx) [[hc]] [[hc_max_workgroup_dim(1024, 1, 1)]]
#else
    hc::parallel_for_each(
        hc::extent<1>(len).tile(256),
        [=](const hc::index<1>& idx) [[hc]]
#endif
        {
          if (idx[0] < len)  // workaround for Carrizo (and Fiji?)
            foo(idx[0] + offset);
        })
        .wait();
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}
};

// MDRangePolicy impl
template <class FunctorType, class... Traits>
class ParallelFor<FunctorType, Kokkos::MDRangePolicy<Traits...>,
                  Kokkos::Experimental::ROCm> {
 private:
  typedef Kokkos::MDRangePolicy<Traits...> Policy;
  using RP = Policy;
  typedef typename Policy::array_index_type array_index_type;
  typedef typename Policy::index_type index_type;
  typedef typename Policy::launch_bounds LaunchBounds;

  const FunctorType m_functor;
  const Policy m_rp;

 public:
  KOKKOS_INLINE_FUNCTION
  void operator()(void) const {
    Kokkos::Impl::Refactor::DeviceIterateTile<Policy::rank, Policy, FunctorType,
                                              typename Policy::work_tag>(
        m_rp, m_functor)
        .exec_range();
  }

  inline void execute() const {
    const array_index_type maxblocks = static_cast<array_index_type>(
        Kokkos::Impl::ROCmTraits::UpperBoundExtentCount);
    if (RP::rank == 2) {
      const dim3 block(m_rp.m_tile[0], m_rp.m_tile[1], 1);
      const dim3 grid(std::min(m_rp.m_upper[0] - m_rp.m_lower[0], maxblocks),
                      std::min(m_rp.m_upper[1] - m_rp.m_lower[1], maxblocks),
                      1);
      ROCmParallelLaunch<ParallelFor, LaunchBounds>(*this, grid, block, 0);
    } else if (RP::rank == 3) {
      const dim3 block(m_rp.m_tile[0], m_rp.m_tile[1], m_rp.m_tile[2]);
      const dim3 grid(std::min(m_rp.m_upper[0] - m_rp.m_lower[0], maxblocks),
                      std::min(m_rp.m_upper[1] - m_rp.m_lower[1], maxblocks),
                      std::min(m_rp.m_upper[2] - m_rp.m_lower[2], maxblocks));
      ROCmParallelLaunch<ParallelFor, LaunchBounds>(*this, grid, block, 0);
    } else if (RP::rank == 4) {
      // id0,id1 encoded within threadIdx.x; id2 to threadIdx.y; id3 to
      // threadIdx.z
      const dim3 block(m_rp.m_tile[0] * m_rp.m_tile[1], m_rp.m_tile[2],
                       m_rp.m_tile[3]);
      const dim3 grid(std::min(m_rp.m_tile_end[0] * m_rp.m_tile_end[1] *
                                   m_rp.m_tile[0] * m_rp.m_tile[1],
                               maxblocks),
                      std::min(m_rp.m_upper[2] - m_rp.m_lower[2], maxblocks),
                      std::min(m_rp.m_upper[3] - m_rp.m_lower[3], maxblocks));
      ROCmParallelLaunch<ParallelFor, LaunchBounds>(*this, grid, block, 0);
    } else if (RP::rank == 5) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4 to
      // threadIdx.z
      const dim3 block(m_rp.m_tile[0] * m_rp.m_tile[1],
                       m_rp.m_tile[2] * m_rp.m_tile[3], m_rp.m_tile[4]);
      const dim3 grid(std::min(m_rp.m_tile_end[0] * m_rp.m_tile_end[1] *
                                   m_rp.m_tile[0] * m_rp.m_tile[1],
                               maxblocks),
                      std::min(m_rp.m_tile_end[2] * m_rp.m_tile_end[3] *
                                   m_rp.m_tile[2] * m_rp.m_tile[3],
                               maxblocks),
                      std::min(m_rp.m_upper[4] - m_rp.m_lower[4], maxblocks));
      ROCmParallelLaunch<ParallelFor, LaunchBounds>(*this, grid, block, 0);
    } else if (RP::rank == 6) {
      // id0,id1 encoded within threadIdx.x; id2,id3 to threadIdx.y; id4,id5 to
      // threadIdx.z
      const dim3 block(m_rp.m_tile[0] * m_rp.m_tile[1],
                       m_rp.m_tile[2] * m_rp.m_tile[3],
                       m_rp.m_tile[4] * m_rp.m_tile[5]);
      const dim3 grid(std::min(m_rp.m_tile_end[0] * m_rp.m_tile_end[1] *
                                   m_rp.m_tile[0] * m_rp.m_tile[1],
                               maxblocks),
                      std::min(m_rp.m_tile_end[2] * m_rp.m_tile_end[3] *
                                   m_rp.m_tile[2] * m_rp.m_tile[3],
                               maxblocks),
                      std::min(m_rp.m_tile_end[4] * m_rp.m_tile_end[5] *
                                   m_rp.m_tile[4] * m_rp.m_tile[5],
                               maxblocks));
      ROCmParallelLaunch<ParallelFor, LaunchBounds>(*this, grid, block, 0);
    } else {
      printf("Kokkos::MDRange Error: Exceeded rank bounds with ROCm\n");
      Kokkos::abort("Aborting");
    }

  }  // end execute

  //  inline
  ParallelFor(const FunctorType& arg_functor, Policy arg_policy)
      : m_functor(arg_functor), m_rp(arg_policy) {}
};

//----------------------------------------------------------------------------

template <class F, class... Traits>
class ParallelFor<F, Kokkos::TeamPolicy<Traits...>,
                  Kokkos::Experimental::ROCm> {
  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::Experimental::ROCm, Traits...>;
  typedef Kokkos::Impl::FunctorValueTraits<F, typename Policy::work_tag>
      ValueTraits;

 public:
  inline ParallelFor(const F& f, const Policy& policy) {
    const auto league_size  = policy.league_size();
    const auto team_size    = policy.team_size();
    const int vector_length = policy.vector_length();
    const auto total_size   = league_size * team_size * vector_length;
    const int scratch_size0 = policy.scratch_size(0, team_size);
    const int scratch_size1 = policy.scratch_size(1, team_size);

    if (total_size == 0) return;

    const auto shared_size = FunctorTeamShmemSize<F>::value(f, team_size);
    char* scratch          = NULL;
    char* shared = (char*)rocm_device_allocate(shared_size * league_size +
                                               scratch_size0 * league_size);
    if (0 < scratch_size1)
      scratch = (char*)rocm_device_allocate(scratch_size1 * league_size);

    hc::extent<1> flat_extent(total_size);

    hc::tiled_extent<1> team_extent =
        flat_extent.tile(vector_length * team_size);
    hc::parallel_for_each(
        team_extent,
        [=](hc::tiled_index<1> idx) [[hc]] {
          rocm_invoke<typename Policy::work_tag>(
              f, typename Policy::member_type(
                     idx, league_size, team_size, shared, shared_size,
                     scratch_size0, scratch, scratch_size1, vector_length));
        })
        .wait();

    if (0 < scratch_size1) rocm_device_free(scratch);
    rocm_device_free(shared);
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}
};

//----------------------------------------------------------------------------

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::RangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::ROCm> {
 public:
  typedef Kokkos::RangePolicy<Traits...> Policy;

  // TODO: Use generic lambdas instead
  struct invoke_fn {
    template <class F, class... Ts>
    KOKKOS_INLINE_FUNCTION void operator()(std::size_t size, F&& f,
                                           hc::tiled_index<1> idx, tile_desc td,
                                           Ts&&... xs) const {
      auto global = idx.global[0];
      if (global < size) f(idx.global[0], static_cast<Ts&&>(xs)...);
    }
  };

  template <class ViewType>
  inline ParallelReduce(
      const FunctorType& f, const Policy& policy, const ViewType& result_view,
      typename std::enable_if<Kokkos::is_view<ViewType>::value &&
                                  !Kokkos::is_reducer_type<ReducerType>::value,
                              void*>::type = NULL) {
    typedef typename Policy::work_tag Tag;
    typedef Kokkos::Impl::FunctorValueTraits<FunctorType, Tag> ValueTraits;
    typedef Kokkos::Impl::FunctorValueInit<FunctorType, Tag> ValueInit;
    typedef typename ValueTraits::reference_type reference_type;

    const auto total_size = policy.end() - policy.begin();

    if (total_size == 0) {
      if (result_view.data()) {
        ValueInit::init(f, result_view.data());
      }
      return;
    }

    Kokkos::Impl::reduce_enqueue<Tag>(
        total_size, f, InvalidType{}, rocm_capture(invoke_fn{}, total_size),
        result_view.data(), result_view.extent(0));
  }

  inline ParallelReduce(const FunctorType& f, Policy policy,
                        const ReducerType& reducer) {
    typedef typename Policy::work_tag Tag;

    typedef Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                               FunctorType, ReducerType>
        ReducerConditional;
    typedef typename ReducerConditional::type ReducerTypeFwd;
    typedef Kokkos::Impl::FunctorValueTraits<FunctorType, Tag> ValueTraits;
    typedef Kokkos::Impl::FunctorValueInit<ReducerType, Tag> ValueInit;

    typedef typename ValueTraits::reference_type reference_type;

    const auto total_size = policy.end() - policy.begin();

    if (total_size == 0) {
      if (reducer.view().data()) {
        ValueInit::init(ReducerConditional::select(f, reducer),
                        reducer.view().data());
      }
      return;
    }

    Kokkos::Impl::reduce_enqueue<Tag>(
        total_size, f, reducer, rocm_capture(invoke_fn{}, total_size),
        reducer.view().data(), reducer.view().extent(0));
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}
};

//----------------------------------------------------------------------------

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::MDRangePolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::ROCm> {
 private:
  typedef Kokkos::MDRangePolicy<Traits...> Policy;
  using RP = Policy;
  typedef typename Policy::array_index_type array_index_type;
  typedef typename Policy::index_type index_type;
  typedef typename Policy::work_tag WorkTag;
  typedef typename Policy::member_type Member;
  typedef typename Policy::launch_bounds LaunchBounds;

  typedef Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                             FunctorType, ReducerType>
      ReducerConditional;
  typedef typename ReducerConditional::type ReducerTypeFwd;
  typedef
      typename Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                                  WorkTag, void>::type WorkTagFwd;

  typedef Kokkos::Impl::FunctorValueTraits<ReducerTypeFwd, WorkTagFwd>
      ValueTraits;
  typedef Kokkos::Impl::FunctorValueInit<ReducerTypeFwd, WorkTagFwd> ValueInit;
  typedef Kokkos::Impl::FunctorValueJoin<ReducerTypeFwd, WorkTagFwd> ValueJoin;

 public:
  typedef typename ValueTraits::pointer_type pointer_type;
  typedef typename ValueTraits::value_type value_type;
  typedef typename ValueTraits::reference_type reference_type;
  typedef FunctorType functor_type;
  typedef Kokkos::Experimental::ROCm::size_type size_type;

  // Algorithmic constraints: blockSize is a power of two AND blockDim.y ==
  // blockDim.z == 1

  const FunctorType m_functor;
  const Policy m_policy;  // used for workrange and nwork
  const ReducerType m_reducer;
  const pointer_type m_result_ptr;
  value_type* m_scratch_space;
  size_type* m_scratch_flags;

  typedef typename Kokkos::Impl::Reduce::DeviceIterateTile<
      Policy::rank, Policy, FunctorType, typename Policy::work_tag,
      reference_type>
      DeviceIteratePattern;

  KOKKOS_INLINE_FUNCTION
  void exec_range(reference_type update) const {
    Kokkos::Impl::Reduce::DeviceIterateTile<Policy::rank, Policy, FunctorType,
                                            typename Policy::work_tag,
                                            reference_type>(m_policy, m_functor,
                                                            update)
        .exec_range();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(void) const { run(); }

  KOKKOS_INLINE_FUNCTION
  void run() const {
    const integral_nonzero_constant<size_type, ValueTraits::StaticValueSize /
                                                   sizeof(value_type)>
        word_count((ValueTraits::value_size(
                       ReducerConditional::select(m_functor, m_reducer))) /
                   sizeof(value_type));
    // pointer to shared data accounts for the reserved space at the start
    value_type* const shared =
        kokkos_impl_rocm_shared_memory<value_type>() + 2 * sizeof(uint64_t);

    {
      reference_type value =
          ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                          shared + threadIdx_y * word_count.value);
      // Number of blocks is bounded so that the reduction can be limited to two
      // passes. Each thread block is given an approximately equal amount of
      // work to perform. Accumulate the values for this block. The accumulation
      // ordering does not match the final pass, but is arithmatically
      // equivalent.

      this->exec_range(value);
    }

    // Reduce with final value at blockDim.y - 1 location.
    // Problem: non power-of-two blockDim

    if (rocm_single_inter_block_reduce_scan<false, ReducerTypeFwd, WorkTagFwd>(
            ReducerConditional::select(m_functor, m_reducer), blockIdx_x,
            gridDim_x, shared, m_scratch_space, m_scratch_flags)) {
      // This is the final block with the final result at the final threads'
      // location
      value_type* const tshared = shared + (blockDim_y - 1) * word_count.value;
      value_type* const global  = m_scratch_space;

      if (threadIdx_y == 0) {
        Kokkos::Impl::FunctorFinal<ReducerTypeFwd, WorkTagFwd>::final(
            ReducerConditional::select(m_functor, m_reducer), tshared);
        //        for ( unsigned i = 0 ; i < word_count.value ; i+=blockDim_y )
        //        { global[i] = tshared[i]; }
        for (unsigned i = 0; i < word_count.value; i++) {
          global[i] = tshared[i];
        }
      }
    }
  }

  // Determine block size constrained by shared memory:
  static inline unsigned local_block_size(const FunctorType& f) {
    unsigned n = ROCmTraits::WavefrontSize * 8;
    while (n &&
           ROCmTraits::SharedMemoryCapacity <
               rocm_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                         WorkTag>(f, n)) {
      n >>= 1;
    }
    return n;
  }

  inline void execute() {
    const int nwork = m_policy.m_num_tiles;
    if (nwork) {
      int block_size = m_policy.m_prod_tile_dims;
      // CONSTRAINT: Algorithm requires block_size >= product of tile dimensions
      // Nearest power of two
      int exponent_pow_two = std::ceil(std::log2((float)block_size));
      block_size           = 1 << (exponent_pow_two);

      m_scratch_space = (value_type*)rocm_internal_scratch_space(
          ValueTraits::value_size(
              ReducerConditional::select(m_functor, m_reducer)) *
          block_size * nwork /* block_size == max block_count */);
      m_scratch_flags = rocm_internal_scratch_flags(sizeof(size_type));
      const dim3 block(1, block_size, 1);
      // Required grid.x <= block.y
      const dim3 grid(nwork, block_size, 1);
      const int shmem =
          rocm_single_inter_block_reduce_scan_shmem<false, FunctorType,
                                                    WorkTag>(m_functor,
                                                             block.y);

      ROCmParallelLaunch<ParallelReduce, LaunchBounds>(
          *this, grid, block, shmem);  // copy to device and execute

      ROCM().fence();

      if (m_result_ptr) {
        const int size = ValueTraits::value_size(
            ReducerConditional::select(m_functor, m_reducer));
        DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace>(
            m_result_ptr, m_scratch_space, size);
      }
    } else {
      if (m_result_ptr) {
        ValueInit::init(ReducerConditional::select(m_functor, m_reducer),
                        m_result_ptr);
      }
    }
  }

  template <class HostViewType>
  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const HostViewType& arg_result,
                 typename std::enable_if<Kokkos::is_view<HostViewType>::value,
                                         void*>::type = NULL)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(InvalidType()),
        m_result_ptr(arg_result.data()),
        m_scratch_space(0),
        m_scratch_flags(0) {}

  ParallelReduce(const FunctorType& arg_functor, const Policy& arg_policy,
                 const ReducerType& reducer)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_reducer(reducer),
        m_result_ptr(reducer.view().data()),
        m_scratch_space(0),
        m_scratch_flags(0) {}
};
//----------------------------------------------------------------------------

template <class FunctorType, class ReducerType, class... Traits>
class ParallelReduce<FunctorType, Kokkos::TeamPolicy<Traits...>, ReducerType,
                     Kokkos::Experimental::ROCm> {
  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::Experimental::ROCm, Traits...>;
  typedef Kokkos::Impl::FunctorValueTraits<FunctorType,
                                           typename Policy::work_tag>
      ValueTraits;

 public:
  struct invoke_fn {
    template <class Create, class F, class... Ts>
    KOKKOS_INLINE_FUNCTION void operator()(Create&& create, F&& f,
                                           hc::tiled_index<1> idx, tile_desc td,
                                           Ts&&... xs) const {
      f(create(idx, td), static_cast<Ts&&>(xs)...);
    }
  };

  template <class ViewType>
  inline ParallelReduce(
      const FunctorType& f, const Policy& policy, const ViewType& result_view,
      typename std::enable_if<Kokkos::is_view<ViewType>::value &&
                                  !Kokkos::is_reducer_type<ReducerType>::value,
                              void*>::type = NULL) {
    const int league_size   = policy.league_size();
    const int team_size     = policy.team_size(f);
    const int vector_length = policy.vector_length();
    const int scratch_size0 = policy.scratch_size(0, team_size);
    const int scratch_size1 = policy.scratch_size(1, team_size);
    const int total_size    = league_size * team_size;

    typedef Kokkos::Impl::FunctorValueInit<FunctorType,
                                           typename Policy::work_tag>
        ValueInit;
    if (total_size == 0) {
      if (result_view.data()) {
        ValueInit::init(f, result_view.data());
      }
      return;
    }

    const int reduce_size = ValueTraits::value_size(f);
    const int shared_size =
        FunctorTeamShmemSize<FunctorType>::value(f, team_size);

    char* shared;
    char* scratch = NULL;

    shared = (char*)rocm_device_allocate(league_size *
                                         (shared_size + scratch_size0));
    if (0 < scratch_size1)
      scratch = (char*)rocm_device_allocate(scratch_size1 * league_size);

    auto create_team_member = [=](hc::tiled_index<1> idx, tile_desc td) {
      return typename Policy::member_type(
          idx, league_size, td.team_size, shared, shared_size, scratch_size0,
          scratch, scratch_size1, vector_length);
    };

    Kokkos::Impl::reduce_enqueue<typename Policy::work_tag>(
        total_size * vector_length, f, InvalidType{},
        rocm_capture(invoke_fn{}, create_team_member),
        result_view.ptr_on_device(), result_view.dimension_0(), team_size,
        vector_length, shared_size);

    if (0 < scratch_size1) rocm_device_free(scratch);
    rocm_device_free(shared);
  }

  inline ParallelReduce(const FunctorType& f, Policy policy,
                        const ReducerType& reducer) {
    const int league_size   = policy.league_size();
    const int team_size     = policy.team_size(f);
    const int vector_length = policy.vector_length();
    const int total_size    = league_size * team_size;

    typedef Kokkos::Impl::FunctorValueInit<ReducerType,
                                           typename Policy::work_tag>
        ValueInit;
    typedef Kokkos::Impl::if_c<std::is_same<InvalidType, ReducerType>::value,
                               FunctorType, ReducerType>
        ReducerConditional;
    if (total_size == 0) {
      if (reducer.view().data()) {
        ValueInit::init(ReducerConditional::select(f, reducer),
                        reducer.view().data());
      }
      return;
    }

    const int reduce_size = ValueTraits::value_size(f);
    const int shared_size =
        FunctorTeamShmemSize<FunctorType>::value(f, team_size);
    const int scratch_size0 = policy.scratch_size(0, team_size);
    const int scratch_size1 = policy.scratch_size(1, team_size);

    char* shared;
    char* scratch = NULL;
    shared        = (char*)rocm_device_allocate((shared_size + scratch_size0) *
                                         league_size);
    if (0 < scratch_size1)
      scratch = (char*)rocm_device_allocate(scratch_size1 * league_size);

    auto create_team_member = [=](hc::tiled_index<1> idx, tile_desc td) {
      return typename Policy::member_type(
          idx, league_size, td.tile_size, shared, shared_size, scratch_size0,
          scratch, scratch_size1, vector_length);
    };

    Kokkos::Impl::reduce_enqueue<typename Policy::work_tag>(
        league_size, f, reducer, rocm_capture(invoke_fn{}, create_team_member),
        reducer.view().data(), reducer.view().extent(0), team_size,
        vector_length, shared_size);

    if (0 < scratch_size1) rocm_device_free(scratch);
    rocm_device_free(shared);
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}
};

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::RangePolicy<Traits...>,
                   Kokkos::Experimental::ROCm> {
 private:
  typedef Kokkos::RangePolicy<Traits...> Policy;
  typedef typename Policy::work_tag Tag;
  typedef Kokkos::Impl::FunctorValueTraits<FunctorType, Tag> ValueTraits;

 public:
  //----------------------------------------

  inline ParallelScan(const FunctorType& f, const Policy& policy) {
    const auto len = policy.end() - policy.begin();

    if (len == 0) return;

    scan_enqueue<Tag>(
        len, f, [](hc::tiled_index<1> idx, int, int) { return idx.global[0]; });
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}

  //----------------------------------------
};

template <class FunctorType, class ReturnType, class... Traits>
class ParallelScanWithTotal<FunctorType, Kokkos::RangePolicy<Traits...>,
                            ReturnType, Kokkos::Experimental::ROCm> {
 private:
  typedef Kokkos::RangePolicy<Traits...> Policy;
  typedef typename Policy::work_tag Tag;
  typedef Kokkos::Impl::FunctorValueTraits<FunctorType, Tag> ValueTraits;

 public:
  //----------------------------------------

  inline ParallelScanWithTotal(const FunctorType& f, const Policy& policy,
                               ReturnType& arg_returnvalue) {
    const auto len = policy.end() - policy.begin();

    if (len == 0) return;

    scan_enqueue<Tag, ReturnType>(
        len, f, arg_returnvalue,
        [](hc::tiled_index<1> idx, int, int) { return idx.global[0]; });
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}

  //----------------------------------------
};

template <class FunctorType, class... Traits>
class ParallelScan<FunctorType, Kokkos::TeamPolicy<Traits...>,
                   Kokkos::Experimental::ROCm> {
 private:
  using Policy =
      Kokkos::Impl::TeamPolicyInternal<Kokkos::Experimental::ROCm, Traits...>;
  typedef typename Policy::work_tag Tag;
  typedef Kokkos::Impl::FunctorValueTraits<FunctorType, Tag> ValueTraits;

 public:
  //----------------------------------------

  inline ParallelScan(const FunctorType& f, const Policy& policy) {
    const auto league_size = policy.league_size();
    const auto team_size   = policy.team_size(f);
    const auto len         = league_size * team_size;

    if (len == 0) return;

    scan_enqueue<Tag>(
        len, f, [&](hc::tiled_index<1> idx, int n_teams, int n_leagues) {
          return typename Policy::member_type(idx, n_leagues, n_teams);
        });
  }

  KOKKOS_INLINE_FUNCTION
  void execute() const {}

  //----------------------------------------
};

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {
template <typename iType>
struct TeamThreadRangeBoundariesStruct<iType, ROCmTeamMember> {
  typedef iType index_type;
  const iType start;
  const iType end;
  const iType increment;
  const ROCmTeamMember& thread;

#if defined(__HCC_ACCELERATOR__)
  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const ROCmTeamMember& thread_,
                                  const iType& count)
      : start(thread_.team_rank()),
        end(count),
        increment(thread_.team_size()),
        thread(thread_) {}
  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const ROCmTeamMember& thread_,
                                  const iType& begin_, const iType& end_)
      : start(begin_ + thread_.team_rank()),
        end(end_),
        increment(thread_.team_size()),
        thread(thread_) {}
#else
  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const ROCmTeamMember& thread_,
                                  const iType& count)
      : start(0), end(count), increment(1), thread(thread_) {}
  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const ROCmTeamMember& thread_,
                                  const iType& begin_, const iType& end_)
      : start(begin_), end(end_), increment(1), thread(thread_) {}
#endif
};

template <typename iType>
struct ThreadVectorRangeBoundariesStruct<iType, ROCmTeamMember> {
  typedef iType index_type;
  const index_type start;
  const index_type end;
  const index_type increment;
  const ROCmTeamMember& thread;

#if defined(__HCC_ACCELERATOR__)
  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const ROCmTeamMember& thread_,
                                    const index_type& count)
      : start(thread_.lindex() % thread_.vector_length()),
        end(count),
        increment(thread_.vector_length()),
        thread(thread_) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const ROCmTeamMember& thread_,
                                    const index_type& arg_begin,
                                    const index_type& arg_end)
      : start(arg_begin + thread_.lindex() % thread_.vector_length()),
        end(arg_end),
        increment(thread_.vector_length()),
        thread(thread_) {}

//    KOKKOS_INLINE_FUNCTION
//    ThreadVectorRangeBoundariesStruct (const index_type& count):
//      start( 0 ),
//      end( count ),
//      increment( 1 )
//    {}
#else
  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const ROCmTeamMember& thread_,
                                    const index_type& count)
      : start(static_cast<index_type>(0)),
        end(count),
        increment(static_cast<index_type>(1)),
        thread(thread_) {}
  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const index_type& count)
      : start(static_cast<index_type>(0)),
        end(count),
        increment(static_cast<index_type>(1)) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const ROCmTeamMember& thread_,
                                    const index_type& arg_begin,
                                    const index_type& arg_end)
      : start(arg_begin),
        end(arg_end),
        increment(static_cast<index_type>(1)),
        thread(thread_) {}
  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const index_type& arg_begin,
                                    const index_type& arg_end)
      : start(arg_begin), end(arg_end), increment(static_cast<index_type>(1)) {}
#endif
};

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamThreadRangeBoundariesStruct<iType, Impl::ROCmTeamMember>
    TeamThreadRange(const Impl::ROCmTeamMember& thread, iType count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::ROCmTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, Impl::ROCmTeamMember>
TeamThreadRange(const Impl::ROCmTeamMember& thread, iType1 begin, iType2 end) {
  typedef typename std::common_type<iType1, iType2>::type iType;
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::ROCmTeamMember>(
      thread, begin, end);
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::ROCmTeamMember>
    ThreadVectorRange(const Impl::ROCmTeamMember& thread, iType count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::ROCmTeamMember>(
      thread, count);
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::ROCmTeamMember>
    ThreadVectorRange(const Impl::ROCmTeamMember& thread, iType arg_begin,
                      iType arg_end) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::ROCmTeamMember>(
      thread, arg_begin, arg_end);
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::ROCmTeamMember> PerTeam(
    const Impl::ROCmTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::ROCmTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::ROCmTeamMember> PerThread(
    const Impl::ROCmTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::ROCmTeamMember>(thread);
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::ROCmTeamMember>& single_struct,
    const FunctorType& lambda) {
  if (single_struct.team_member.vector_rank() == 0) lambda();
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::ROCmTeamMember>& single_struct,
    const FunctorType& lambda) {
  if ((single_struct.team_member.lindex() == 0)) lambda();
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::ROCmTeamMember>& single_struct,
    const FunctorType& lambda, ValueType& val) {
#if defined(ROCM15)
  // 1.5 needs this more proper restriction on which work units run
  if (single_struct.team_member.vector_rank() == 0) lambda(val);
  val = shfl(val, 0, single_struct.team_member.vector_length());
#else
  // but older compilers are fine with this (TestTeamVector::Test<
  // Kokkos::Experimental::ROCm >(4))
  lambda(val);
#endif
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::ROCmTeamMember>& single_struct,
    const FunctorType& lambda, ValueType& val) {
  if (single_struct.team_member.lindex() == 0) lambda(val);
  single_struct.team_member.team_broadcast(val, 0);
}

}  // namespace Kokkos

namespace Kokkos {

/** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team.
 * This functionality requires C++11 support.*/
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::ROCmTeamMember>&
        loop_boundaries,
    const Lambda& lambda) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment)
    lambda(i);
}

/** \brief  Inter-thread thread range parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team
 * and a summation of val is performed and put into result. This functionality
 * requires C++11 support.*/
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!Kokkos::is_reducer<ValueType>::value>::type
    parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                        iType, Impl::ROCmTeamMember>& loop_boundaries,
                    const Lambda& lambda, ValueType& result) {
  Kokkos::Sum<ValueType> reducer(result);
  reducer.init(reducer.reference());

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, reducer.reference());
  }
  loop_boundaries.thread.team_reduce(reducer);
}

/** \brief  Inter-thread thread range parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team
 * and a summation of val is performed and put into result. This functionality
 * requires C++11 support.*/
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<Kokkos::is_reducer<ReducerType>::value>::type
    parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                        iType, Impl::ROCmTeamMember>& loop_boundaries,
                    const Lambda& lambda, ReducerType const& reducer) {
  reducer.init(reducer.reference());

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, reducer.reference());
  }
  loop_boundaries.thread.team_reduce(reducer);
}

/** \brief  Intra-thread thread range parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a reduction of val is performed using JoinType(ValueType& val, const
 * ValueType& update) and put into init_result. The input value of init_result
 * is used as initializer for temporary variables of ValueType. Therefore the
 * input value should be the neutral element with respect to the join operation
 * (e.g. '0 for +-' or '1 for *'). This functionality requires C++11 support.*/
template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::ROCmTeamMember>&
        loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& result) {
#if defined(ROCM15)
  ValueType tmp = result;
  //  Simpler code works with ROCM1.5
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, tmp);
  }
  result = loop_boundaries.thread.team_reduce(tmp, join);
#else
  // this workaround freezes up with ROCM1.5, but needed for earlier compilers
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i, tmp);
    join(result, tmp);
  }
  result = loop_boundaries.thread.team_reduce(result, join);
#endif
  //  Impl::rocm_intra_workgroup_reduction( loop_boundaries.thread,
  //  result,join); Impl::rocm_inter_workgroup_reduction(
  //  loop_boundaries.thread, result,join);
}

}  // namespace Kokkos

namespace Kokkos {
/** \brief  Intra-thread vector parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread.
 * This functionality requires C++11 support.*/
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::ROCmTeamMember>&
        loop_boundaries,
    const Lambda& lambda) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment)
    lambda(i);
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a summation of val is performed and put into result. This functionality
 * requires C++11 support.*/
template <typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!Kokkos::is_reducer<ValueType>::value>::type
    parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<
                        iType, Impl::ROCmTeamMember>& loop_boundaries,
                    const Lambda& lambda, ValueType& result) {
  result = ValueType();

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    ValueType tmp = ValueType();
    lambda(i, tmp);
    result += tmp;
  }
  result =
      loop_boundaries.thread.thread_reduce(result, Impl::JoinAdd<ValueType>());
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a reduction of val is performed using JoinType(ValueType& val, const
 * ValueType& update) and put into init_result. The input value of init_result
 * is used as initializer for temporary variables of ValueType. Therefore the
 * input value should be the neutral element with respect to the join operation
 * (e.g. '0 for +-' or '1 for *'). This functionality requires C++11 support.*/
template <typename iType, class Lambda, typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::ROCmTeamMember>&
        loop_boundaries,
    const Lambda& lambda, const JoinType& join, ValueType& result) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
    loop_boundaries.thread.team_barrier();
  }
  result = loop_boundaries.thread.thread_reduce(result, join);
}

/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a summation of val is performed and put into result. This functionality
 * requires C++11 support.*/
template <typename iType, class Lambda, typename ReducerType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<Kokkos::is_reducer<ReducerType>::value>::type
    parallel_reduce(const Impl::ThreadVectorRangeBoundariesStruct<
                        iType, Impl::ROCmTeamMember>& loop_boundaries,
                    const Lambda& lambda, ReducerType const& reducer) {
  reducer.init(reducer.reference());

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, reducer.reference());
  }
  loop_boundaries.thread.vector_reduce(reducer);
}
/** \brief  Intra-thread vector parallel_reduce. Executes lambda(iType i,
 * ValueType & val) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all vector lanes of the the calling thread
 * and a reduction of val is performed using JoinType(ValueType& val, const
 * ValueType& update) and put into init_result. The input value of init_result
 * is used as initializer for temporary variables of ValueType. Therefore the
 * input value should be the neutral element with respect to the join operation
 * (e.g. '0 for +-' or '1 for *'). This functionality requires C++11 support.*/
template <typename iType, class Lambda, typename ReducerType, class JoinType>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::ROCmTeamMember>&
        loop_boundaries,
    const Lambda& lambda, const JoinType& join, ReducerType const& reducer) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, reducer.reference());
    loop_boundaries.thread.team_barrier();
  }
  reducer.reference() =
      loop_boundaries.thread.thread_reduce(reducer.reference(), join);
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
 * lanes. This functionality requires C++11 support.*/
template <typename iType, class FunctorType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::ROCmTeamMember>&
        loop_boundaries,
    const FunctorType& lambda) {
  typedef Kokkos::Impl::FunctorValueTraits<FunctorType, void> ValueTraits;
  typedef typename ValueTraits::value_type value_type;

  value_type val          = value_type();
  const int vector_length = loop_boundaries.thread.vector_length();
  const int vector_rank   = loop_boundaries.thread.vector_rank();

  iType end = ((loop_boundaries.end + vector_length - 1) / vector_length) *
              vector_length;
  value_type accum = value_type();

  for (int i = vector_rank; i < end; i += vector_length) {
    value_type val = 0;

    // First acquire per-lane contributions:
    if (i < loop_boundaries.end) lambda(i, val, false);

    value_type sval = val;

    // Bottom up inclusive scan in triangular pattern
    // where each thread is the root of a reduction tree
    // from the zeroth "lane" to itself.
    //  [t] += [t-1] if t >= 1
    //  [t] += [t-2] if t >= 2
    //  [t] += [t-4] if t >= 4
    //  ...

    for (int j = 1; j < vector_length; j <<= 1) {
      value_type tmp = 0;
      tmp            = shfl_up(sval, j, vector_length);
      if (j <= vector_rank) {
        sval += tmp;
      }
    }

    // Include accumulation and remove value for exclusive scan:
    val = accum + sval - val;

    // Provide exclusive scan value:
    if (i < loop_boundaries.end) lambda(i, val, true);

    // Accumulate the last value in the inclusive scan:
    sval = shfl(sval, vector_length - 1, vector_length);
    accum += sval;
  }
}

}  // namespace Kokkos
