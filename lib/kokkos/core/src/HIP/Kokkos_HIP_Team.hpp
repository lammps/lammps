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

#ifndef KOKKOS_HIP_TEAM_HPP
#define KOKKOS_HIP_TEAM_HPP

#include <Kokkos_Macros.hpp>

#if defined(__HIPCC__)

#include <utility>
#include <Kokkos_Parallel.hpp>

#include <HIP/Kokkos_HIP_KernelLaunch.hpp>
#include <HIP/Kokkos_HIP_ReduceScan.hpp>
#include <HIP/Kokkos_HIP_BlockSize_Deduction.hpp>
#include <Kokkos_Vectorization.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <typename Type>
struct HIPJoinFunctor {
  typedef Type value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type& update,
                   volatile const value_type& input) {
    update += input;
  }
};

/**\brief  Team member_type passed to TeamPolicy or TeamTask closures.
 *
 *  HIP thread blocks for team closures are dimensioned as:
 *    hipBlockDim_x == number of "vector lanes" per "thread"
 *    hipBlockDim_y == number of "threads" per team
 *    hipBlockDim_z == number of teams in a block
 *  where
 *    A set of teams exactly fill a warp OR a team is the whole block
 *      ( 0 == WarpSize % ( hipBlockDim_x * hipBlockDim_y ) )
 *      OR
 *      ( 1 == hipBlockDim_z )

 *  Thus when 1 < hipBlockDim_z the team is warp-synchronous
 *  and __syncthreads should not be called in team collectives.
 *
 *  When multiple teams are mapped onto a single block then the
 *  total available shared memory must be partitioned among teams.
 */
class HIPTeamMember {
 public:
  using execution_space      = Kokkos::Experimental::HIP;
  using scratch_memory_space = execution_space::scratch_memory_space;

 private:
  mutable void* m_team_reduce;
  scratch_memory_space m_team_shared;
  int m_team_reduce_size;
  int m_league_rank;
  int m_league_size;

 public:
  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_shmem() const {
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

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size; }
  KOKKOS_INLINE_FUNCTION int team_rank() const {
#ifdef __HIP_DEVICE_COMPILE__
    return hipThreadIdx_y;
#else
    return 0;
#endif
  }

  KOKKOS_INLINE_FUNCTION int team_size() const {
#ifdef __HIP_DEVICE_COMPILE__
    return hipBlockDim_y;
#else
    return 1;
#endif
  }

  KOKKOS_INLINE_FUNCTION void team_barrier() const {
#ifdef __HIP_DEVICE_COMPILE__
    if (1 == hipBlockDim_z)
      __syncthreads();  // team == block
    else
      __threadfence_block();  // team <= warp
#endif
  }

  //--------------------------------------------------------------------------

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(ValueType& val,
                                             const int& thread_id) const {
#ifdef __HIP_DEVICE_COMPILE__
    if (1 == hipBlockDim_z) {  // team == block
      __syncthreads();
      // Wait for shared data write until all threads arrive here
      if (hipThreadIdx_x == 0u &&
          hipThreadIdx_y == static_cast<uint32_t>(thread_id)) {
        *(reinterpret_cast<ValueType*>(m_team_reduce)) = val;
      }
      __syncthreads();  // Wait for shared data read until root thread writes
      val = *(reinterpret_cast<ValueType*>(m_team_reduce));
    } else {               // team <= warp
      ValueType tmp(val);  // input might not be a register variable
      ::Kokkos::Experimental::Impl::in_place_shfl(
          val, tmp, hipBlockDim_x * thread_id, hipBlockDim_x * hipBlockDim_y);
    }
#else
    (void)val;
    (void)thread_id;
#endif
  }

  template <class Closure, class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(Closure const& f, ValueType& val,
                                             const int& thread_id) const {
#ifdef __HIP_DEVICE_COMPILE__
    f(val);

    if (1 == hipBlockDim_z) {  // team == block
      __syncthreads();
      // Wait for shared data write until all threads arrive here
      if (hipThreadIdx_x == 0u &&
          hipThreadIdx_y == static_cast<uint32_t>(thread_id)) {
        *(reinterpret_cast<ValueType*>(m_team_reduce)) = val;
      }
      __syncthreads();  // Wait for shared data read until root thread writes
      val = *(reinterpret_cast<ValueType*>(m_team_reduce));
    } else {               // team <= warp
      ValueType tmp(val);  // input might not be a register variable
      ::Kokkos::Experimental::Impl::in_place_shfl(
          val, tmp, hipBlockDim_x * thread_id, hipBlockDim_x * hipBlockDim_y);
    }
#else
    (void)f;
    (void)val;
    (void)thread_id;
#endif
  }

  //--------------------------------------------------------------------------
  /**\brief  Reduction across a team
   *
   *  Mapping of teams onto blocks:
   *    hipBlockDim_x  is "vector lanes"
   *    hipBlockDim_y  is team "threads"
   *    hipBlockDim_z  is number of teams per block
   *
   *  Requires:
   *    hipBlockDim_x is power two
   *    hipBlockDim_x <= HIPTraits::WarpSize
   *    ( 0 == HIPTraits::WarpSize % ( hipBlockDim_x * hipBlockDim_y )
   *      OR
   *    ( 1 == hipBlockDim_z )
   */
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
                  typename ReducerType::value_type& value) const noexcept {
#ifdef __HIP_DEVICE_COMPILE__
    hip_intra_block_reduction(reducer, value, hipBlockDim_y);
#else
    (void)reducer;
    (void)value;
#endif
  }

  //--------------------------------------------------------------------------
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
  KOKKOS_INLINE_FUNCTION Type team_scan(const Type& value,
                                        Type* const global_accum) const {
#ifdef __HIP_DEVICE_COMPILE__
    Type* const base_data = reinterpret_cast<Type*>(m_team_reduce);

    __syncthreads();  // Don't write in to shared data until all threads have
                      // entered this function

    if (0 == hipThreadIdx_y) {
      base_data[0] = 0;
    }

    base_data[hipThreadIdx_y + 1] = value;

    Impl::hip_intra_block_reduce_scan<true, Impl::HIPJoinFunctor<Type>, void>(
        Impl::HIPJoinFunctor<Type>(), base_data + 1);

    if (global_accum) {
      if (hipBlockDim_y == hipThreadIdx_y + 1) {
        base_data[hipBlockDim_y] =
            atomic_fetch_add(global_accum, base_data[hipBlockDim_y]);
      }
      __syncthreads();  // Wait for atomic
      base_data[hipThreadIdx_y] += base_data[hipBlockDim_y];
    }

    return base_data[hipThreadIdx_y];
#else
    (void)value;
    (void)global_accum;
    return Type();
#endif
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

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION static
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      vector_reduce(ReducerType const& reducer) {
    vector_reduce(reducer, reducer.reference());
  }

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION static
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      vector_reduce(ReducerType const& reducer,
                    typename ReducerType::value_type& value) {
#ifdef __HIP_DEVICE_COMPILE__
    if (hipBlockDim_x == 1) return;

    // Intra vector lane shuffle reduction:
    typename ReducerType::value_type tmp(value);
    typename ReducerType::value_type tmp2 = tmp;

    int constexpr warp_size = ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
    unsigned mask =
        hipBlockDim_x == warp_size
            ? 0xffffffff
            : ((1 << hipBlockDim_x) - 1)
                  << ((hipThreadIdx_y % (warp_size / hipBlockDim_x)) *
                      hipBlockDim_x);

    for (int i = hipBlockDim_x; (i >>= 1);) {
      ::Kokkos::Experimental::Impl::in_place_shfl_down(tmp2, tmp, i,
                                                       hipBlockDim_x, mask);
      if (static_cast<int>(hipThreadIdx_x) < i) {
        reducer.join(tmp, tmp2);
      }
    }

    // Broadcast from root lane to all other lanes.
    // Cannot use "butterfly" algorithm to avoid the broadcast
    // because floating point summation is not associative
    // and thus different threads could have different results.

    ::Kokkos::Experimental::Impl::in_place_shfl(tmp2, tmp, 0, hipBlockDim_x,
                                                mask);
    value               = tmp2;
    reducer.reference() = tmp2;
#else
    (void)reducer;
    (void)value;
#endif
  }

  //--------------------------------------------------------------------------
  /**\brief  Global reduction across all blocks
   *
   *  Return !0 if reducer contains the final value
   */
  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION static
      typename std::enable_if<is_reducer<ReducerType>::value, int>::type
      global_reduce(ReducerType const& reducer, int* const global_scratch_flags,
                    void* const global_scratch_space, void* const shmem,
                    int const shmem_size) {
#ifdef __HIP_COMPILE_DEVICE__

    typedef typename ReducerType::value_type value_type;
    typedef value_type volatile* pointer_type;

    // Number of shared memory entries for the reduction:
    const int nsh = shmem_size / sizeof(value_type);

    // Number of HIP threads in the block, rank within the block
    const int nid = hipBlockDim_x * hipBlockDim_y * hipBlockDim_z;
    const int tid =
        hipThreadIdx_x +
        hipBlockDim_x * (hipThreadIdx_y + hipBlockDim_y * hipThreadIdx_z);

    // Reduces within block using all available shared memory
    // Contributes if it is the root "vector lane"

    // wn == number of warps in the block
    // wx == which lane within the warp
    // wy == which warp within the block

    const int wn =
        (nid + HIPTraits::WarpIndexMask) >> HIPTraits::WarpIndexShift;
    const int wx = tid & HIPTraits::WarpIndexMask;
    const int wy = tid >> HIPTraits::WarpIndexShift;

    //------------------------
    {  // Intra warp shuffle reduction from contributing CUDA threads

      value_type tmp(reducer.reference());

      int constexpr warp_size =
          ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
      for (int i = warp_size; static_cast<int>(hipBlockDim_x) <= (i >>= 1);) {
        Impl::in_place_shfl_down(reducer.reference(), tmp, i, warp_size);

        // Root of each vector lane reduces "thread" contribution
        if (0 == hipThreadIdx_x && wx < i) {
          reducer.join(&tmp, reducer.data());
        }
      }

      // Reduce across warps using shared memory.
      // Number of warps may not be power of two.

      __syncthreads();  // Wait before shared data write

      // Number of shared memory entries for the reduction
      // is at most one per warp
      const int nentry = wn < nsh ? wn : nsh;

      if (0 == wx && wy < nentry) {
        // Root thread of warp 'wy' has warp's value to contribute
        (reinterpret_cast<value_type*>(shmem))[wy] = tmp;
      }

      __syncthreads();  // Wait for write to be visible to block

      // When more warps than shared entries
      // then warps must take turns joining their contribution
      // to the designated shared memory entry.
      for (int i = nentry; i < wn; i += nentry) {
        const int k = wy - i;

        if (0 == wx && i <= wy && k < nentry) {
          // Root thread of warp 'wy' has warp's value to contribute
          reducer.join((reinterpret_cast<value_type*>(shmem)) + k, &tmp);
        }

        __syncthreads();  // Wait for write to be visible to block
      }

      // One warp performs the inter-warp reduction:

      if (0 == wy) {
        // Start fan-in at power of two covering nentry

        for (int i = (1 << (32 - __clz(nentry - 1))); (i >>= 1);) {
          const int k = wx + i;
          if (wx < i && k < nentry) {
            reducer.join((reinterpret_cast<pointer_type>(shmem)) + wx,
                         (reinterpret_cast<pointer_type>(shmem)) + k);
            __threadfence_block();  // Wait for write to be visible to warp
          }
        }
      }
    }
    //------------------------
    {  // Write block's value to global_scratch_memory

      int last_block = 0;

      if (0 == wx) {
        reducer.copy((reinterpret_cast<pointer_type>(global_scratch_space)) +
                         hipBlockIdx_x * reducer.length(),
                     reducer.data());

        __threadfence();  // Wait until global write is visible.

        last_block = static_cast<int>(hipGridDim_x) ==
                     1 + Kokkos::atomic_fetch_add(global_scratch_flags, 1);

        // If last block then reset count
        if (last_block) *global_scratch_flags = 0;
      }

      // FIXME hip does not support __syncthreads_or so we need to do it by hand
      // last_block = __syncthreads_or(last_block);

      __shared__ int last_block_shared;
      if (last_block) last_block_shared = last_block;
      __threadfence_block();

      if (!last_block_shared) return 0;
    }
    //------------------------
    // Last block reads global_scratch_memory into shared memory.

    const int nentry = nid < hipGridDim_x
                           ? (nid < nsh ? nid : nsh)
                           : (hipGridDim_x < nsh ? hipGridDim_x : nsh);

    // nentry = min( nid , nsh , gridDim.x )

    // whole block reads global memory into shared memory:

    if (tid < nentry) {
      const int offset = tid * reducer.length();

      reducer.copy(
          (reinterpret_cast<pointer_type>(shmem)) + offset,
          (reinterpret_cast<pointer_type>(global_scratch_space)) + offset);

      for (int i = nentry + tid; i < static_cast<int>(hipGridDim_x);
           i += nentry) {
        reducer.join((reinterpret_cast<pointer_type>(shmem)) + offset,
                     (reinterpret_cast<pointer_type>(global_scratch_space)) +
                         i * reducer.length());
      }
    }

    __syncthreads();  // Wait for writes to be visible to block

    if (0 == wy) {
      // Iterate to reduce shared memory to single warp fan-in size

      int constexpr warp_size =
          ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
      const int nreduce = warp_size < nentry ? warp_size : nentry;

      if (wx < nreduce && nreduce < nentry) {
        for (int i = nreduce + wx; i < nentry; i += nreduce) {
          reducer.join(((pointer_type)shmem) + wx, ((pointer_type)shmem) + i);
        }
        __threadfence_block();  // Wait for writes to be visible to warp
      }

      // Start fan-in at power of two covering nentry

      for (int i = (1 << (warp_size - __clz(nreduce - 1))); (i >>= 1);) {
        const int k = wx + i;
        if (wx < i && k < nreduce) {
          reducer.join((reinterpret_cast<pointer_type>(shmem)) + wx,
                       (reinterpret_cast<pointer_type>(shmem)) + k);
          __threadfence_block();  // Wait for writes to be visible to warp
        }
      }

      if (0 == wx) {
        reducer.copy(reducer.data(), reinterpret_cast<pointer_type>(shmem));
        return 1;
      }
    }
    return 0;

#else
    (void)reducer;
    (void)global_scratch_flags;
    (void)shmem;
    (void)global_scratch_space;
    (void)shmem_size;
    return 0;
#endif
  }

  //----------------------------------------
  // Private for the driver

  KOKKOS_INLINE_FUNCTION
  HIPTeamMember(void* shared, const int shared_begin, const int shared_size,
                void* scratch_level_1_ptr, const int scratch_level_1_size,
                const int arg_league_rank, const int arg_league_size)
      : m_team_reduce(shared),
        m_team_shared(((char*)shared) + shared_begin, shared_size,
                      scratch_level_1_ptr, scratch_level_1_size),
        m_team_reduce_size(shared_begin),
        m_league_rank(arg_league_rank),
        m_league_size(arg_league_size) {}

 public:
  // Declare to avoid unused private member warnings which are trigger
  // when SFINAE excludes the member function which uses these variables
  // Making another class a friend also surpresses these warnings
  bool impl_avoid_sfinae_warning() const noexcept {
    return m_team_reduce_size > 0 && m_team_reduce != nullptr;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <typename iType>
struct TeamThreadRangeBoundariesStruct<iType, HIPTeamMember> {
  typedef iType index_type;
  const HIPTeamMember& member;
  const iType start;
  const iType end;

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const HIPTeamMember& thread_, iType count)
      : member(thread_), start(0), end(count) {}

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const HIPTeamMember& thread_, iType begin_,
                                  iType end_)
      : member(thread_), start(begin_), end(end_) {}
};

template <typename iType>
struct TeamVectorRangeBoundariesStruct<iType, HIPTeamMember> {
  typedef iType index_type;
  const HIPTeamMember& member;
  const iType start;
  const iType end;

  KOKKOS_INLINE_FUNCTION
  TeamVectorRangeBoundariesStruct(const HIPTeamMember& thread_,
                                  const iType& count)
      : member(thread_), start(0), end(count) {}

  KOKKOS_INLINE_FUNCTION
  TeamVectorRangeBoundariesStruct(const HIPTeamMember& thread_,
                                  const iType& begin_, const iType& end_)
      : member(thread_), start(begin_), end(end_) {}
};

template <typename iType>
struct ThreadVectorRangeBoundariesStruct<iType, HIPTeamMember> {
  typedef iType index_type;
  const index_type start;
  const index_type end;

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const HIPTeamMember, index_type count)
      : start(static_cast<index_type>(0)), end(count) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(index_type count)
      : start(static_cast<index_type>(0)), end(count) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const HIPTeamMember, index_type arg_begin,
                                    index_type arg_end)
      : start(arg_begin), end(arg_end) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(index_type arg_begin, index_type arg_end)
      : start(arg_begin), end(arg_end) {}
};

}  // namespace Impl

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HIPTeamMember>
    TeamThreadRange(const Impl::HIPTeamMember& thread, iType count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HIPTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, Impl::HIPTeamMember>
TeamThreadRange(const Impl::HIPTeamMember& thread, iType1 begin, iType2 end) {
  typedef typename std::common_type<iType1, iType2>::type iType;
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HIPTeamMember>(
      thread, iType(begin), iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>
    TeamVectorRange(const Impl::HIPTeamMember& thread, const iType& count) {
  return Impl::TeamVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamVectorRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, Impl::HIPTeamMember>
TeamVectorRange(const Impl::HIPTeamMember& thread, const iType1& begin,
                const iType2& end) {
  typedef typename std::common_type<iType1, iType2>::type iType;
  return Impl::TeamVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>(
      thread, iType(begin), iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>
    ThreadVectorRange(const Impl::HIPTeamMember& thread, iType count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>(
      thread, count);
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>
    ThreadVectorRange(const Impl::HIPTeamMember& thread, iType arg_begin,
                      iType arg_end) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>(
      thread, arg_begin, arg_end);
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::HIPTeamMember> PerTeam(
    const Impl::HIPTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::HIPTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::HIPTeamMember> PerThread(
    const Impl::HIPTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::HIPTeamMember>(thread);
}

//----------------------------------------------------------------------------

/** \brief  Inter-thread parallel_for.
 *
 *  Executes closure(iType i) for each i=[0..N).
 *
 * The range [0..N) is mapped to all threads of the the calling thread team.
 */
template <typename iType, class Closure>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::HIPTeamMember>&
        loop_boundaries,
    const Closure& closure) {
#ifdef __HIP_DEVICE_COMPILE__
  for (iType i = loop_boundaries.start + hipThreadIdx_y;
       i < loop_boundaries.end; i += hipBlockDim_y)
    closure(i);
#else
  (void)loop_boundaries;
  (void)closure;
#endif
}

//----------------------------------------------------------------------------

/** \brief  Inter-thread parallel_reduce with a reducer.
 *
 *  Executes closure(iType i, ValueType & val) for each i=[0..N)
 *
 *  The range [0..N) is mapped to all threads of the
 *  calling thread team and a summation of val is
 *  performed and put into result.
 */
template <typename iType, class Closure, class ReducerType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<Kokkos::is_reducer<ReducerType>::value>::type
    parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                        iType, Impl::HIPTeamMember>& loop_boundaries,
                    const Closure& closure, const ReducerType& reducer) {
#ifdef __HIP_DEVICE_COMPILE__
  typename ReducerType::value_type value;
  reducer.init(value);

  for (iType i = loop_boundaries.start + hipThreadIdx_y;
       i < loop_boundaries.end; i += hipBlockDim_y) {
    closure(i, value);
  }

  loop_boundaries.member.team_reduce(reducer, value);
#else
  (void)loop_boundaries;
  (void)closure;
  (void)reducer;
#endif
}

/** \brief  Inter-thread parallel_reduce assuming summation.
 *
 *  Executes closure(iType i, ValueType & val) for each i=[0..N)
 *
 *  The range [0..N) is mapped to all threads of the
 *  calling thread team and a summation of val is
 *  performed and put into result.
 */
template <typename iType, class Closure, typename ValueType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!Kokkos::is_reducer<ValueType>::value>::type
    parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                        iType, Impl::HIPTeamMember>& loop_boundaries,
                    const Closure& closure, ValueType& result) {
#ifdef __HIP_DEVICE_COMPILE__
  ValueType val;
  Kokkos::Sum<ValueType> reducer(val);

  reducer.init(reducer.reference());

  for (iType i = loop_boundaries.start + hipThreadIdx_y;
       i < loop_boundaries.end; i += hipBlockDim_y) {
    closure(i, val);
  }

  loop_boundaries.member.team_reduce(reducer, val);
  result = reducer.reference();
#else
  (void)loop_boundaries;
  (void)closure;
  (void)result;
#endif
}

template <typename iType, class Closure>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>&
        loop_boundaries,
    const Closure& closure) {
#ifdef __HIP_DEVICE_COMPILE__
  for (iType i = loop_boundaries.start + hipThreadIdx_y * hipBlockDim_x +
                 hipThreadIdx_x;
       i < loop_boundaries.end; i += hipBlockDim_y * hipBlockDim_x)
    closure(i);
#else
  (void)loop_boundaries;
  (void)closure;
#endif
}

template <typename iType, class Closure, class ReducerType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<Kokkos::is_reducer<ReducerType>::value>::type
    parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                        iType, Impl::HIPTeamMember>& loop_boundaries,
                    const Closure& closure, const ReducerType& reducer) {
#ifdef __HIP_DEVICE_COMPILE__
  typename ReducerType::value_type value;
  reducer.init(value);

  for (iType i = loop_boundaries.start + hipThreadIdx_y * hipBlockDim_x +
                 hipThreadIdx_x;
       i < loop_boundaries.end; i += hipBlockDim_y * hipBlockDim_x) {
    closure(i, value);
  }

  loop_boundaries.member.vector_reduce(reducer, value);
  loop_boundaries.member.team_reduce(reducer, value);
#else
  (void)loop_boundaries;
  (void)closure;
  (void)reducer;
#endif
}

template <typename iType, class Closure, typename ValueType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!Kokkos::is_reducer<ValueType>::value>::type
    parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                        iType, Impl::HIPTeamMember>& loop_boundaries,
                    const Closure& closure, ValueType& result) {
#ifdef __HIP_DEVICE_COMPILE__
  ValueType val;
  Kokkos::Sum<ValueType> reducer(val);

  reducer.init(reducer.reference());

  for (iType i = loop_boundaries.start + hipThreadIdx_y * hipBlockDim_x +
                 hipThreadIdx_x;
       i < loop_boundaries.end; i += hipBlockDim_y * hipBlockDim_x) {
    closure(i, val);
  }

  loop_boundaries.member.vector_reduce(reducer);
  loop_boundaries.member.team_reduce(reducer);
  result = reducer.reference();
#else
  (void)loop_boundaries;
  (void)closure;
  (void)result;
#endif
}

//----------------------------------------------------------------------------

/** \brief  Intra-thread vector parallel_for.
 *
 *  Executes closure(iType i) for each i=[0..N)
 *
 * The range [0..N) is mapped to all vector lanes of the the calling thread.
 */
template <typename iType, class Closure>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>&
        loop_boundaries,
    const Closure& closure) {
#ifdef __HIP_DEVICE_COMPILE__
  for (iType i = loop_boundaries.start + hipThreadIdx_x;
       i < loop_boundaries.end; i += hipBlockDim_x) {
    closure(i);
  }
#else
  (void)loop_boundaries;
  (void)closure;
#endif
}

//----------------------------------------------------------------------------

/** \brief  Intra-thread vector parallel_reduce.
 *
 *  Calls closure(iType i, ValueType & val) for each i=[0..N).
 *
 *  The range [0..N) is mapped to all vector lanes of
 *  the calling thread and a reduction of val is performed using +=
 *  and output into result.
 *
 *  The identity value for the += operator is assumed to be the default
 *  constructed value.
 */
template <typename iType, class Closure, class ReducerType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<is_reducer<ReducerType>::value>::type
    parallel_reduce(Impl::ThreadVectorRangeBoundariesStruct<
                        iType, Impl::HIPTeamMember> const& loop_boundaries,
                    Closure const& closure, ReducerType const& reducer) {
#ifdef __HIP_DEVICE_COMPILE__
  reducer.init(reducer.reference());

  for (iType i = loop_boundaries.start + hipThreadIdx_x;
       i < loop_boundaries.end; i += hipBlockDim_x) {
    closure(i, reducer.reference());
  }

  Impl::HIPTeamMember::vector_reduce(reducer);
#else
  (void)loop_boundaries;
  (void)closure;
  (void)reducer;
#endif
}

/** \brief  Intra-thread vector parallel_reduce.
 *
 *  Calls closure(iType i, ValueType & val) for each i=[0..N).
 *
 *  The range [0..N) is mapped to all vector lanes of
 *  the calling thread and a reduction of val is performed using +=
 *  and output into result.
 *
 *  The identity value for the += operator is assumed to be the default
 *  constructed value.
 */
template <typename iType, class Closure, typename ValueType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!is_reducer<ValueType>::value>::type
    parallel_reduce(Impl::ThreadVectorRangeBoundariesStruct<
                        iType, Impl::HIPTeamMember> const& loop_boundaries,
                    Closure const& closure, ValueType& result) {
#ifdef __HIP_DEVICE_COMPILE__
  result = ValueType();

  for (iType i = loop_boundaries.start + hipThreadIdx_x;
       i < loop_boundaries.end; i += hipBlockDim_x) {
    closure(i, result);
  }

  Impl::HIPTeamMember::vector_reduce(Kokkos::Sum<ValueType>(result));
#else
  (void)loop_boundaries;
  (void)closure;
  (void)result;
#endif
}

//----------------------------------------------------------------------------

/** \brief  Intra-thread vector parallel exclusive prefix sum.
 *
 *  Executes closure(iType i, ValueType & val, bool final) for each i=[0..N)
 *
 *  The range [0..N) is mapped to all vector lanes in the
 *  thread and a scan operation is performed.
 *  The last call to closure has final == true.
 */
template <typename iType, class Closure>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::HIPTeamMember>&
        loop_boundaries,
    const Closure& closure) {
#ifdef __HIP_DEVICE_COMPILE__
  // Extract value_type from closure

  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, Closure>::value_type;

  // Loop through boundaries by vector-length chunks
  // must scan at each iteration

  value_type accum = 0;

  // All thread "lanes" must loop the same number of times.
  // Determine an loop end for all thread "lanes."
  // Requires:
  //   hipBlockDim_x is power of two and thus
  //     ( end % hipBlockDim_x ) == ( end & ( hipBlockDim_x - 1 ) )
  //   1 <= hipBlockDim_x <= HIPTraits::WarpSize

  int constexpr warp_size = ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
  const int mask          = hipBlockDim_x - 1;
  const unsigned active_mask =
      blockDim.x == warp_size
          ? 0xffffffff
          : ((1 << hipBlockDim_x) - 1)
                << (hipThreadIdx_y % (warp_size / hipBlockDim_x)) *
                       hipBlockDim_x;
  const int rem = loop_boundaries.end & mask;  // == end % hipBlockDim_x
  const int end = loop_boundaries.end + (rem ? hipBlockDim_x - rem : 0);

  for (int i = hipThreadIdx_x; i < end; i += hipBlockDim_x) {
    value_type val = 0;

    // First acquire per-lane contributions:
    if (i < loop_boundaries.end) closure(i, val, false);

    value_type sval = val;

    // Bottom up inclusive scan in triangular pattern
    // where each HIP thread is the root of a reduction tree
    // from the zeroth "lane" to itself.
    //  [t] += [t-1] if t >= 1
    //  [t] += [t-2] if t >= 2
    //  [t] += [t-4] if t >= 4
    //  ...

    for (int j = 1; j < static_cast<int>(hipBlockDim_x); j <<= 1) {
      value_type tmp = 0;
      ::Kokkos::Experimental::Impl::in_place_shfl_up(
          tmp, sval, j, hipBlockDim_x, active_mask);
      if (j <= static_cast<int>(hipThreadIdx_x)) {
        sval += tmp;
      }
    }

    // Include accumulation and remove value for exclusive scan:
    val = accum + sval - val;

    // Provide exclusive scan value:
    if (i < loop_boundaries.end) closure(i, val, true);

    // Accumulate the last value in the inclusive scan:
    ::Kokkos::Experimental::Impl::in_place_shfl(sval, sval, mask, blockDim.x,
                                                active_mask);

    accum += sval;
  }
#else
  (void)loop_boundaries;
  (void)closure;
#endif
}

}  // namespace Kokkos

namespace Kokkos {

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::HIPTeamMember>&,
    const FunctorType& lambda) {
#ifdef __HIP_DEVICE_COMPILE__
  if (hipThreadIdx_x == 0) lambda();
#else
  (void)lambda;
#endif
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::HIPTeamMember>&,
    const FunctorType& lambda) {
#ifdef __HIP_DEVICE_COMPILE__
  if (hipThreadIdx_x == 0 && hipThreadIdx_y == 0) lambda();
#else
  (void)lambda;
#endif
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::HIPTeamMember>&,
    const FunctorType& lambda, ValueType& val) {
#ifdef __HIP_DEVICE_COMPILE__
  int constexpr warp_size = ::Kokkos::Experimental::Impl::HIPTraits::WarpSize;
  if (hipThreadIdx_x == 0) lambda(val);
  unsigned mask = hipBlockDim_x == warp_size
                      ? 0xffffffff
                      : ((1 << hipBlockDim_x) - 1)
                            << ((hipThreadIdx_y % (warp_size / hipBlockDim_x)) *
                                hipBlockDim_x);
  ::Kokkos::Experimental::Impl::in_place_shfl(val, val, 0, hipBlockDim_x, mask);
#else
  (void)lambda;
  (void)val;
#endif
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::HIPTeamMember>& single_struct,
    const FunctorType& lambda, ValueType& val) {
  (void)single_struct;
  (void)lambda;
  (void)val;
#ifdef __HIP_DEVICE_COMPILE__
  if (hipThreadIdx_x == 0 && hipThreadIdx_y == 0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val, 0);
#else
  (void)single_struct;
  (void)lambda;
  (void)val;
#endif
}

}  // namespace Kokkos

#endif /* defined( __HIPCC__ ) */

#endif /* #ifndef KOKKOS_CUDA_TEAM_HPP */
