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

#ifndef KOKKOS_CUDA_TEAM_HPP
#define KOKKOS_CUDA_TEAM_HPP

#include <algorithm>

#include <Kokkos_Macros.hpp>

/* only compile this file if CUDA is enabled for Kokkos */
#if defined(KOKKOS_ENABLE_CUDA)

#include <utility>
#include <Kokkos_Parallel.hpp>

#include <Cuda/Kokkos_Cuda_KernelLaunch.hpp>
#include <Cuda/Kokkos_Cuda_ReduceScan.hpp>
#include <Cuda/Kokkos_Cuda_BlockSize_Deduction.hpp>
#include <Kokkos_Vectorization.hpp>

#include <impl/Kokkos_Tools.hpp>
#include <typeinfo>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <typename Type>
struct CudaJoinFunctor {
  using value_type = Type;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type& update, const value_type& input) {
    update += input;
  }
};

/**\brief  Team member_type passed to TeamPolicy or TeamTask closures.
 *
 *  Cuda thread blocks for team closures are dimensioned as:
 *    blockDim.x == number of "vector lanes" per "thread"
 *    blockDim.y == number of "threads" per team
 *    blockDim.z == number of teams in a block
 *  where
 *    A set of teams exactly fill a warp OR a team is the whole block
 *      ( 0 == WarpSize % ( blockDim.x * blockDim.y ) )
 *      OR
 *      ( 1 == blockDim.z )
 *
 *  Thus when 1 < blockDim.z the team is warp-synchronous
 *  and __syncthreads should not be called in team collectives.
 *
 *  When multiple teams are mapped onto a single block then the
 *  total available shared memory must be partitioned among teams.
 */
class CudaTeamMember {
 public:
  using execution_space      = Kokkos::Cuda;
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
    KOKKOS_IF_ON_DEVICE((return threadIdx.y;))
    KOKKOS_IF_ON_HOST((return 0;))
  }

  KOKKOS_INLINE_FUNCTION int team_size() const {
    KOKKOS_IF_ON_DEVICE((return blockDim.y;))
    KOKKOS_IF_ON_HOST((return 1;))
  }

  KOKKOS_INLINE_FUNCTION void team_barrier() const {
    KOKKOS_IF_ON_DEVICE((
        if (1 == blockDim.z) { __syncthreads(); }  // team == block
        else { __threadfence_block(); }            // team <= warp
        ))
  }

  //--------------------------------------------------------------------------

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(ValueType& val,
                                             const int& thread_id) const {
    (void)val;
    (void)thread_id;
    KOKKOS_IF_ON_DEVICE((
        if (1 == blockDim.z) {  // team == block
          __syncthreads();
          // Wait for shared data write until all threads arrive here
          if (threadIdx.x == 0u && threadIdx.y == (uint32_t)thread_id) {
            *((ValueType*)m_team_reduce) = val;
          }
          __syncthreads();  // Wait for shared data read until root thread
                            // writes
          val = *((ValueType*)m_team_reduce);
        } else {               // team <= warp
          ValueType tmp(val);  // input might not be a register variable
          Impl::in_place_shfl(val, tmp, blockDim.x * thread_id,
                              blockDim.x * blockDim.y);
        }))
  }

  template <class Closure, class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(Closure const& f, ValueType& val,
                                             const int& thread_id) const {
    (void)f;
    (void)val;
    (void)thread_id;
    KOKKOS_IF_ON_DEVICE((
        f(val);

        if (1 == blockDim.z) {  // team == block
          __syncthreads();
          // Wait for shared data write until all threads arrive here
          if (threadIdx.x == 0u && threadIdx.y == (uint32_t)thread_id) {
            *((ValueType*)m_team_reduce) = val;
          }
          __syncthreads();  // Wait for shared data read until root thread
                            // writes
          val = *((ValueType*)m_team_reduce);
        } else {               // team <= warp
          ValueType tmp(val);  // input might not be a register variable
          Impl::in_place_shfl(val, tmp, blockDim.x * thread_id,
                              blockDim.x * blockDim.y);
        }))
  }

  //--------------------------------------------------------------------------
  /**\brief  Reduction across a team
   *
   *  Mapping of teams onto blocks:
   *    blockDim.x  is "vector lanes"
   *    blockDim.y  is team "threads"
   *    blockDim.z  is number of teams per block
   *
   *  Requires:
   *    blockDim.x is power two
   *    blockDim.x <= CudaTraits::WarpSize
   *    ( 0 == CudaTraits::WarpSize % ( blockDim.x * blockDim.y )
   *      OR
   *    ( 1 == blockDim.z )
   */
  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<is_reducer<ReducerType>::value>
  team_reduce(ReducerType const& reducer) const noexcept {
    team_reduce(reducer, reducer.reference());
  }

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<is_reducer<ReducerType>::value>
  team_reduce(ReducerType const& reducer,
              typename ReducerType::value_type& value) const noexcept {
    (void)reducer;
    (void)value;
    KOKKOS_IF_ON_DEVICE(
        (typename Impl::FunctorAnalysis<Impl::FunctorPatternInterface::REDUCE,
                                        TeamPolicy<Cuda>, ReducerType>::Reducer
             wrapped_reducer(&reducer);
         cuda_intra_block_reduction(value, wrapped_reducer, blockDim.y);
         reducer.reference() = value;))
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
    KOKKOS_IF_ON_DEVICE((
        Type* const base_data = (Type*)m_team_reduce;

        __syncthreads();  // Don't write in to shared data until all threads
                          // have entered this function

        if (0 == threadIdx.y) { base_data[0] = 0; }

        base_data[threadIdx.y + 1] = value;
        Impl::CudaJoinFunctor<Type> cuda_join_functor;
        typename Impl::FunctorAnalysis<
            Impl::FunctorPatternInterface::SCAN, TeamPolicy<Cuda>,
            Impl::CudaJoinFunctor<Type>>::Reducer reducer(&cuda_join_functor);
        Impl::cuda_intra_block_reduce_scan<true>(reducer, base_data + 1);

        if (global_accum) {
          if (blockDim.y == threadIdx.y + 1) {
            base_data[blockDim.y] =
                atomic_fetch_add(global_accum, base_data[blockDim.y]);
          }
          __syncthreads();  // Wait for atomic
          base_data[threadIdx.y] += base_data[blockDim.y];
        }

        return base_data[threadIdx.y];))

    KOKKOS_IF_ON_HOST(((void)value; (void)global_accum; return Type();))
  }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template <typename Type>
  KOKKOS_INLINE_FUNCTION Type team_scan(const Type& value) const {
    return this->template team_scan<Type>(value, nullptr);
  }

  //----------------------------------------

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION static std::enable_if_t<is_reducer<ReducerType>::value>
  vector_reduce(ReducerType const& reducer) {
    vector_reduce(reducer, reducer.reference());
  }

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION static std::enable_if_t<is_reducer<ReducerType>::value>
  vector_reduce(ReducerType const& reducer,
                typename ReducerType::value_type& value) {
    (void)reducer;
    (void)value;
    KOKKOS_IF_ON_DEVICE(
        (if (blockDim.x == 1) return;

         // Intra vector lane shuffle reduction:
         typename ReducerType::value_type tmp(value);
         typename ReducerType::value_type tmp2 = tmp;

         unsigned mask =
             blockDim.x == 32
                 ? 0xffffffff
                 : ((1 << blockDim.x) - 1)
                       << ((threadIdx.y % (32 / blockDim.x)) * blockDim.x);

         for (int i = blockDim.x; (i >>= 1);) {
           Impl::in_place_shfl_down(tmp2, tmp, i, blockDim.x, mask);
           if ((int)threadIdx.x < i) {
             reducer.join(tmp, tmp2);
           }
         }

         // Broadcast from root lane to all other lanes.
         // Cannot use "butterfly" algorithm to avoid the broadcast
         // because floating point summation is not associative
         // and thus different threads could have different results.

         Impl::in_place_shfl(tmp2, tmp, 0, blockDim.x, mask);
         value = tmp2; reducer.reference() = tmp2;))
  }

  //----------------------------------------
  // Private for the driver

  KOKKOS_INLINE_FUNCTION
  CudaTeamMember(void* shared, const size_t shared_begin,
                 const size_t shared_size, void* scratch_level_1_ptr,
                 const size_t scratch_level_1_size, const int arg_league_rank,
                 const int arg_league_size)
      : m_team_reduce(shared),
        m_team_shared(static_cast<char*>(shared) + shared_begin, shared_size,
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
struct TeamThreadRangeBoundariesStruct<iType, CudaTeamMember> {
  using index_type = iType;
  const CudaTeamMember& member;
  const iType start;
  const iType end;

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const CudaTeamMember& thread_, iType count)
      : member(thread_), start(0), end(count) {}

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const CudaTeamMember& thread_, iType begin_,
                                  iType end_)
      : member(thread_), start(begin_), end(end_) {}
};

template <typename iType>
struct TeamVectorRangeBoundariesStruct<iType, CudaTeamMember> {
  using index_type = iType;
  const CudaTeamMember& member;
  const iType start;
  const iType end;

  KOKKOS_INLINE_FUNCTION
  TeamVectorRangeBoundariesStruct(const CudaTeamMember& thread_,
                                  const iType& count)
      : member(thread_), start(0), end(count) {}

  KOKKOS_INLINE_FUNCTION
  TeamVectorRangeBoundariesStruct(const CudaTeamMember& thread_,
                                  const iType& begin_, const iType& end_)
      : member(thread_), start(begin_), end(end_) {}
};

template <typename iType>
struct ThreadVectorRangeBoundariesStruct<iType, CudaTeamMember> {
  using index_type = iType;
  const index_type start;
  const index_type end;

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const CudaTeamMember, index_type count)
      : start(static_cast<index_type>(0)), end(count) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(index_type count)
      : start(static_cast<index_type>(0)), end(count) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const CudaTeamMember, index_type arg_begin,
                                    index_type arg_end)
      : start(arg_begin), end(arg_end) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(index_type arg_begin, index_type arg_end)
      : start(arg_begin), end(arg_end) {}
};

}  // namespace Impl

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamThreadRangeBoundariesStruct<iType, Impl::CudaTeamMember>
    TeamThreadRange(const Impl::CudaTeamMember& thread, iType count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::CudaTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::CudaTeamMember>
TeamThreadRange(const Impl::CudaTeamMember& thread, iType1 begin, iType2 end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::CudaTeamMember>(
      thread, iType(begin), iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamVectorRangeBoundariesStruct<iType, Impl::CudaTeamMember>
    TeamVectorRange(const Impl::CudaTeamMember& thread, const iType& count) {
  return Impl::TeamVectorRangeBoundariesStruct<iType, Impl::CudaTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamVectorRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::CudaTeamMember>
TeamVectorRange(const Impl::CudaTeamMember& thread, const iType1& begin,
                const iType2& end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::TeamVectorRangeBoundariesStruct<iType, Impl::CudaTeamMember>(
      thread, iType(begin), iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::CudaTeamMember>
    ThreadVectorRange(const Impl::CudaTeamMember& thread, iType count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::CudaTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<
    std::common_type_t<iType1, iType2>, Impl::CudaTeamMember>
ThreadVectorRange(const Impl::CudaTeamMember& thread, iType1 arg_begin,
                  iType2 arg_end) {
  using iType = std::common_type_t<iType1, iType2>;
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::CudaTeamMember>(
      thread, iType(arg_begin), iType(arg_end));
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::CudaTeamMember> PerTeam(
    const Impl::CudaTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::CudaTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::CudaTeamMember> PerThread(
    const Impl::CudaTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::CudaTeamMember>(thread);
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
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::CudaTeamMember>&
        loop_boundaries,
    const Closure& closure) {
  (void)loop_boundaries;
  (void)closure;
  KOKKOS_IF_ON_DEVICE(
      (for (iType i = loop_boundaries.start + threadIdx.y;
            i < loop_boundaries.end; i += blockDim.y) { closure(i); }))
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
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::CudaTeamMember>& loop_boundaries,
                const Closure& closure, const ReducerType& reducer) {
  (void)loop_boundaries;
  (void)closure;
  (void)reducer;
  KOKKOS_IF_ON_DEVICE(
      (typename ReducerType::value_type value;

       reducer.init(value);

       for (iType i = loop_boundaries.start + threadIdx.y;
            i < loop_boundaries.end; i += blockDim.y) { closure(i, value); }

       loop_boundaries.member.team_reduce(reducer, value);))
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
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer<ValueType>::value>
parallel_reduce(const Impl::TeamThreadRangeBoundariesStruct<
                    iType, Impl::CudaTeamMember>& loop_boundaries,
                const Closure& closure, ValueType& result) {
  (void)loop_boundaries;
  (void)closure;
  (void)result;
  KOKKOS_IF_ON_DEVICE(
      (ValueType val; Kokkos::Sum<ValueType> reducer(val);

       reducer.init(reducer.reference());

       for (iType i = loop_boundaries.start + threadIdx.y;
            i < loop_boundaries.end; i += blockDim.y) { closure(i, val); }

       loop_boundaries.member.team_reduce(reducer, val);
       result = reducer.reference();))
}

template <typename iType, class Closure>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamVectorRangeBoundariesStruct<iType, Impl::CudaTeamMember>&
        loop_boundaries,
    const Closure& closure) {
  (void)loop_boundaries;
  (void)closure;
  KOKKOS_IF_ON_DEVICE((for (iType i = loop_boundaries.start +
                                      threadIdx.y * blockDim.x + threadIdx.x;
                            i < loop_boundaries.end;
                            i += blockDim.y * blockDim.x) { closure(i); }))
}

template <typename iType, class Closure, class ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                    iType, Impl::CudaTeamMember>& loop_boundaries,
                const Closure& closure, const ReducerType& reducer) {
  (void)loop_boundaries;
  (void)closure;
  (void)reducer;
  KOKKOS_IF_ON_DEVICE((typename ReducerType::value_type value;
                       reducer.init(value);

                       for (iType i = loop_boundaries.start +
                                      threadIdx.y * blockDim.x + threadIdx.x;
                            i < loop_boundaries.end;
                            i += blockDim.y * blockDim.x) { closure(i, value); }

                       loop_boundaries.member.vector_reduce(reducer, value);
                       loop_boundaries.member.team_reduce(reducer, value);))
}

template <typename iType, class Closure, typename ValueType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<!Kokkos::is_reducer<ValueType>::value>
parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                    iType, Impl::CudaTeamMember>& loop_boundaries,
                const Closure& closure, ValueType& result) {
  (void)loop_boundaries;
  (void)closure;
  (void)result;
  KOKKOS_IF_ON_DEVICE((ValueType val; Kokkos::Sum<ValueType> reducer(val);

                       reducer.init(reducer.reference());

                       for (iType i = loop_boundaries.start +
                                      threadIdx.y * blockDim.x + threadIdx.x;
                            i < loop_boundaries.end;
                            i += blockDim.y * blockDim.x) { closure(i, val); }

                       loop_boundaries.member.vector_reduce(reducer);
                       loop_boundaries.member.team_reduce(reducer);
                       result = reducer.reference();))
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
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::CudaTeamMember>&
        loop_boundaries,
    const Closure& closure) {
  (void)loop_boundaries;
  (void)closure;
  KOKKOS_IF_ON_DEVICE((
      for (iType i = loop_boundaries.start + threadIdx.x;
           i < loop_boundaries.end; i += blockDim.x) { closure(i); }

      __syncwarp(blockDim.x == 32
                     ? 0xffffffff
                     : ((1 << blockDim.x) - 1)
                           << (threadIdx.y % (32 / blockDim.x)) * blockDim.x);))
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
KOKKOS_INLINE_FUNCTION std::enable_if_t<is_reducer<ReducerType>::value>
parallel_reduce(Impl::ThreadVectorRangeBoundariesStruct<
                    iType, Impl::CudaTeamMember> const& loop_boundaries,
                Closure const& closure, ReducerType const& reducer) {
  (void)loop_boundaries;
  (void)closure;
  (void)reducer;
  KOKKOS_IF_ON_DEVICE((

      reducer.init(reducer.reference());

      for (iType i = loop_boundaries.start + threadIdx.x;
           i < loop_boundaries.end;
           i += blockDim.x) { closure(i, reducer.reference()); }

      Impl::CudaTeamMember::vector_reduce(reducer);

      ))
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
KOKKOS_INLINE_FUNCTION std::enable_if_t<!is_reducer<ValueType>::value>
parallel_reduce(Impl::ThreadVectorRangeBoundariesStruct<
                    iType, Impl::CudaTeamMember> const& loop_boundaries,
                Closure const& closure, ValueType& result) {
  (void)loop_boundaries;
  (void)closure;
  (void)result;
  KOKKOS_IF_ON_DEVICE(
      (result = ValueType();

       for (iType i = loop_boundaries.start + threadIdx.x;
            i < loop_boundaries.end; i += blockDim.x) { closure(i, result); }

       Impl::CudaTeamMember::vector_reduce(Kokkos::Sum<ValueType>(result));

       ))
}

//----------------------------------------------------------------------------

/** \brief  Inter-thread parallel exclusive prefix sum.
 *
 *  Executes closure(iType i, ValueType & val, bool final) for each i=[0..N)
 *
 *  The range [0..N) is mapped to each rank in the team (whose global rank is
 *  less than N) and a scan operation is performed. The last call to closure has
 *  final == true.
 */
// This is the same code as in HIP and largely the same as in OpenMPTarget
template <typename iType, typename FunctorType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::CudaTeamMember>&
        loop_bounds,
    const FunctorType& lambda) {
  // Extract value_type from lambda
  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void,
      FunctorType>::value_type;

  const auto start     = loop_bounds.start;
  const auto end       = loop_bounds.end;
  auto& member         = loop_bounds.member;
  const auto team_size = member.team_size();
  const auto team_rank = member.team_rank();
  const auto nchunk    = (end - start + team_size - 1) / team_size;
  value_type accum     = 0;
  // each team has to process one or more chunks of the prefix scan
  for (iType i = 0; i < nchunk; ++i) {
    auto ii = start + i * team_size + team_rank;
    // local accumulation for this chunk
    value_type local_accum = 0;
    // user updates value with prefix value
    if (ii < loop_bounds.end) lambda(ii, local_accum, false);
    // perform team scan
    local_accum = member.team_scan(local_accum);
    // add this blocks accum to total accumulation
    auto val = accum + local_accum;
    // user updates their data with total accumulation
    if (ii < loop_bounds.end) lambda(ii, val, true);
    // the last value needs to be propogated to next chunk
    if (team_rank == team_size - 1) accum = val;
    // broadcast last value to rest of the team
    member.team_broadcast(accum, team_size - 1);
  }
}

//----------------------------------------------------------------------------

/** \brief  Intra-thread vector parallel scan with reducer.
 *
 *  Executes closure(iType i, ValueType & val, bool final) for each i=[0..N)
 *
 *  The range [0..N) is mapped to all vector lanes in the
 *  thread and a scan operation is performed.
 *  The last call to closure has final == true.
 */
template <typename iType, class Closure, typename ReducerType>
KOKKOS_INLINE_FUNCTION std::enable_if_t<Kokkos::is_reducer<ReducerType>::value>
parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<
                  iType, Impl::CudaTeamMember>& loop_boundaries,
              const Closure& closure, const ReducerType& reducer) {
  (void)loop_boundaries;
  (void)closure;
  (void)reducer;
  KOKKOS_IF_ON_DEVICE((

      using value_type = typename ReducerType::value_type;

      value_type accum;

      reducer.init(accum);

      const value_type identity = accum;

      // Loop through boundaries by vector-length chunks
      // must scan at each iteration

      // All thread "lanes" must loop the same number of times.
      // Determine an loop end for all thread "lanes."
      // Requires:
      //   blockDim.x is power of two and thus
      //     ( end % blockDim.x ) == ( end & ( blockDim.x - 1 ) )
      //   1 <= blockDim.x <= CudaTraits::WarpSize

      const int mask = blockDim.x - 1;
      const unsigned active_mask =
          blockDim.x == 32
              ? 0xffffffff
              : ((1 << blockDim.x) - 1)
                    << (threadIdx.y % (32 / blockDim.x)) * blockDim.x;
      const int rem = loop_boundaries.end & mask;  // == end % blockDim.x
      const int end = loop_boundaries.end + (rem ? blockDim.x - rem : 0);

      for (int i = threadIdx.x; i < end; i += blockDim.x) {
        value_type val = identity;

        // First acquire per-lane contributions.
        // This sets i's val to i-1's contribution
        // to make the latter in_place_shfl_up an
        // exclusive scan -- the final accumulation
        // of i's val will be included in the second
        // closure call later.
        if (i < loop_boundaries.end && threadIdx.x > 0) {
          closure(i - 1, val, false);
        }

        // Bottom up exclusive scan in triangular pattern
        // where each CUDA thread is the root of a reduction tree
        // from the zeroth "lane" to itself.
        //  [t] += [t-1] if t >= 1
        //  [t] += [t-2] if t >= 2
        //  [t] += [t-4] if t >= 4
        //  ...
        //  This differs from the non-reducer overload, where an inclusive scan
        //  was implemented, because in general the binary operator cannot be
        //  inverted and we would not be able to remove the inclusive
        //  contribution by inversion.
        for (int j = 1; j < (int)blockDim.x; j <<= 1) {
          value_type tmp = identity;
          Impl::in_place_shfl_up(tmp, val, j, blockDim.x, active_mask);
          if (j <= (int)threadIdx.x) {
            reducer.join(val, tmp);
          }
        }

        // Include accumulation
        reducer.join(val, accum);

        // Update i's contribution into the val
        // and add it to accum for next round
        if (i < loop_boundaries.end) closure(i, val, true);
        Impl::in_place_shfl(accum, val, mask, blockDim.x, active_mask);
      }

      ))
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
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::CudaTeamMember>&
        loop_boundaries,
    const Closure& closure) {
  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, Closure>::value_type;
  value_type dummy;
  parallel_scan(loop_boundaries, closure, Kokkos::Sum<value_type>(dummy));
}

}  // namespace Kokkos

namespace Kokkos {

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::CudaTeamMember>&,
    const FunctorType& lambda) {
  (void)lambda;
  KOKKOS_IF_ON_DEVICE((
      if (threadIdx.x == 0) { lambda(); }

      __syncwarp(blockDim.x == 32
                     ? 0xffffffff
                     : ((1 << blockDim.x) - 1)
                           << (threadIdx.y % (32 / blockDim.x)) * blockDim.x);))
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::CudaTeamMember>&,
    const FunctorType& lambda) {
  (void)lambda;
  KOKKOS_IF_ON_DEVICE((
      if (threadIdx.x == 0 && threadIdx.y == 0) { lambda(); }

      __syncwarp(blockDim.x == 32
                     ? 0xffffffff
                     : ((1 << blockDim.x) - 1)
                           << (threadIdx.y % (32 / blockDim.x)) * blockDim.x);))
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::CudaTeamMember>&,
    const FunctorType& lambda, ValueType& val) {
  (void)lambda;
  (void)val;
  KOKKOS_IF_ON_DEVICE(
      (if (threadIdx.x == 0) { lambda(val); }

       unsigned mask =
           blockDim.x == 32
               ? 0xffffffff
               : ((1 << blockDim.x) - 1)
                     << ((threadIdx.y % (32 / blockDim.x)) * blockDim.x);

       Impl::in_place_shfl(val, val, 0, blockDim.x, mask);))
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::CudaTeamMember>& single_struct,
    const FunctorType& lambda, ValueType& val) {
  (void)single_struct;
  (void)lambda;
  (void)val;
  KOKKOS_IF_ON_DEVICE(
      (if (threadIdx.x == 0 && threadIdx.y == 0) { lambda(val); }

       single_struct.team_member.team_broadcast(val, 0);))
}

}  // namespace Kokkos

#endif /* defined(KOKKOS_ENABLE_CUDA) */

#endif /* #ifndef KOKKOS_CUDA_TEAM_HPP */
