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

#ifndef KOKKOS_SYCL_TEAM_HPP
#define KOKKOS_SYCL_TEAM_HPP

#include <Kokkos_Macros.hpp>

#ifdef KOKKOS_ENABLE_SYCL

#include <utility>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/**\brief  Team member_type passed to TeamPolicy or TeamTask closures.
 */
class SYCLTeamMember {
 public:
  using execution_space      = Kokkos::Experimental::SYCL;
  using scratch_memory_space = execution_space::scratch_memory_space;

 private:
  mutable void* m_team_reduce;
  scratch_memory_space m_team_shared;
  int m_team_reduce_size;
  sycl::nd_item<2> m_item;

 public:
  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_shmem() const {
    return m_team_shared.set_team_thread_mode(0, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& team_scratch(
      const int level) const {
    return m_team_shared.set_team_thread_mode(level, 1, 0);
  }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space& thread_scratch(
      const int level) const {
    return m_team_shared.set_team_thread_mode(level, team_size(), team_rank());
  }

  KOKKOS_INLINE_FUNCTION int league_rank() const {
    return m_item.get_group_linear_id();
  }
  KOKKOS_INLINE_FUNCTION int league_size() const {
    return m_item.get_group_range(1);
  }
  KOKKOS_INLINE_FUNCTION int team_rank() const {
    return m_item.get_local_id(0);
  }
  KOKKOS_INLINE_FUNCTION int team_size() const {
    return m_item.get_local_range(0);
  }
  KOKKOS_INLINE_FUNCTION void team_barrier() const {
    sycl::group_barrier(m_item.get_group());
  }

  KOKKOS_INLINE_FUNCTION const sycl::nd_item<2>& item() const { return m_item; }

  //--------------------------------------------------------------------------

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_arithmetic_v<ValueType>>
  team_broadcast(ValueType& val, const int thread_id) const {
    val = sycl::group_broadcast(m_item.get_group(), val,
                                sycl::id<2>(thread_id, 0));
  }

  // FIXME_SYCL remove/adapt this overload once the Intel oneAPI implementation
  // is conforming to the SYCL2020 standard (allowing trivially-copyable types)
  template <class ValueType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<!std::is_arithmetic_v<ValueType>>
  team_broadcast(ValueType& val, const int thread_id) const {
    // Wait for shared data write until all threads arrive here
    sycl::group_barrier(m_item.get_group());
    if (m_item.get_local_id(1) == 0 &&
        static_cast<int>(m_item.get_local_id(0)) == thread_id) {
      *static_cast<ValueType*>(m_team_reduce) = val;
    }
    // Wait for shared data read until root thread writes
    sycl::group_barrier(m_item.get_group());
    val = *(static_cast<ValueType*>(m_team_reduce));
  }

  template <class Closure, class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(Closure const& f, ValueType& val,
                                             const int thread_id) const {
    f(val);
    team_broadcast(val, thread_id);
  }

  //--------------------------------------------------------------------------
  /**\brief  Reduction across a team
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
    using value_type = typename ReducerType::value_type;

    auto sg                       = m_item.get_sub_group();
    const auto sub_group_range    = sg.get_local_range()[0];
    const auto vector_range       = m_item.get_local_range(1);
    const unsigned int team_size_ = team_size();
    const unsigned int team_rank_ = team_rank();

    // First combine the values in the same subgroup
    for (unsigned int shift = 1; vector_range * shift < sub_group_range;
         shift <<= 1) {
      const value_type tmp = sg.shuffle_down(value, vector_range * shift);
      if (team_rank_ + shift < team_size_) reducer.join(value, tmp);
    }
    value = sg.shuffle(value, 0);

    // We need to chunk up the whole reduction because we might not have
    // allocated enough memory.
    const auto n_subgroups = sg.get_group_range()[0];
    const unsigned int maximum_work_range =
        std::min<int>(m_team_reduce_size / sizeof(value_type), n_subgroups);

    const auto id_in_sg  = sg.get_local_id()[0];
    auto reduction_array = static_cast<value_type*>(m_team_reduce);

    // Load values into the first maximum_work_range values of the reduction
    // array in chunks. This means that only sub groups with an id in the
    // corresponding chunk load values.
    const auto group_id = sg.get_group_id()[0];
    if (id_in_sg == 0 && group_id < maximum_work_range)
      reduction_array[group_id] = value;
    sycl::group_barrier(m_item.get_group());

    for (unsigned int start = maximum_work_range; start < n_subgroups;
         start += maximum_work_range) {
      if (id_in_sg == 0 && group_id >= start &&
          group_id <
              std::min<unsigned int>(start + maximum_work_range, n_subgroups))
        reducer.join(reduction_array[group_id - start], value);
      sycl::group_barrier(m_item.get_group());
    }

    // Let the first subgroup do the final reduction
    if (group_id == 0) {
      const auto local_range = sg.get_local_range()[0];
      auto result =
          reduction_array[id_in_sg < maximum_work_range ? id_in_sg : 0];
      // In case the maximum_work_range is larger than the range of the first
      // subgroup, we first combine the items with a higher index.
      for (unsigned int offset = local_range; offset < maximum_work_range;
           offset += local_range)
        if (id_in_sg + offset < maximum_work_range)
          reducer.join(result, reduction_array[id_in_sg + offset]);
      sycl::group_barrier(sg);

      // Now do the actual subgroup reduction.
      const auto min_range =
          std::min<unsigned int>(maximum_work_range, local_range);
      for (unsigned int stride = 1; stride < min_range; stride <<= 1) {
        const auto tmp = sg.shuffle_down(result, stride);
        if (id_in_sg + stride < min_range) reducer.join(result, tmp);
      }
      if (id_in_sg == 0) reduction_array[0] = result;
    }
    sycl::group_barrier(m_item.get_group());

    reducer.reference() = reduction_array[0];
    // Make sure that the reduction array hasn't been modified in the meantime.
    m_item.barrier(sycl::access::fence_space::local_space);
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
  KOKKOS_INLINE_FUNCTION Type team_scan(const Type& input_value,
                                        Type* const global_accum) const {
    Type value                 = input_value;
    auto sg                    = m_item.get_sub_group();
    const auto sub_group_range = sg.get_local_range()[0];
    const auto vector_range    = m_item.get_local_range(1);
    const auto id_in_sg        = sg.get_local_id()[0];

    // First combine the values in the same subgroup
    for (unsigned int stride = 1; vector_range * stride < sub_group_range;
         stride <<= 1) {
      auto tmp = sg.shuffle_up(value, vector_range * stride);
      if (id_in_sg >= vector_range * stride) value += tmp;
    }

    const auto n_active_subgroups = sg.get_group_range()[0];
    const auto base_data          = static_cast<Type*>(m_team_reduce);
    if (static_cast<int>(n_active_subgroups * sizeof(Type)) >
        m_team_reduce_size)
      Kokkos::abort("Not implemented!");

    const auto group_id = sg.get_group_id()[0];
    if (id_in_sg == sub_group_range - 1) base_data[group_id] = value;
    sycl::group_barrier(m_item.get_group());

    // scan subgroup results using the first subgroup
    if (n_active_subgroups > 1) {
      if (group_id == 0) {
        const auto n_rounds =
            (n_active_subgroups + sub_group_range - 1) / sub_group_range;
        for (unsigned int round = 0; round < n_rounds; ++round) {
          const auto idx         = id_in_sg + round * sub_group_range;
          const auto upper_bound = std::min(
              sub_group_range, n_active_subgroups - round * sub_group_range);
          auto local_value = base_data[idx];
          for (unsigned int stride = 1; stride < upper_bound; stride <<= 1) {
            auto tmp = sg.shuffle_up(local_value, stride);
            if (id_in_sg >= stride) {
              if (idx < n_active_subgroups)
                local_value += tmp;
              else
                local_value = tmp;
            }
          }
          base_data[idx] = local_value;
          if (round > 0)
            base_data[idx] += base_data[round * sub_group_range - 1];
          if (round + 1 < n_rounds) sycl::group_barrier(sg);
        }
      }
      sycl::group_barrier(m_item.get_group());
    }
    auto total = base_data[n_active_subgroups - 1];

    const auto update = sg.shuffle_up(value, vector_range);
    Type intermediate = (group_id > 0 ? base_data[group_id - 1] : 0) +
                        (id_in_sg >= vector_range ? update : 0);

    if (global_accum) {
      if (id_in_sg == sub_group_range - 1 &&
          group_id == n_active_subgroups - 1) {
        base_data[n_active_subgroups - 1] =
            atomic_fetch_add(global_accum, total);
      }
      sycl::group_barrier(m_item.get_group());  // Wait for atomic
      intermediate += base_data[n_active_subgroups - 1];
    }
    // Make sure that the reduction array hasn't been modified in the meantime.
    m_item.barrier(sycl::access::fence_space::local_space);

    return intermediate;
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
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      vector_reduce(ReducerType const& reducer) const {
    vector_reduce(reducer, reducer.reference());
  }

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      vector_reduce(ReducerType const& reducer,
                    typename ReducerType::value_type& value) const {
    const auto tidx1   = m_item.get_local_id(1);
    const auto grange1 = m_item.get_local_range(1);

    const auto sg = m_item.get_sub_group();

    if (grange1 == 1) return;

    // Intra vector lane shuffle reduction:
    typename ReducerType::value_type tmp(value);
    typename ReducerType::value_type tmp2 = tmp;

    for (int i = grange1; (i >>= 1);) {
      tmp2 = sg.shuffle_down(tmp, i);
      if (static_cast<int>(tidx1) < i) {
        reducer.join(tmp, tmp2);
      }
    }

    // Broadcast from root lane to all other lanes.
    // Cannot use "butterfly" algorithm to avoid the broadcast
    // because floating point summation is not associative
    // and thus different threads could have different results.

    tmp2  = sg.shuffle(tmp, (sg.get_local_id() / grange1) * grange1);
    value = tmp2;
    reducer.reference() = tmp2;
  }

  //----------------------------------------
  // Private for the driver

  KOKKOS_INLINE_FUNCTION
  SYCLTeamMember(void* shared, const int shared_begin, const int shared_size,
                 void* scratch_level_1_ptr, const int scratch_level_1_size,
                 const sycl::nd_item<2> item)
      : m_team_reduce(shared),
        m_team_shared(static_cast<char*>(shared) + shared_begin, shared_size,
                      scratch_level_1_ptr, scratch_level_1_size),
        m_team_reduce_size(shared_begin),
        m_item(item) {}

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
struct TeamThreadRangeBoundariesStruct<iType, SYCLTeamMember> {
  using index_type = iType;
  const SYCLTeamMember& member;
  const iType start;
  const iType end;

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const SYCLTeamMember& thread_, iType count)
      : member(thread_), start(0), end(count) {}

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const SYCLTeamMember& thread_, iType begin_,
                                  iType end_)
      : member(thread_), start(begin_), end(end_) {}
};

template <typename iType>
struct TeamVectorRangeBoundariesStruct<iType, SYCLTeamMember> {
  using index_type = iType;
  const SYCLTeamMember& member;
  const iType start;
  const iType end;

  KOKKOS_INLINE_FUNCTION
  TeamVectorRangeBoundariesStruct(const SYCLTeamMember& thread_,
                                  const iType& count)
      : member(thread_), start(0), end(count) {}

  KOKKOS_INLINE_FUNCTION
  TeamVectorRangeBoundariesStruct(const SYCLTeamMember& thread_,
                                  const iType& begin_, const iType& end_)
      : member(thread_), start(begin_), end(end_) {}
};

template <typename iType>
struct ThreadVectorRangeBoundariesStruct<iType, SYCLTeamMember> {
  using index_type = iType;
  const SYCLTeamMember& member;
  const index_type start;
  const index_type end;

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const SYCLTeamMember& thread,
                                    index_type count)
      : member(thread), start(static_cast<index_type>(0)), end(count) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct(const SYCLTeamMember& thread,
                                    index_type arg_begin, index_type arg_end)
      : member(thread), start(arg_begin), end(arg_end) {}
};

}  // namespace Impl

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamThreadRangeBoundariesStruct<iType, Impl::SYCLTeamMember>
    TeamThreadRange(const Impl::SYCLTeamMember& thread, iType count) {
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::SYCLTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamThreadRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, Impl::SYCLTeamMember>
TeamThreadRange(const Impl::SYCLTeamMember& thread, iType1 begin, iType2 end) {
  using iType = typename std::common_type<iType1, iType2>::type;
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::SYCLTeamMember>(
      thread, iType(begin), iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::TeamVectorRangeBoundariesStruct<iType, Impl::SYCLTeamMember>
    TeamVectorRange(const Impl::SYCLTeamMember& thread, const iType& count) {
  return Impl::TeamVectorRangeBoundariesStruct<iType, Impl::SYCLTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::TeamVectorRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, Impl::SYCLTeamMember>
TeamVectorRange(const Impl::SYCLTeamMember& thread, const iType1& begin,
                const iType2& end) {
  using iType = typename std::common_type<iType1, iType2>::type;
  return Impl::TeamVectorRangeBoundariesStruct<iType, Impl::SYCLTeamMember>(
      thread, iType(begin), iType(end));
}

template <typename iType>
KOKKOS_INLINE_FUNCTION
    Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::SYCLTeamMember>
    ThreadVectorRange(const Impl::SYCLTeamMember& thread, iType count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::SYCLTeamMember>(
      thread, count);
}

template <typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION Impl::ThreadVectorRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, Impl::SYCLTeamMember>
ThreadVectorRange(const Impl::SYCLTeamMember& thread, iType1 arg_begin,
                  iType2 arg_end) {
  using iType = typename std::common_type<iType1, iType2>::type;
  return Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::SYCLTeamMember>(
      thread, iType(arg_begin), iType(arg_end));
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::SYCLTeamMember> PerTeam(
    const Impl::SYCLTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::SYCLTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::SYCLTeamMember> PerThread(
    const Impl::SYCLTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::SYCLTeamMember>(thread);
}

//----------------------------------------------------------------------------

/** \brief  Inter-thread parallel_for.
 *
 *  Executes closure(iType i) for each i=[0..N).
 *
 * The range [0..N) is mapped to all threads of the calling thread team.
 */
template <typename iType, class Closure>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::SYCLTeamMember>&
        loop_boundaries,
    const Closure& closure) {
  for (iType i = loop_boundaries.start +
                 loop_boundaries.member.item().get_local_id(0);
       i < loop_boundaries.end;
       i += loop_boundaries.member.item().get_local_range(0))
    closure(i);
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
                        iType, Impl::SYCLTeamMember>& loop_boundaries,
                    const Closure& closure, const ReducerType& reducer) {
  typename ReducerType::value_type value;
  reducer.init(value);

  for (iType i = loop_boundaries.start +
                 loop_boundaries.member.item().get_local_id(0);
       i < loop_boundaries.end;
       i += loop_boundaries.member.item().get_local_range(0)) {
    closure(i, value);
  }

  loop_boundaries.member.team_reduce(reducer, value);
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
                        iType, Impl::SYCLTeamMember>& loop_boundaries,
                    const Closure& closure, ValueType& result) {
  ValueType val;
  Kokkos::Sum<ValueType> reducer(val);

  reducer.init(reducer.reference());

  for (iType i = loop_boundaries.start +
                 loop_boundaries.member.item().get_local_id(0);
       i < loop_boundaries.end;
       i += loop_boundaries.member.item().get_local_range(0)) {
    closure(i, val);
  }

  loop_boundaries.member.team_reduce(reducer, val);
  result = reducer.reference();
}

/** \brief  Inter-thread parallel exclusive prefix sum.
 *
 *  Executes closure(iType i, ValueType & val, bool final) for each i=[0..N)
 *
 *  The range [0..N) is mapped to each rank in the team (whose global rank is
 *  less than N) and a scan operation is performed. The last call to closure has
 *  final == true.
 */
// This is the same code as in CUDA and largely the same as in OpenMPTarget
template <typename iType, typename FunctorType>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::SYCLTeamMember>&
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

template <typename iType, class Closure>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamVectorRangeBoundariesStruct<iType, Impl::SYCLTeamMember>&
        loop_boundaries,
    const Closure& closure) {
  const iType tidx0 = loop_boundaries.member.item().get_local_id(0);
  const iType tidx1 = loop_boundaries.member.item().get_local_id(1);

  const iType grange0 = loop_boundaries.member.item().get_local_range(0);
  const iType grange1 = loop_boundaries.member.item().get_local_range(1);

  for (iType i = loop_boundaries.start + tidx0 * grange1 + tidx1;
       i < loop_boundaries.end; i += grange0 * grange1)
    closure(i);
}

template <typename iType, class Closure, class ReducerType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<Kokkos::is_reducer<ReducerType>::value>::type
    parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                        iType, Impl::SYCLTeamMember>& loop_boundaries,
                    const Closure& closure, const ReducerType& reducer) {
  typename ReducerType::value_type value;
  reducer.init(value);

  const iType tidx0 = loop_boundaries.member.item().get_local_id(0);
  const iType tidx1 = loop_boundaries.member.item().get_local_id(1);

  const iType grange0 = loop_boundaries.member.item().get_local_range(0);
  const iType grange1 = loop_boundaries.member.item().get_local_range(1);

  for (iType i = loop_boundaries.start + tidx0 * grange1 + tidx1;
       i < loop_boundaries.end; i += grange0 * grange1)
    closure(i, value);

  loop_boundaries.member.vector_reduce(reducer, value);
  loop_boundaries.member.team_reduce(reducer, value);
}

template <typename iType, class Closure, typename ValueType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!Kokkos::is_reducer<ValueType>::value>::type
    parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                        iType, Impl::SYCLTeamMember>& loop_boundaries,
                    const Closure& closure, ValueType& result) {
  ValueType val;
  Kokkos::Sum<ValueType> reducer(val);

  reducer.init(reducer.reference());

  const iType tidx0 = loop_boundaries.member.item().get_local_id(0);
  const iType tidx1 = loop_boundaries.member.item().get_local_id(1);

  const iType grange0 = loop_boundaries.member.item().get_local_range(0);
  const iType grange1 = loop_boundaries.member.item().get_local_range(1);

  for (iType i = loop_boundaries.start + tidx0 * grange1 + tidx1;
       i < loop_boundaries.end; i += grange0 * grange1)
    closure(i, val);

  loop_boundaries.member.vector_reduce(reducer);
  loop_boundaries.member.team_reduce(reducer);
  result = reducer.reference();
}

//----------------------------------------------------------------------------

/** \brief  Intra-thread vector parallel_for.
 *
 *  Executes closure(iType i) for each i=[0..N)
 *
 * The range [0..N) is mapped to all vector lanes of the calling thread.
 */
template <typename iType, class Closure>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::SYCLTeamMember>&
        loop_boundaries,
    const Closure& closure) {
  const iType tidx1   = loop_boundaries.member.item().get_local_id(1);
  const iType grange1 = loop_boundaries.member.item().get_local_range(1);

  for (iType i = loop_boundaries.start + tidx1; i < loop_boundaries.end;
       i += grange1)
    closure(i);

  // FIXME_SYCL We only should fence active threads here but this not yet
  // available in the compiler. We need https://github.com/intel/llvm/pull/4904
  // or https://github.com/intel/llvm/pull/4903 for that. The current
  // implementation leads to a deadlock only for SYCL+CUDA if not all threads in
  // a subgroup see this barrier. For SYCL on Intel GPUs, the subgroup barrier
  // is essentially a no-op (only a memory fence), though.
  sycl::group_barrier(loop_boundaries.member.item().get_sub_group());
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
                        iType, Impl::SYCLTeamMember> const& loop_boundaries,
                    Closure const& closure, ReducerType const& reducer) {
  reducer.init(reducer.reference());

  const iType tidx1   = loop_boundaries.member.item().get_local_id(1);
  const iType grange1 = loop_boundaries.member.item().get_local_range(1);

  for (iType i = loop_boundaries.start + tidx1; i < loop_boundaries.end;
       i += grange1)
    closure(i, reducer.reference());

  loop_boundaries.member.vector_reduce(reducer);
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
                        iType, Impl::SYCLTeamMember> const& loop_boundaries,
                    Closure const& closure, ValueType& result) {
  result = ValueType();

  const iType tidx1 = loop_boundaries.member.item().get_local_id(1);
  const int grange1 = loop_boundaries.member.item().get_local_range(1);

  for (iType i = loop_boundaries.start + tidx1; i < loop_boundaries.end;
       i += grange1)
    closure(i, result);

  loop_boundaries.member.vector_reduce(Kokkos::Sum<ValueType>(result));
}

//----------------------------------------------------------------------------

/** \brief  Intra-thread vector parallel exclusive prefix sum with reducer.
 *
 *  Executes closure(iType i, ValueType & val, bool final) for each i=[0..N)
 *
 *  The range [0..N) is mapped to all vector lanes in the
 *  thread and a scan operation is performed.
 *  The last call to closure has final == true.
 */
template <typename iType, class Closure, typename ReducerType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<Kokkos::is_reducer<ReducerType>::value>::type
    parallel_scan(const Impl::ThreadVectorRangeBoundariesStruct<
                      iType, Impl::SYCLTeamMember>& loop_boundaries,
                  const Closure& closure, const ReducerType& reducer) {
  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, Closure>::value_type;

  value_type accum;
  reducer.init(accum);
  const value_type identity = accum;

  // Loop through boundaries by vector-length chunks must scan at each iteration

  // All thread "lanes" must loop the same number of times.
  // Determine an loop end for all thread "lanes."
  // Requires:
  //   grange1 is power of two and thus
  //     ( end % grange1 ) == ( end & ( grange1 - 1 ) )
  //   1 <= grange1 <= sub_group size

  const iType tidx1   = loop_boundaries.member.item().get_local_id(1);
  const iType grange1 = loop_boundaries.member.item().get_local_range(1);

  const int mask          = grange1 - 1;
  const int rem           = loop_boundaries.end & mask;  // == end % grange1
  const int end           = loop_boundaries.end + (rem ? grange1 - rem : 0);
  const auto sg           = loop_boundaries.member.item().get_sub_group();
  const int vector_offset = (sg.get_local_id() / grange1) * grange1;

  for (int i = tidx1; i < end; i += grange1) {
    value_type val = identity;

    // First acquire per-lane contributions.
    // This sets i's val to i-1's contribution to make the latter shfl_up an
    // exclusive scan -- the final accumulation of i's val will be included in
    // the second closure call later.
    if (i < loop_boundaries.end && tidx1 > 0) closure(i - 1, val, false);

    // Bottom up exclusive scan in triangular pattern where each SYCL thread is
    // the root of a reduction tree from the zeroth "lane" to itself.
    //  [t] += [t-1] if t >= 1
    //  [t] += [t-2] if t >= 2
    //  [t] += [t-4] if t >= 4
    //  ...
    for (int j = 1; j < static_cast<int>(grange1); j <<= 1) {
      value_type tmp = sg.shuffle_up(val, j);
      if (j <= static_cast<int>(tidx1)) {
        reducer.join(val, tmp);
      }
    }

    // Include accumulation
    reducer.join(val, accum);

    // Update i's contribution into the val and add it to accum for next round
    if (i < loop_boundaries.end) closure(i, val, true);
    accum = sg.shuffle(val, mask + vector_offset);
  }
}

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
    const Impl::ThreadVectorRangeBoundariesStruct<iType, Impl::SYCLTeamMember>&
        loop_boundaries,
    const Closure& closure) {
  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, Closure>::value_type;
  value_type dummy;
  parallel_scan(loop_boundaries, closure, Kokkos::Sum<value_type>{dummy});
}

}  // namespace Kokkos

namespace Kokkos {

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::SYCLTeamMember>& single_struct,
    const FunctorType& lambda) {
  if (single_struct.team_member.item().get_local_id(1) == 0) lambda();
}

template <class FunctorType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::SYCLTeamMember>& single_struct,
    const FunctorType& lambda) {
  if (single_struct.team_member.item().get_local_linear_id() == 0) lambda();
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::SYCLTeamMember>& single_struct,
    const FunctorType& lambda, ValueType& val) {
  const sycl::nd_item<2> item = single_struct.team_member.item();
  const auto grange1          = item.get_local_range(1);
  const auto sg               = item.get_sub_group();
  if (item.get_local_id(1) == 0) lambda(val);
  val = sg.shuffle(val, (sg.get_local_id() / grange1) * grange1);
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::SYCLTeamMember>& single_struct,
    const FunctorType& lambda, ValueType& val) {
  if (single_struct.team_member.item().get_local_linear_id() == 0) lambda(val);
  single_struct.team_member.team_broadcast(val, 0);
}

}  // namespace Kokkos

#endif

#endif /* #ifndef KOKKOS_SYCL_TEAM_HPP */
