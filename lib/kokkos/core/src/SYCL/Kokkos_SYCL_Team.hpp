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
    // FIXME_SYCL needs to be revised for vector_length>1.
    return m_item.get_group_range(0);
  }
  KOKKOS_INLINE_FUNCTION int team_rank() const {
    return m_item.get_local_linear_id();
  }
  KOKKOS_INLINE_FUNCTION int team_size() const {
    // FIXME_SYCL needs to be revised for vector_length>1.
    return m_item.get_local_range(0);
  }
  KOKKOS_INLINE_FUNCTION void team_barrier() const { m_item.barrier(); }

  KOKKOS_INLINE_FUNCTION const sycl::nd_item<2>& item() const { return m_item; }

  //--------------------------------------------------------------------------

  template <class ValueType>
  KOKKOS_INLINE_FUNCTION void team_broadcast(ValueType& val,
                                             const int thread_id) const {
    // Wait for shared data write until all threads arrive here
    m_item.barrier(sycl::access::fence_space::local_space);
    if (m_item.get_local_id(1) == 0 &&
        static_cast<int>(m_item.get_local_id(0)) == thread_id) {
      *static_cast<ValueType*>(m_team_reduce) = val;
    }
    // Wait for shared data read until root thread writes
    m_item.barrier(sycl::access::fence_space::local_space);
    val = *static_cast<ValueType*>(m_team_reduce);
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

    // We need to chunk up the whole reduction because we might not have
    // allocated enough memory.
    const int maximum_work_range =
        std::min<int>(m_team_reduce_size / sizeof(value_type), team_size());

    int smaller_power_of_two = 1;
    while ((smaller_power_of_two << 1) < maximum_work_range)
      smaller_power_of_two <<= 1;

    const int idx        = team_rank();
    auto reduction_array = static_cast<value_type*>(m_team_reduce);

    // Load values into the first maximum_work_range values of the reduction
    // array in chunks. This means that only threads with an id in the
    // corresponding chunk load values and the reduction is always done by the
    // first smaller_power_of_two threads.
    if (idx < maximum_work_range) reduction_array[idx] = value;
    m_item.barrier(sycl::access::fence_space::local_space);

    for (int start = maximum_work_range; start < team_size();
         start += maximum_work_range) {
      if (idx >= start &&
          idx < std::min(start + maximum_work_range, team_size()))
        reducer.join(reduction_array[idx - start], value);
      m_item.barrier(sycl::access::fence_space::local_space);
    }

    for (int stride = smaller_power_of_two; stride > 0; stride >>= 1) {
      if (idx < stride && idx + stride < maximum_work_range)
        reducer.join(reduction_array[idx], reduction_array[idx + stride]);
      m_item.barrier(sycl::access::fence_space::local_space);
    }
    reducer.reference() = reduction_array[0];
    m_item.barrier(sycl::access::fence_space::local_space);
  }

  // FIXME_SYCL move somewhere else and combine with other places that do
  // parallel_scan
  // Exclusive scan returning the total sum.
  // n is required to be a power of two and
  // temp must point to an array containing the data to be processed
  // The accumulated value is returned.
  template <typename Type>
  static Type prescan(sycl::nd_item<2> m_item, Type* temp, int n) {
    int thid = m_item.get_local_id(0);

    // First do a reduction saving intermediate results
    for (int stride = 1; stride < n; stride <<= 1) {
      auto idx = 2 * stride * (thid + 1) - 1;
      if (idx < n) temp[idx] += temp[idx - stride];
      m_item.barrier(sycl::access::fence_space::local_space);
    }

    Type total_sum = temp[n - 1];
    m_item.barrier(sycl::access::fence_space::local_space);

    // clear the last element so we get an exclusive scan
    if (thid == 0) temp[n - 1] = Type{};
    m_item.barrier(sycl::access::fence_space::local_space);

    // Now add the intermediate results to the remaining items again
    for (int stride = n / 2; stride > 0; stride >>= 1) {
      auto idx = 2 * stride * (thid + 1) - 1;
      if (idx < n) {
        Type dummy         = temp[idx - stride];
        temp[idx - stride] = temp[idx];
        temp[idx] += dummy;
      }
      m_item.barrier(sycl::access::fence_space::local_space);
    }

    return total_sum;
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
    // We need to chunk up the whole reduction because we might not have
    // allocated enough memory.
    const int maximum_work_range =
        std::min<int>(m_team_reduce_size / sizeof(Type), team_size());

    int not_greater_power_of_two = 1;
    while ((not_greater_power_of_two << 1) < maximum_work_range + 1)
      not_greater_power_of_two <<= 1;

    Type intermediate;
    Type total{};

    const int idx        = team_rank();
    const auto base_data = static_cast<Type*>(m_team_reduce);

    // Load values into the first not_greater_power_of_two values of the
    // reduction array in chunks. This means that only threads with an id in the
    // corresponding chunk load values and the reduction is always done by the
    // first not_greater_power_of_two threads.
    for (int start = 0; start < team_size();
         start += not_greater_power_of_two) {
      m_item.barrier(sycl::access::fence_space::local_space);
      if (idx >= start && idx < start + not_greater_power_of_two) {
        base_data[idx - start] = value;
      }
      m_item.barrier(sycl::access::fence_space::local_space);

      const Type partial_total =
          prescan(m_item, base_data, not_greater_power_of_two);
      if (idx >= start && idx < start + not_greater_power_of_two)
        intermediate = base_data[idx - start] + total;
      if (start == 0)
        total = partial_total;
      else
        total += partial_total;
    }

    if (global_accum) {
      if (team_size() == idx + 1) {
        base_data[team_size()] = atomic_fetch_add(global_accum, total);
      }
      m_item.barrier();  // Wait for atomic
      intermediate += base_data[team_size()];
    }

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
  KOKKOS_INLINE_FUNCTION static
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      vector_reduce(ReducerType const& reducer) {
    vector_reduce(reducer, reducer.reference());
  }

  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION static
      typename std::enable_if<is_reducer<ReducerType>::value>::type
      vector_reduce(ReducerType const& /*reducer*/,
                    typename ReducerType::value_type& /*value*/) {
    // FIXME_SYCL
    Kokkos::abort("Not implemented!");
  }

  //--------------------------------------------------------------------------
  /**\brief  Global reduction across all blocks
   *
   *  Return !0 if reducer contains the final value
   */
  template <typename ReducerType>
  KOKKOS_INLINE_FUNCTION static
      typename std::enable_if<is_reducer<ReducerType>::value, int>::type
      global_reduce(ReducerType const& /*reducer*/,
                    int* const /*global_scratch_flags*/,
                    void* const /*global_scratch_space*/, void* const /*shmem*/,
                    int const /*shmem_size*/) {
    // FIXME_SYCL
    Kokkos::abort("Not implemented!");
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
  // FIXME_SYCL Fix for vector_length>1.
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

  // FIXME_SYCL Fix for vector_length>1.
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

  // FIXME_SYCL Fix for vector_length>1.
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
  // FIXME_SYCL adapt for vector_length != 1
  for (iType i = loop_boundaries.start +
                 loop_boundaries.member.item().get_local_id(0);
       i < loop_boundaries.end;
       i += loop_boundaries.member.item().get_local_range(0))
    closure(i);
}

template <typename iType, class Closure, class ReducerType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<Kokkos::is_reducer<ReducerType>::value>::type
    parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                        iType, Impl::SYCLTeamMember>& loop_boundaries,
                    const Closure& closure, const ReducerType& reducer) {
  // FIXME_SYCL adapt for vector_length != 1
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

template <typename iType, class Closure, typename ValueType>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!Kokkos::is_reducer<ValueType>::value>::type
    parallel_reduce(const Impl::TeamVectorRangeBoundariesStruct<
                        iType, Impl::SYCLTeamMember>& loop_boundaries,
                    const Closure& closure, ValueType& result) {
  // FIXME_SYCL adapt for vector_length != 1
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
  // FIXME_SYC: adapt for vector_length!=1
  for (auto i = loop_boundaries.start; i != loop_boundaries.end; ++i)
    closure(i);
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
  // FIXME_SYCL adapt for vector_length != 1
  reducer.init(reducer.reference());

  for (iType i = loop_boundaries.start; i < loop_boundaries.end; ++i) {
    closure(i, reducer.reference());
  }
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
  // FIXME_SYCL adapt for vector_length != 1
  result = ValueType();

  for (iType i = loop_boundaries.start; i < loop_boundaries.end; ++i) {
    closure(i, result);
  }
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
  // FIXME_SYCL modify for vector_length!=1
  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, Closure>::value_type;

  value_type accum;
  reducer.init(accum);

  for (iType i = loop_boundaries.start; i < loop_boundaries.end; ++i) {
    closure(i, accum, true);
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
  if (single_struct.team_member.team_rank() == 0) lambda();
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::SYCLTeamMember>& single_struct,
    const FunctorType& lambda, ValueType& val) {
  if (single_struct.team_member.item().get_local_id(1) == 0) lambda(val);
}

template <class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::SYCLTeamMember>& single_struct,
    const FunctorType& lambda, ValueType& val) {
  if (single_struct.team_member.team_rank() == 0) lambda(val);
}

}  // namespace Kokkos

#endif

#endif /* #ifndef KOKKOS_SYCL_TEAM_HPP */
