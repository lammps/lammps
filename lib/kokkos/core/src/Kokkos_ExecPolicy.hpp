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

#ifndef KOKKOS_EXECPOLICY_HPP
#define KOKKOS_EXECPOLICY_HPP

#include <Kokkos_Core_fwd.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Tags.hpp>
#include <impl/Kokkos_AnalyzePolicy.hpp>
#include <Kokkos_Concepts.hpp>
#include <iostream>
#if defined(KOKKOS_ENABLE_PROFILING)
#include <typeinfo>
#endif  // KOKKOS_ENABLE_PROFILING

//----------------------------------------------------------------------------

namespace Kokkos {

struct ParallelForTag {};
struct ParallelScanTag {};
struct ParallelReduceTag {};

struct ChunkSize {
  int value;
  ChunkSize(int value_) : value(value_) {}
};

/** \brief  Execution policy for work over a range of an integral type.
 *
 * Valid template argument options:
 *
 *  With a specified execution space:
 *    < ExecSpace , WorkTag , { IntConst | IntType } >
 *    < ExecSpace , WorkTag , void >
 *    < ExecSpace , { IntConst | IntType } , void >
 *    < ExecSpace , void , void >
 *
 *  With the default execution space:
 *    < WorkTag , { IntConst | IntType } , void >
 *    < WorkTag , void , void >
 *    < { IntConst | IntType } , void , void >
 *    < void , void , void >
 *
 *  IntType  is a fundamental integral type
 *  IntConst is an Impl::integral_constant< IntType , Blocking >
 *
 *  Blocking is the granularity of partitioning the range among threads.
 */
template <class... Properties>
class RangePolicy : public Impl::PolicyTraits<Properties...> {
 public:
  typedef Impl::PolicyTraits<Properties...> traits;

 private:
  typename traits::execution_space m_space;
  typename traits::index_type m_begin;
  typename traits::index_type m_end;
  typename traits::index_type m_granularity;
  typename traits::index_type m_granularity_mask;

  template <class... OtherProperties>
  friend class RangePolicy;

 public:
  //! Tag this class as an execution policy
  typedef RangePolicy execution_policy;
  typedef typename traits::index_type member_type;
  typedef typename traits::index_type index_type;

  KOKKOS_INLINE_FUNCTION const typename traits::execution_space& space() const {
    return m_space;
  }
  KOKKOS_INLINE_FUNCTION member_type begin() const { return m_begin; }
  KOKKOS_INLINE_FUNCTION member_type end() const { return m_end; }

  // TODO: find a better workaround for Clangs weird instantiation order
  // This thing is here because of an instantiation error, where the RangePolicy
  // is inserted into FunctorValue Traits, which tries decltype on the operator.
  // It tries to do this even though the first argument of parallel for clearly
  // doesn't match.
  void operator()(const int&) const {}

  RangePolicy(const RangePolicy&) = default;
  RangePolicy(RangePolicy&&)      = default;

  template <class... OtherProperties>
  RangePolicy(const RangePolicy<OtherProperties...> p) {
    m_space            = p.m_space;
    m_begin            = p.m_begin;
    m_end              = p.m_end;
    m_granularity      = p.m_granularity;
    m_granularity_mask = p.m_granularity_mask;
  }

  inline RangePolicy() : m_space(), m_begin(0), m_end(0) {}

  /** \brief  Total range */
  inline RangePolicy(const typename traits::execution_space& work_space,
                     const member_type work_begin, const member_type work_end)
      : m_space(work_space),
        m_begin(work_begin < work_end ? work_begin : 0),
        m_end(work_begin < work_end ? work_end : 0),
        m_granularity(0),
        m_granularity_mask(0) {
    set_auto_chunk_size();
  }

  /** \brief  Total range */
  inline RangePolicy(const member_type work_begin, const member_type work_end)
      : RangePolicy(typename traits::execution_space(), work_begin, work_end) {
    set_auto_chunk_size();
  }

  /** \brief  Total range */
  template <class... Args>
  inline RangePolicy(const typename traits::execution_space& work_space,
                     const member_type work_begin, const member_type work_end,
                     Args... args)
      : m_space(work_space),
        m_begin(work_begin < work_end ? work_begin : 0),
        m_end(work_begin < work_end ? work_end : 0),
        m_granularity(0),
        m_granularity_mask(0) {
    set_auto_chunk_size();
    set(args...);
  }

  /** \brief  Total range */
  template <class... Args>
  inline RangePolicy(const member_type work_begin, const member_type work_end,
                     Args... args)
      : RangePolicy(typename traits::execution_space(), work_begin, work_end) {
    set_auto_chunk_size();
    set(args...);
  }

 private:
  inline void set() {}

 public:
  template <class... Args>
  inline void set(Args...) {
    static_assert(
        0 == sizeof...(Args),
        "Kokkos::RangePolicy: unhandled constructor arguments encountered.");
  }

  template <class... Args>
  inline void set(const ChunkSize& chunksize, Args... args) {
    m_granularity      = chunksize.value;
    m_granularity_mask = m_granularity - 1;
  }

 public:
  /** \brief return chunk_size */
  inline member_type chunk_size() const { return m_granularity; }

  /** \brief set chunk_size to a discrete value*/
  inline RangePolicy set_chunk_size(int chunk_size_) const {
    RangePolicy p        = *this;
    p.m_granularity      = chunk_size_;
    p.m_granularity_mask = p.m_granularity - 1;
    return p;
  }

 private:
  /** \brief finalize chunk_size if it was set to AUTO*/
  inline void set_auto_chunk_size() {
    typename traits::index_type concurrency =
        traits::execution_space::concurrency();
    if (concurrency == 0) concurrency = 1;

    if (m_granularity > 0) {
      if (!Impl::is_integral_power_of_two(m_granularity))
        Kokkos::abort("RangePolicy blocking granularity must be power of two");
    }

    member_type new_chunk_size = 1;
    while (new_chunk_size * 100 * concurrency < m_end - m_begin)
      new_chunk_size *= 2;
    if (new_chunk_size < 128) {
      new_chunk_size = 1;
      while ((new_chunk_size * 40 * concurrency < m_end - m_begin) &&
             (new_chunk_size < 128))
        new_chunk_size *= 2;
    }
    m_granularity      = new_chunk_size;
    m_granularity_mask = m_granularity - 1;
  }

 public:
  /** \brief  Subrange for a partition's rank and size.
   *
   *  Typically used to partition a range over a group of threads.
   */
  struct WorkRange {
    typedef typename RangePolicy::work_tag work_tag;
    typedef typename RangePolicy::member_type member_type;

    KOKKOS_INLINE_FUNCTION member_type begin() const { return m_begin; }
    KOKKOS_INLINE_FUNCTION member_type end() const { return m_end; }

    /** \brief  Subrange for a partition's rank and size.
     *
     *  Typically used to partition a range over a group of threads.
     */
    KOKKOS_INLINE_FUNCTION
    WorkRange(const RangePolicy& range, const int part_rank,
              const int part_size)
        : m_begin(0), m_end(0) {
      if (part_size) {
        // Split evenly among partitions, then round up to the granularity.
        const member_type work_part =
            ((((range.end() - range.begin()) + (part_size - 1)) / part_size) +
             range.m_granularity_mask) &
            ~member_type(range.m_granularity_mask);

        m_begin = range.begin() + work_part * part_rank;
        m_end   = m_begin + work_part;

        if (range.end() < m_begin) m_begin = range.end();
        if (range.end() < m_end) m_end = range.end();
      }
    }

   private:
    member_type m_begin;
    member_type m_end;
    WorkRange();
    WorkRange& operator=(const WorkRange&);
  };
};

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

namespace Impl {

template <class ExecSpace, class... Properties>
class TeamPolicyInternal : public Impl::PolicyTraits<Properties...> {
 private:
  typedef Impl::PolicyTraits<Properties...> traits;

 public:
  typedef typename traits::index_type index_type;

  //----------------------------------------
  /** \brief  Query maximum team size for a given functor.
   *
   *  This size takes into account execution space concurrency limitations and
   *  scratch memory space limitations for reductions, team reduce/scan, and
   *  team shared memory.
   *
   *  This function only works for single-operator functors.
   *  With multi-operator functors it cannot be determined
   *  which operator will be called.
   */
  template <class FunctorType>
  static int team_size_max(const FunctorType&);

  /** \brief  Query recommended team size for a given functor.
   *
   *  This size takes into account execution space concurrency limitations and
   *  scratch memory space limitations for reductions, team reduce/scan, and
   *  team shared memory.
   *
   *  This function only works for single-operator functors.
   *  With multi-operator functors it cannot be determined
   *  which operator will be called.
   */
  template <class FunctorType>
  static int team_size_recommended(const FunctorType&);

  template <class FunctorType>
  static int team_size_recommended(const FunctorType&, const int&);

  template <class FunctorType>
  int team_size_recommended(const FunctorType& functor,
                            const int vector_length);

  //----------------------------------------
  /** \brief  Construct policy with the given instance of the execution space */
  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request, int team_size_request,
                     int vector_length_request = 1);

  TeamPolicyInternal(const typename traits::execution_space&,
                     int league_size_request, const Kokkos::AUTO_t&,
                     int vector_length_request = 1);

  /** \brief  Construct policy with the default instance of the execution space
   */
  TeamPolicyInternal(int league_size_request, int team_size_request,
                     int vector_length_request = 1);

  TeamPolicyInternal(int league_size_request, const Kokkos::AUTO_t&,
                     int vector_length_request = 1);

  /*  TeamPolicyInternal( int league_size_request , int team_size_request );

    TeamPolicyInternal( int league_size_request , const Kokkos::AUTO_t & );*/

  /** \brief  The actual league size (number of teams) of the policy.
   *
   *  This may be smaller than the requested league size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int league_size() const;

  /** \brief  The actual team size (number of threads per team) of the policy.
   *
   *  This may be smaller than the requested team size due to limitations
   *  of the execution space.
   */
  KOKKOS_INLINE_FUNCTION int team_size() const;

  inline typename traits::index_type chunk_size() const;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  inline TeamPolicyInternal set_chunk_size(int chunk_size) const;
#else
  inline TeamPolicyInternal& set_chunk_size(int chunk_size);
#endif

  /** \brief  Parallel execution of a functor calls the functor once with
   *          each member of the execution policy.
   */
  struct member_type {
    /** \brief  Handle to the currently executing team shared scratch memory */
    KOKKOS_INLINE_FUNCTION
    typename traits::execution_space::scratch_memory_space team_shmem() const;

    /** \brief  Rank of this team within the league of teams */
    KOKKOS_INLINE_FUNCTION int league_rank() const;

    /** \brief  Number of teams in the league */
    KOKKOS_INLINE_FUNCTION int league_size() const;

    /** \brief  Rank of this thread within this team */
    KOKKOS_INLINE_FUNCTION int team_rank() const;

    /** \brief  Number of threads in this team */
    KOKKOS_INLINE_FUNCTION int team_size() const;

    /** \brief  Barrier among the threads of this team */
    KOKKOS_INLINE_FUNCTION void team_barrier() const;

    /** \brief  Intra-team reduction. Returns join of all values of the team
     * members. */
    template <class JoinOp>
    KOKKOS_INLINE_FUNCTION typename JoinOp::value_type team_reduce(
        const typename JoinOp::value_type, const JoinOp&) const;

    /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
     *
     *  The highest rank thread can compute the reduction total as
     *    reduction_total = dev.team_scan( value ) + value ;
     */
    template <typename Type>
    KOKKOS_INLINE_FUNCTION Type team_scan(const Type& value) const;

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
                                          Type* const global_accum) const;
  };
};

struct PerTeamValue {
  int value;
  PerTeamValue(int arg);
};

struct PerThreadValue {
  int value;
  PerThreadValue(int arg);
};

template <class iType, class... Args>
struct ExtractVectorLength {
  static inline iType value(
      typename std::enable_if<std::is_integral<iType>::value, iType>::type val,
      Args...) {
    return val;
  }
  static inline
      typename std::enable_if<!std::is_integral<iType>::value, int>::type
      value(
          typename std::enable_if<!std::is_integral<iType>::value, iType>::type,
          Args...) {
    return 1;
  }
};

template <class iType, class... Args>
inline typename std::enable_if<std::is_integral<iType>::value, iType>::type
extract_vector_length(iType val, Args...) {
  return val;
}

template <class iType, class... Args>
inline typename std::enable_if<!std::is_integral<iType>::value, int>::type
extract_vector_length(iType, Args...) {
  return 1;
}

}  // namespace Impl

Impl::PerTeamValue PerTeam(const int& arg);
Impl::PerThreadValue PerThread(const int& arg);

struct ScratchRequest {
  int level;

  int per_team;
  int per_thread;

  inline ScratchRequest(const int& level_,
                        const Impl::PerTeamValue& team_value) {
    level      = level_;
    per_team   = team_value.value;
    per_thread = 0;
  }

  inline ScratchRequest(const int& level_,
                        const Impl::PerThreadValue& thread_value) {
    level      = level_;
    per_team   = 0;
    per_thread = thread_value.value;
    ;
  }

  inline ScratchRequest(const int& level_, const Impl::PerTeamValue& team_value,
                        const Impl::PerThreadValue& thread_value) {
    level      = level_;
    per_team   = team_value.value;
    per_thread = thread_value.value;
    ;
  }

  inline ScratchRequest(const int& level_,
                        const Impl::PerThreadValue& thread_value,
                        const Impl::PerTeamValue& team_value) {
    level      = level_;
    per_team   = team_value.value;
    per_thread = thread_value.value;
    ;
  }
};

/** \brief  Execution policy for parallel work over a league of teams of
 * threads.
 *
 *  The work functor is called for each thread of each team such that
 *  the team's member threads are guaranteed to be concurrent.
 *
 *  The team's threads have access to team shared scratch memory and
 *  team collective operations.
 *
 *  If the WorkTag is non-void then the first calling argument of the
 *  work functor's parentheses operator is 'const WorkTag &'.
 *  This allows a functor to have multiple work member functions.
 *
 *  Order of template arguments does not matter, since the implementation
 *  uses variadic templates. Each and any of the template arguments can
 *  be omitted.
 *
 *  Possible Template arguments and their default values:
 *    ExecutionSpace (DefaultExecutionSpace): where to execute code. Must be
 * enabled. WorkTag (none): Tag which is used as the first argument for the
 * functor operator. Schedule<Type> (Schedule<Static>): Scheduling Policy
 * (Dynamic, or Static). IndexType<Type> (IndexType<ExecutionSpace::size_type>:
 * Integer Index type used to iterate over the Index space.
 *    LaunchBounds<unsigned,unsigned> Launch Bounds for CUDA compilation,
 *    default of LaunchBounds<0,0> indicates no launch bounds specified.
 */
template <class... Properties>
class TeamPolicy
    : public Impl::TeamPolicyInternal<
          typename Impl::PolicyTraits<Properties...>::execution_space,
          Properties...> {
  typedef Impl::TeamPolicyInternal<
      typename Impl::PolicyTraits<Properties...>::execution_space,
      Properties...>
      internal_policy;

  template <class... OtherProperties>
  friend class TeamPolicy;

 public:
  typedef Impl::PolicyTraits<Properties...> traits;

  typedef TeamPolicy execution_policy;

  TeamPolicy& operator=(const TeamPolicy&) = default;

  /** \brief  Construct policy with the given instance of the execution space */
  TeamPolicy(const typename traits::execution_space& space_,
             int league_size_request, int team_size_request,
             int vector_length_request = 1)
      : internal_policy(space_, league_size_request, team_size_request,
                        vector_length_request) {
    first_arg = false;
  }

  TeamPolicy(const typename traits::execution_space& space_,
             int league_size_request, const Kokkos::AUTO_t&,
             int vector_length_request = 1)
      : internal_policy(space_, league_size_request, Kokkos::AUTO(),
                        vector_length_request) {
    first_arg = false;
  }

  /** \brief  Construct policy with the default instance of the execution space
   */
  TeamPolicy(int league_size_request, int team_size_request,
             int vector_length_request = 1)
      : internal_policy(league_size_request, team_size_request,
                        vector_length_request) {
    first_arg = false;
  }

  TeamPolicy(int league_size_request, const Kokkos::AUTO_t&,
             int vector_length_request = 1)
      : internal_policy(league_size_request, Kokkos::AUTO(),
                        vector_length_request) {
    first_arg = false;
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  /** \brief  Construct policy with the given instance of the execution space */
  template <class... Args>
  TeamPolicy(const typename traits::execution_space&, int league_size_request,
             int team_size_request, int vector_length_request, Args... args)
      : internal_policy(typename traits::execution_space(), league_size_request,
                        team_size_request, vector_length_request) {
    first_arg = false;
    set(args...);
  }

  template <class... Args>
  TeamPolicy(const typename traits::execution_space&, int league_size_request,
             const Kokkos::AUTO_t&, int vector_length_request, Args... args)
      : internal_policy(typename traits::execution_space(), league_size_request,
                        Kokkos::AUTO(), vector_length_request) {
    first_arg = false;
    set(args...);
  }

  /** \brief  Construct policy with the default instance of the execution space
   */
  template <class... Args>
  TeamPolicy(int league_size_request, int team_size_request,
             int vector_length_request, Args... args)
      : internal_policy(league_size_request, team_size_request,
                        vector_length_request) {
    first_arg = false;
    set(args...);
  }

  template <class... Args>
  TeamPolicy(int league_size_request, const Kokkos::AUTO_t&,
             int vector_length_request, Args... args)
      : internal_policy(league_size_request, Kokkos::AUTO(),
                        vector_length_request) {
    first_arg = false;
    set(args...);
  }

  /** \brief  Construct policy with the given instance of the execution space */
  template <class... Args>
  TeamPolicy(const typename traits::execution_space&, int league_size_request,
             int team_size_request, Args... args)
      : internal_policy(typename traits::execution_space(), league_size_request,
                        team_size_request,
                        Kokkos::Impl::extract_vector_length<Args...>(args...)) {
    first_arg = true;
    set(args...);
  }

  template <class... Args>
  TeamPolicy(const typename traits::execution_space&, int league_size_request,
             const Kokkos::AUTO_t&, Args... args)
      : internal_policy(typename traits::execution_space(), league_size_request,
                        Kokkos::AUTO(),
                        Kokkos::Impl::extract_vector_length<Args...>(args...)) {
    first_arg = true;
    set(args...);
  }

  /** \brief  Construct policy with the default instance of the execution space
   */
  template <class... Args>
  TeamPolicy(int league_size_request, int team_size_request, Args... args)
      : internal_policy(league_size_request, team_size_request,
                        Kokkos::Impl::extract_vector_length<Args...>(args...)) {
    first_arg = true;
    set(args...);
  }

  template <class... Args>
  TeamPolicy(int league_size_request, const Kokkos::AUTO_t&, Args... args)
      : internal_policy(league_size_request, Kokkos::AUTO(),
                        Kokkos::Impl::extract_vector_length<Args...>(args...)) {
    first_arg = true;
    set(args...);
  }
#endif

  template <class... OtherProperties>
  TeamPolicy(const TeamPolicy<OtherProperties...> p) : internal_policy(p) {
    first_arg = p.first_arg;
  }

 private:
  bool first_arg;
  TeamPolicy(const internal_policy& p) : internal_policy(p) {
    first_arg = false;
  }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  inline void set() {}
#endif

 public:
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  template <class... Args>
  inline void set(Args...) {
    static_assert(
        0 == sizeof...(Args),
        "Kokkos::TeamPolicy: unhandled constructor arguments encountered.");
  }

  template <class iType, class... Args>
  inline typename std::enable_if<std::is_integral<iType>::value>::type set(
      iType, Args... args) {
    if (first_arg) {
      first_arg = false;
      set(args...);
    } else {
      first_arg = false;
      Kokkos::Impl::throw_runtime_exception(
          "Kokkos::TeamPolicy: integer argument to constructor in illegal "
          "place.");
    }
  }

  template <class... Args>
  inline void set(const ChunkSize& chunksize, Args... args) {
    first_arg = false;
    internal_policy::internal_set_chunk_size(chunksize.value);
    set(args...);
  }

  template <class... Args>
  inline void set(const ScratchRequest& scr_request, Args... args) {
    first_arg = false;
    internal_policy::internal_set_scratch_size(
        scr_request.level, Impl::PerTeamValue(scr_request.per_team),
        Impl::PerThreadValue(scr_request.per_thread));
    set(args...);
  }

  inline TeamPolicy set_chunk_size(int chunk) const {
    return TeamPolicy(internal_policy::set_chunk_size(chunk));
  }

  inline TeamPolicy set_scratch_size(const int& level,
                                     const Impl::PerTeamValue& per_team) const {
    return TeamPolicy(internal_policy::set_scratch_size(level, per_team));
  }
  inline TeamPolicy set_scratch_size(
      const int& level, const Impl::PerThreadValue& per_thread) const {
    return TeamPolicy(internal_policy::set_scratch_size(level, per_thread));
  }
  inline TeamPolicy set_scratch_size(
      const int& level, const Impl::PerTeamValue& per_team,
      const Impl::PerThreadValue& per_thread) const {
    return TeamPolicy(
        internal_policy::set_scratch_size(level, per_team, per_thread));
  }
  inline TeamPolicy set_scratch_size(const int& level,
                                     const Impl::PerThreadValue& per_thread,
                                     const Impl::PerTeamValue& per_team) const {
    return TeamPolicy(
        internal_policy::set_scratch_size(level, per_team, per_thread));
  }

#else
  inline TeamPolicy& set_chunk_size(int chunk) {
    static_assert(std::is_same<decltype(internal_policy::set_chunk_size(chunk)),
                               internal_policy&>::value,
                  "internal set_chunk_size should return a reference");
    return static_cast<TeamPolicy&>(internal_policy::set_chunk_size(chunk));
  }

  inline TeamPolicy& set_scratch_size(const int& level,
                                      const Impl::PerTeamValue& per_team) {
    static_assert(std::is_same<decltype(internal_policy::set_scratch_size(
                                   level, per_team)),
                               internal_policy&>::value,
                  "internal set_chunk_size should return a reference");
    return static_cast<TeamPolicy&>(
        internal_policy::set_scratch_size(level, per_team));
  }
  inline TeamPolicy& set_scratch_size(const int& level,
                                      const Impl::PerThreadValue& per_thread) {
    return static_cast<TeamPolicy&>(
        internal_policy::set_scratch_size(level, per_thread));
  }
  inline TeamPolicy& set_scratch_size(const int& level,
                                      const Impl::PerTeamValue& per_team,
                                      const Impl::PerThreadValue& per_thread) {
    return static_cast<TeamPolicy&>(
        internal_policy::set_scratch_size(level, per_team, per_thread));
  }
  inline TeamPolicy& set_scratch_size(const int& level,
                                      const Impl::PerThreadValue& per_thread,
                                      const Impl::PerTeamValue& per_team) {
    return static_cast<TeamPolicy&>(
        internal_policy::set_scratch_size(level, per_team, per_thread));
  }
#endif
};

namespace Impl {

template <typename iType, class TeamMemberType>
struct TeamThreadRangeBoundariesStruct {
 private:
  KOKKOS_INLINE_FUNCTION static iType ibegin(const iType& arg_begin,
                                             const iType& arg_end,
                                             const iType& arg_rank,
                                             const iType& arg_size) {
    return arg_begin +
           ((arg_end - arg_begin + arg_size - 1) / arg_size) * arg_rank;
  }

  KOKKOS_INLINE_FUNCTION static iType iend(const iType& arg_begin,
                                           const iType& arg_end,
                                           const iType& arg_rank,
                                           const iType& arg_size) {
    const iType end_ =
        arg_begin +
        ((arg_end - arg_begin + arg_size - 1) / arg_size) * (arg_rank + 1);
    return end_ < arg_end ? end_ : arg_end;
  }

 public:
  typedef iType index_type;
  const iType start;
  const iType end;
  enum { increment = 1 };
  const TeamMemberType& thread;

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const TeamMemberType& arg_thread,
                                  const iType& arg_end)
      : start(
            ibegin(0, arg_end, arg_thread.team_rank(), arg_thread.team_size())),
        end(iend(0, arg_end, arg_thread.team_rank(), arg_thread.team_size())),
        thread(arg_thread) {}

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct(const TeamMemberType& arg_thread,
                                  const iType& arg_begin, const iType& arg_end)
      : start(ibegin(arg_begin, arg_end, arg_thread.team_rank(),
                     arg_thread.team_size())),
        end(iend(arg_begin, arg_end, arg_thread.team_rank(),
                 arg_thread.team_size())),
        thread(arg_thread) {}
};

template <typename iType, class TeamMemberType>
struct TeamVectorRangeBoundariesStruct {
 private:
  KOKKOS_INLINE_FUNCTION static iType ibegin(const iType& arg_begin,
                                             const iType& arg_end,
                                             const iType& arg_rank,
                                             const iType& arg_size) {
    return arg_begin +
           ((arg_end - arg_begin + arg_size - 1) / arg_size) * arg_rank;
  }

  KOKKOS_INLINE_FUNCTION static iType iend(const iType& arg_begin,
                                           const iType& arg_end,
                                           const iType& arg_rank,
                                           const iType& arg_size) {
    const iType end_ =
        arg_begin +
        ((arg_end - arg_begin + arg_size - 1) / arg_size) * (arg_rank + 1);
    return end_ < arg_end ? end_ : arg_end;
  }

 public:
  typedef iType index_type;
  const iType start;
  const iType end;
  enum { increment = 1 };
  const TeamMemberType& thread;

  KOKKOS_INLINE_FUNCTION
  TeamVectorRangeBoundariesStruct(const TeamMemberType& arg_thread,
                                  const iType& arg_end)
      : start(
            ibegin(0, arg_end, arg_thread.team_rank(), arg_thread.team_size())),
        end(iend(0, arg_end, arg_thread.team_rank(), arg_thread.team_size())),
        thread(arg_thread) {}

  KOKKOS_INLINE_FUNCTION
  TeamVectorRangeBoundariesStruct(const TeamMemberType& arg_thread,
                                  const iType& arg_begin, const iType& arg_end)
      : start(ibegin(arg_begin, arg_end, arg_thread.team_rank(),
                     arg_thread.team_size())),
        end(iend(arg_begin, arg_end, arg_thread.team_rank(),
                 arg_thread.team_size())),
        thread(arg_thread) {}
};

template <typename iType, class TeamMemberType>
struct ThreadVectorRangeBoundariesStruct {
  typedef iType index_type;
  const index_type start;
  const index_type end;
  enum { increment = 1 };

  KOKKOS_INLINE_FUNCTION
  constexpr ThreadVectorRangeBoundariesStruct(const TeamMemberType,
                                              const index_type& count) noexcept
      : start(static_cast<index_type>(0)), end(count) {}

  KOKKOS_INLINE_FUNCTION
  constexpr ThreadVectorRangeBoundariesStruct(const index_type& count) noexcept
      : start(static_cast<index_type>(0)), end(count) {}

  KOKKOS_INLINE_FUNCTION
  constexpr ThreadVectorRangeBoundariesStruct(
      const TeamMemberType, const index_type& arg_begin,
      const index_type& arg_end) noexcept
      : start(static_cast<index_type>(arg_begin)), end(arg_end) {}

  KOKKOS_INLINE_FUNCTION
  constexpr ThreadVectorRangeBoundariesStruct(
      const index_type& arg_begin, const index_type& arg_end) noexcept
      : start(static_cast<index_type>(arg_begin)), end(arg_end) {}
};

template <class TeamMemberType>
struct ThreadSingleStruct {
  const TeamMemberType& team_member;
  KOKKOS_INLINE_FUNCTION
  ThreadSingleStruct(const TeamMemberType& team_member_)
      : team_member(team_member_) {}
};

template <class TeamMemberType>
struct VectorSingleStruct {
  const TeamMemberType& team_member;
  KOKKOS_INLINE_FUNCTION
  VectorSingleStruct(const TeamMemberType& team_member_)
      : team_member(team_member_) {}
};

}  // namespace Impl

/** \brief  Execution policy for parallel work over a threads within a team.
 *
 *  The range is split over all threads in a team. The Mapping scheme depends on
 * the architecture. This policy is used together with a parallel pattern as a
 * nested layer within a kernel launched with the TeamPolicy. This variant
 * expects a single count. So the range is (0,count].
 */
template <typename iType, class TeamMemberType, class _never_use_this_overload>
KOKKOS_INLINE_FUNCTION_DELETED
    Impl::TeamThreadRangeBoundariesStruct<iType, TeamMemberType>
    TeamThreadRange(const TeamMemberType&, const iType& count) = delete;

/** \brief  Execution policy for parallel work over a threads within a team.
 *
 *  The range is split over all threads in a team. The Mapping scheme depends on
 * the architecture. This policy is used together with a parallel pattern as a
 * nested layer within a kernel launched with the TeamPolicy. This variant
 * expects a begin and end. So the range is (begin,end].
 */
template <typename iType1, typename iType2, class TeamMemberType,
          class _never_use_this_overload>
KOKKOS_INLINE_FUNCTION_DELETED Impl::TeamThreadRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, TeamMemberType>
TeamThreadRange(const TeamMemberType&, const iType1& begin,
                const iType2& end) = delete;

/** \brief  Execution policy for parallel work over a threads within a team.
 *
 *  The range is split over all threads in a team. The Mapping scheme depends on
 * the architecture. This policy is used together with a parallel pattern as a
 * nested layer within a kernel launched with the TeamPolicy. This variant
 * expects a single count. So the range is (0,count].
 */
template <typename iType, class TeamMemberType, class _never_use_this_overload>
KOKKOS_INLINE_FUNCTION_DELETED
    Impl::TeamThreadRangeBoundariesStruct<iType, TeamMemberType>
    TeamVectorRange(const TeamMemberType&, const iType& count) = delete;

/** \brief  Execution policy for parallel work over a threads within a team.
 *
 *  The range is split over all threads in a team. The Mapping scheme depends on
 * the architecture. This policy is used together with a parallel pattern as a
 * nested layer within a kernel launched with the TeamPolicy. This variant
 * expects a begin and end. So the range is (begin,end].
 */
template <typename iType1, typename iType2, class TeamMemberType,
          class _never_use_this_overload>
KOKKOS_INLINE_FUNCTION_DELETED Impl::TeamThreadRangeBoundariesStruct<
    typename std::common_type<iType1, iType2>::type, TeamMemberType>
TeamVectorRange(const TeamMemberType&, const iType1& begin,
                const iType2& end) = delete;

/** \brief  Execution policy for a vector parallel loop.
 *
 *  The range is split over all vector lanes in a thread. The Mapping scheme
 * depends on the architecture. This policy is used together with a parallel
 * pattern as a nested layer within a kernel launched with the TeamPolicy. This
 * variant expects a single count. So the range is (0,count].
 */
template <typename iType, class TeamMemberType, class _never_use_this_overload>
KOKKOS_INLINE_FUNCTION_DELETED
    Impl::ThreadVectorRangeBoundariesStruct<iType, TeamMemberType>
    ThreadVectorRange(const TeamMemberType&, const iType& count) = delete;

template <typename iType, class TeamMemberType, class _never_use_this_overload>
KOKKOS_INLINE_FUNCTION_DELETED
    Impl::ThreadVectorRangeBoundariesStruct<iType, TeamMemberType>
    ThreadVectorRange(const TeamMemberType&, const iType& arg_begin,
                      const iType& arg_end) = delete;

#if defined(KOKKOS_ENABLE_PROFILING)
namespace Impl {

template <typename FunctorType, typename TagType,
          bool HasTag = !std::is_same<TagType, void>::value>
struct ParallelConstructName;

template <typename FunctorType, typename TagType>
struct ParallelConstructName<FunctorType, TagType, true> {
  ParallelConstructName(std::string const& label) : label_ref(label) {
    if (label.empty()) {
      default_name = std::string(typeid(FunctorType).name()) + "/" +
                     typeid(TagType).name();
    }
  }
  std::string const& get() {
    return (label_ref.empty()) ? default_name : label_ref;
  }
  std::string const& label_ref;
  std::string default_name;
};

template <typename FunctorType, typename TagType>
struct ParallelConstructName<FunctorType, TagType, false> {
  ParallelConstructName(std::string const& label) : label_ref(label) {
    if (label.empty()) {
      default_name = std::string(typeid(FunctorType).name());
    }
  }
  std::string const& get() {
    return (label_ref.empty()) ? default_name : label_ref;
  }
  std::string const& label_ref;
  std::string default_name;
};

}  // namespace Impl
#endif /* defined KOKKOS_ENABLE_PROFILING */

}  // namespace Kokkos

namespace Kokkos {
namespace Experimental {

namespace Impl {
template <class Property, class Policy>
struct PolicyPropertyAdaptor;

template <unsigned long P, class... Properties>
struct PolicyPropertyAdaptor<WorkItemProperty::ImplWorkItemProperty<P>,
                             RangePolicy<Properties...>> {
  typedef RangePolicy<Properties...> policy_in_t;
  typedef RangePolicy<typename policy_in_t::traits::execution_space,
                      typename policy_in_t::traits::schedule_type,
                      typename policy_in_t::traits::work_tag,
                      typename policy_in_t::traits::index_type,
                      typename policy_in_t::traits::iteration_pattern,
                      typename policy_in_t::traits::launch_bounds,
                      WorkItemProperty::ImplWorkItemProperty<P>>
      policy_out_t;
};

template <unsigned long P, class... Properties>
struct PolicyPropertyAdaptor<WorkItemProperty::ImplWorkItemProperty<P>,
                             TeamPolicy<Properties...>> {
  typedef TeamPolicy<Properties...> policy_in_t;
  typedef TeamPolicy<typename policy_in_t::traits::execution_space,
                     typename policy_in_t::traits::schedule_type,
                     typename policy_in_t::traits::work_tag,
                     typename policy_in_t::traits::index_type,
                     typename policy_in_t::traits::iteration_pattern,
                     typename policy_in_t::traits::launch_bounds,
                     WorkItemProperty::ImplWorkItemProperty<P>>
      policy_out_t;
};
}  // namespace Impl

template <class PolicyType, unsigned long P>
constexpr typename Impl::PolicyPropertyAdaptor<
    WorkItemProperty::ImplWorkItemProperty<P>, PolicyType>::policy_out_t
require(const PolicyType p, WorkItemProperty::ImplWorkItemProperty<P>) {
  return typename Impl::PolicyPropertyAdaptor<
      WorkItemProperty::ImplWorkItemProperty<P>, PolicyType>::policy_out_t(p);
}
}  // namespace Experimental
}  // namespace Kokkos
#endif /* #define KOKKOS_EXECPOLICY_HPP */
