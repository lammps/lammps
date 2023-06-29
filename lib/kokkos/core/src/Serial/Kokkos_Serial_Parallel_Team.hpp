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

#ifndef KOKKO_SERIAL_PARALLEL_TEAM_HPP
#define KOKKO_SERIAL_PARALLEL_TEAM_HPP

#include <Kokkos_Parallel.hpp>

namespace Kokkos {
namespace Impl {

/*
 * < Kokkos::Serial , WorkArgTag >
 * < WorkArgTag , Impl::enable_if< std::is_same< Kokkos::Serial ,
 * Kokkos::DefaultExecutionSpace >::value >::type >
 *
 */
template <class... Properties>
class TeamPolicyInternal<Kokkos::Serial, Properties...>
    : public PolicyTraits<Properties...> {
 private:
  size_t m_team_scratch_size[2];
  size_t m_thread_scratch_size[2];
  int m_league_size;
  int m_chunk_size;

 public:
  //! Tag this class as a kokkos execution policy
  using execution_policy = TeamPolicyInternal;

  using traits = PolicyTraits<Properties...>;

  //! Execution space of this execution policy:
  using execution_space = Kokkos::Serial;

  const typename traits::execution_space& space() const {
    static typename traits::execution_space m_space;
    return m_space;
  }

  template <class ExecSpace, class... OtherProperties>
  friend class TeamPolicyInternal;

  template <class... OtherProperties>
  TeamPolicyInternal(
      const TeamPolicyInternal<Kokkos::Serial, OtherProperties...>& p) {
    m_league_size            = p.m_league_size;
    m_team_scratch_size[0]   = p.m_team_scratch_size[0];
    m_thread_scratch_size[0] = p.m_thread_scratch_size[0];
    m_team_scratch_size[1]   = p.m_team_scratch_size[1];
    m_thread_scratch_size[1] = p.m_thread_scratch_size[1];
    m_chunk_size             = p.m_chunk_size;
  }

  //----------------------------------------

  template <class FunctorType>
  int team_size_max(const FunctorType&, const ParallelForTag&) const {
    return 1;
  }
  template <class FunctorType>
  int team_size_max(const FunctorType&, const ParallelReduceTag&) const {
    return 1;
  }
  template <class FunctorType, class ReducerType>
  int team_size_max(const FunctorType&, const ReducerType&,
                    const ParallelReduceTag&) const {
    return 1;
  }
  template <class FunctorType>
  int team_size_recommended(const FunctorType&, const ParallelForTag&) const {
    return 1;
  }
  template <class FunctorType>
  int team_size_recommended(const FunctorType&,
                            const ParallelReduceTag&) const {
    return 1;
  }
  template <class FunctorType, class ReducerType>
  int team_size_recommended(const FunctorType&, const ReducerType&,
                            const ParallelReduceTag&) const {
    return 1;
  }

  //----------------------------------------

  inline int team_size() const { return 1; }
  inline bool impl_auto_team_size() const { return false; }
  inline bool impl_auto_vector_length() const { return false; }
  inline void impl_set_team_size(size_t) {}
  inline void impl_set_vector_length(size_t) {}
  inline int league_size() const { return m_league_size; }
  inline size_t scratch_size(const int& level, int = 0) const {
    return m_team_scratch_size[level] + m_thread_scratch_size[level];
  }

  inline int impl_vector_length() const { return 1; }
  inline static int vector_length_max() {
    return 1024;
  }  // Use arbitrary large number, is meant as a vectorizable length

  inline static int scratch_size_max(int level) {
    return (level == 0 ? 1024 * 32 : 20 * 1024 * 1024);
  }
  /** \brief  Specify league size, request team size */
  TeamPolicyInternal(const execution_space&, int league_size_request,
                     int team_size_request, int /* vector_length_request */ = 1)
      : m_team_scratch_size{0, 0},
        m_thread_scratch_size{0, 0},
        m_league_size(league_size_request),
        m_chunk_size(32) {
    if (team_size_request > 1)
      Kokkos::abort("Kokkos::abort: Requested Team Size is too large!");
  }

  TeamPolicyInternal(const execution_space& space, int league_size_request,
                     const Kokkos::AUTO_t& /**team_size_request*/,
                     int vector_length_request = 1)
      : TeamPolicyInternal(space, league_size_request, -1,
                           vector_length_request) {}

  TeamPolicyInternal(const execution_space& space, int league_size_request,
                     const Kokkos::AUTO_t& /* team_size_request */
                     ,
                     const Kokkos::AUTO_t& /* vector_length_request */
                     )
      : TeamPolicyInternal(space, league_size_request, -1, -1) {}

  TeamPolicyInternal(const execution_space& space, int league_size_request,
                     int team_size_request,
                     const Kokkos::AUTO_t& /* vector_length_request */
                     )
      : TeamPolicyInternal(space, league_size_request, team_size_request, -1) {}

  TeamPolicyInternal(int league_size_request,
                     const Kokkos::AUTO_t& team_size_request,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(),
                           league_size_request, team_size_request,
                           vector_length_request) {}

  TeamPolicyInternal(int league_size_request,
                     const Kokkos::AUTO_t& team_size_request,
                     const Kokkos::AUTO_t& vector_length_request)
      : TeamPolicyInternal(typename traits::execution_space(),
                           league_size_request, team_size_request,
                           vector_length_request) {}
  TeamPolicyInternal(int league_size_request, int team_size_request,
                     const Kokkos::AUTO_t& vector_length_request)
      : TeamPolicyInternal(typename traits::execution_space(),
                           league_size_request, team_size_request,
                           vector_length_request) {}

  TeamPolicyInternal(int league_size_request, int team_size_request,
                     int vector_length_request = 1)
      : TeamPolicyInternal(typename traits::execution_space(),
                           league_size_request, team_size_request,
                           vector_length_request) {}

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

  using member_type = Impl::HostThreadTeamMember<Kokkos::Serial>;
};

template <class FunctorType, class... Properties>
class ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                  Kokkos::Serial> {
 private:
  enum { TEAM_REDUCE_SIZE = 512 };

  using Policy = TeamPolicyInternal<Kokkos::Serial, Properties...>;
  using Member = typename Policy::member_type;

  const FunctorType m_functor;
  const Policy m_policy;
  const int m_league;
  const size_t m_shared;

  template <class TagType>
  inline std::enable_if_t<std::is_void<TagType>::value> exec(
      HostThreadTeamData& data) const {
    for (int ileague = 0; ileague < m_league; ++ileague) {
      m_functor(Member(data, ileague, m_league));
    }
  }

  template <class TagType>
  inline std::enable_if_t<!std::is_void<TagType>::value> exec(
      HostThreadTeamData& data) const {
    const TagType t{};
    for (int ileague = 0; ileague < m_league; ++ileague) {
      m_functor(t, Member(data, ileague, m_league));
    }
  }

 public:
  inline void execute() const {
    const size_t pool_reduce_size  = 0;  // Never shrinks
    const size_t team_reduce_size  = TEAM_REDUCE_SIZE;
    const size_t team_shared_size  = m_shared;
    const size_t thread_local_size = 0;  // Never shrinks

    auto* internal_instance = m_policy.space().impl_internal_space_instance();
    // Need to lock resize_thread_team_data
    std::lock_guard<std::mutex> lock(
        internal_instance->m_thread_team_data_mutex);
    internal_instance->resize_thread_team_data(
        pool_reduce_size, team_reduce_size, team_shared_size,
        thread_local_size);

    this->template exec<typename Policy::work_tag>(
        internal_instance->m_thread_team_data);
  }

  ParallelFor(const FunctorType& arg_functor, const Policy& arg_policy)
      : m_functor(arg_functor),
        m_policy(arg_policy),
        m_league(arg_policy.league_size()),
        m_shared(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                 FunctorTeamShmemSize<FunctorType>::value(arg_functor, 1)) {}
};

/*--------------------------------------------------------------------------*/

template <class CombinedFunctorReducerType, class... Properties>
class ParallelReduce<CombinedFunctorReducerType,
                     Kokkos::TeamPolicy<Properties...>, Kokkos::Serial> {
 private:
  static constexpr int TEAM_REDUCE_SIZE = 512;

  using Policy      = TeamPolicyInternal<Kokkos::Serial, Properties...>;
  using FunctorType = typename CombinedFunctorReducerType::functor_type;
  using ReducerType = typename CombinedFunctorReducerType::reducer_type;

  using Member  = typename Policy::member_type;
  using WorkTag = typename Policy::work_tag;

  using pointer_type   = typename ReducerType::pointer_type;
  using reference_type = typename ReducerType::reference_type;

  const CombinedFunctorReducerType m_functor_reducer;
  const Policy m_policy;
  const int m_league;
  pointer_type m_result_ptr;
  size_t m_shared;

  template <class TagType>
  inline std::enable_if_t<std::is_void<TagType>::value> exec(
      HostThreadTeamData& data, reference_type update) const {
    for (int ileague = 0; ileague < m_league; ++ileague) {
      m_functor_reducer.get_functor()(Member(data, ileague, m_league), update);
    }
  }

  template <class TagType>
  inline std::enable_if_t<!std::is_void<TagType>::value> exec(
      HostThreadTeamData& data, reference_type update) const {
    const TagType t{};

    for (int ileague = 0; ileague < m_league; ++ileague) {
      m_functor_reducer.get_functor()(t, Member(data, ileague, m_league),
                                      update);
    }
  }

 public:
  inline void execute() const {
    const size_t pool_reduce_size =
        m_functor_reducer.get_reducer().value_size();

    const size_t team_reduce_size  = TEAM_REDUCE_SIZE;
    const size_t team_shared_size  = m_shared;
    const size_t thread_local_size = 0;  // Never shrinks

    auto* internal_instance = m_policy.space().impl_internal_space_instance();
    // Need to lock resize_thread_team_data
    std::lock_guard<std::mutex> lock(
        internal_instance->m_thread_team_data_mutex);
    internal_instance->resize_thread_team_data(
        pool_reduce_size, team_reduce_size, team_shared_size,
        thread_local_size);

    pointer_type ptr =
        m_result_ptr
            ? m_result_ptr
            : pointer_type(
                  internal_instance->m_thread_team_data.pool_reduce_local());

    reference_type update = m_functor_reducer.get_reducer().init(ptr);

    this->template exec<WorkTag>(internal_instance->m_thread_team_data, update);

    m_functor_reducer.get_reducer().final(ptr);
  }

  template <class ViewType>
  ParallelReduce(const CombinedFunctorReducerType& arg_functor_reducer,
                 const Policy& arg_policy, const ViewType& arg_result)
      : m_functor_reducer(arg_functor_reducer),
        m_policy(arg_policy),
        m_league(arg_policy.league_size()),
        m_result_ptr(arg_result.data()),
        m_shared(arg_policy.scratch_size(0) + arg_policy.scratch_size(1) +
                 FunctorTeamShmemSize<FunctorType>::value(
                     m_functor_reducer.get_functor(), 1)) {
    static_assert(Kokkos::is_view<ViewType>::value,
                  "Reduction result on Kokkos::Serial must be a Kokkos::View");

    static_assert(
        Kokkos::Impl::MemorySpaceAccess<typename ViewType::memory_space,
                                        Kokkos::HostSpace>::accessible,
        "Kokkos::Serial reduce result must be a View accessible from "
        "HostSpace");
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
