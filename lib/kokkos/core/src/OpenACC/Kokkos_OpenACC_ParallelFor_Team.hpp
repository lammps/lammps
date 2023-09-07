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

#ifndef KOKKOS_OPENACC_PARALLEL_FOR_TEAM_HPP
#define KOKKOS_OPENACC_PARALLEL_FOR_TEAM_HPP

#include <OpenACC/Kokkos_OpenACC_Team.hpp>
#include <OpenACC/Kokkos_OpenACC_FunctorAdapter.hpp>

#ifdef KOKKOS_ENABLE_OPENACC_COLLAPSE_HIERARCHICAL_CONSTRUCTS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Hierarchical Parallelism -> Team level implementation
template <class FunctorType, class... Properties>
class Kokkos::Impl::ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                                Kokkos::Experimental::OpenACC> {
 private:
  using Policy = Kokkos::Impl::TeamPolicyInternal<Kokkos::Experimental::OpenACC,
                                                  Properties...>;
  Kokkos::Experimental::Impl::FunctorAdapter<
      FunctorType, Policy, Kokkos::Experimental::Impl::RoutineClause::seq>
      m_functor;
  using Member = typename Policy::member_type;

  const Policy m_policy;

 public:
  inline void execute() const {
    auto league_size   = m_policy.league_size();
    auto team_size     = m_policy.team_size();
    auto vector_length = m_policy.impl_vector_length();

    auto const a_functor(m_functor);

#pragma acc parallel loop gang vector num_gangs(league_size) \
    vector_length(team_size* vector_length) copyin(a_functor)
    for (int i = 0; i < league_size * team_size * vector_length; i++) {
      int league_id = i / (team_size * vector_length);
      typename Policy::member_type team(league_id, league_size, team_size,
                                        vector_length);
      a_functor(team);
    }
  }

  inline ParallelFor(const FunctorType& arg_functor, Policy arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

namespace Kokkos {

// Hierarchical Parallelism -> Team thread level implementation
#pragma acc routine seq
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::OpenACCTeamMember>&
        loop_boundaries,
    const Lambda& lambda) {
  iType j_start =
      loop_boundaries.team.team_rank() / loop_boundaries.team.vector_length();
  iType j_end  = loop_boundaries.end;
  iType j_step = loop_boundaries.team.team_size();
  if (j_start >= loop_boundaries.start) {
#pragma acc loop seq
    for (iType j = j_start; j < j_end; j += j_step) {
      lambda(j);
    }
  }
}

// Hierarchical Parallelism -> Thread vector level implementation
#pragma acc routine seq
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::OpenACCTeamMember>& loop_boundaries,
    const Lambda& lambda) {
  iType j_start =
      loop_boundaries.team.team_rank() % loop_boundaries.team.vector_length();
  iType j_end  = loop_boundaries.end;
  iType j_step = loop_boundaries.team.vector_length();
  if (j_start >= loop_boundaries.start) {
#pragma acc loop seq
    for (iType j = j_start; j < j_end; j += j_step) {
      lambda(j);
    }
  }
}

// Hierarchical Parallelism -> Team vector level implementation
#pragma acc routine seq
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamVectorRangeBoundariesStruct<iType, Impl::OpenACCTeamMember>&
        loop_boundaries,
    const Lambda& lambda) {
  iType j_start =
      loop_boundaries.team.team_rank() % loop_boundaries.team.vector_length();
  iType j_end  = loop_boundaries.end;
  iType j_step = loop_boundaries.team.vector_length();
  if (j_start >= loop_boundaries.start) {
#pragma acc loop seq
    for (iType j = j_start; j < j_end; j += j_step) {
      lambda(j);
    }
  }
}

}  // namespace Kokkos

#else  // KOKKOS_ENABLE_OPENACC_COLLAPSE_HIERARCHICAL_CONSTRUCTS

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Hierarchical Parallelism -> Team level implementation
template <class FunctorType, class... Properties>
class Kokkos::Impl::ParallelFor<FunctorType, Kokkos::TeamPolicy<Properties...>,
                                Kokkos::Experimental::OpenACC> {
 private:
  using Policy = Kokkos::Impl::TeamPolicyInternal<Kokkos::Experimental::OpenACC,
                                                  Properties...>;
  Kokkos::Experimental::Impl::FunctorAdapter<FunctorType, Policy, worker>
      m_functor;
  using Member = typename Policy::member_type;

  const Policy m_policy;

 public:
  inline void execute() const {
    auto league_size   = m_policy.league_size();
    auto team_size     = m_policy.team_size();
    auto vector_length = m_policy.impl_vector_length();

    auto const a_functor(m_functor);

#pragma acc parallel loop gang num_gangs(league_size) num_workers(team_size) \
    vector_length(vector_length) copyin(a_functor)
    for (int i = 0; i < league_size; i++) {
      int league_id = i;
      typename Policy::member_type team(league_id, league_size, team_size,
                                        vector_length);
      a_functor(team);
    }
  }

  inline ParallelFor(const FunctorType& arg_functor, Policy arg_policy)
      : m_functor(arg_functor), m_policy(arg_policy) {}
};

namespace Kokkos {

// Hierarchical Parallelism -> Team thread level implementation
#pragma acc routine worker
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<iType, Impl::OpenACCTeamMember>&
        loop_boundaries,
    const Lambda& lambda) {
#pragma acc loop worker
  for (iType j = loop_boundaries.start; j < loop_boundaries.end; j++) {
    lambda(j);
  }
}

// Hierarchical Parallelism -> Thread vector level implementation
#pragma acc routine vector
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::OpenACCTeamMember>& loop_boundaries,
    const Lambda& lambda) {
#pragma acc loop vector
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i);
  }
}

// Hierarchical Parallelism -> Team vector level implementation
#pragma acc routine vector
template <typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamVectorRangeBoundariesStruct<iType, Impl::OpenACCTeamMember>&
        loop_boundaries,
    const Lambda& lambda) {
#pragma acc loop vector
  for (iType i = loop_boundaries.start; i < loop_boundaries.end; i++) {
    lambda(i);
  }
}

}  // namespace Kokkos

#endif /* #ifdef KOKKOS_ENABLE_OPENACC_COLLAPSE_HIERARCHICAL_CONSTRUCTS */

#endif /* #ifndef KOKKOS_OPENACC_PARALLEL_FOR_TEAM_HPP */
