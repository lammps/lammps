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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_IMPL_TEAMMDPOLICY_HPP
#define KOKKOS_IMPL_TEAMMDPOLICY_HPP

namespace Kokkos {

namespace Impl {

// Tag class to choose the nested loop specialization
//   - LastNestLevel means call the actual closure
//   - ParThread means use TeamThreadRange
//   - ParVector means use ThreadVectorRange
template <TeamMDRangeLastNestLevel LastNestLevel,
          TeamMDRangeParThread ParThread, TeamMDRangeParVector ParVector>
struct TeamMDRangeMode {
  static constexpr TeamMDRangeLastNestLevel last_nest_level = LastNestLevel;
  static constexpr TeamMDRangeParThread par_thread          = ParThread;
  static constexpr TeamMDRangeParVector par_vector          = ParVector;
};

// Tag class to keep track of the loop nest level and where to deploy thread and
// vector parallelism
//   - Rank is Kokkos::Rank<TotalNestLevel, Iter>
//     - total_nest_level is the total number of loop nests
//     - iter is whether to go forward or backward through ranks (i.e. the
//       iteration order for MDRangePolicy)
//   - ParThreadNestLevel is the nesting level on which to deploy thread
//   parallelism
//   - ParVectorNestLevel is the nesting level on which to deploy vector
//   parallelism
//   - CurrentNestLevel is the nest level of interest
template <typename Rank, int ParThreadNestLevel, int ParVectorNestLevel,
          int CurrentNestLevel>
struct TeamMDRangeNestingTracker {
  using NestLevelType                                  = int;
  static constexpr Iterate iter                        = Rank::outer_direction;
  static constexpr NestLevelType total_nest_level      = Rank::rank;
  static constexpr NestLevelType par_thread_nest_level = ParThreadNestLevel;
  static constexpr NestLevelType par_vector_nest_level = ParVectorNestLevel;
  static constexpr NestLevelType current_nest_level    = CurrentNestLevel;

  // We have to recursively process ranks [0..total_nest_level-1]
  using RangeMode =
      TeamMDRangeMode<(iter == Iterate::Right)
                          ? static_cast<TeamMDRangeLastNestLevel>(
                                current_nest_level == total_nest_level)
                          : static_cast<TeamMDRangeLastNestLevel>(
                                current_nest_level == -1),
                      static_cast<TeamMDRangeParThread>(current_nest_level ==
                                                        par_thread_nest_level),
                      static_cast<TeamMDRangeParVector>(current_nest_level ==
                                                        par_vector_nest_level)>;
};

// Structs to determine on which nested level parallelization happens.
//   - Rank is Kokkos::Rank<TotalNestLevel, Iter>
//     - TotalNestLevel is the total number of loop nests
//     - Iter is whether to go forward or backward through ranks (i.e. the
//       iteration order for MDRangePolicy)
//   - ThreadAndVector determines whether both vector and thread parallelism is
//     in use
template <typename Rank, TeamMDRangeThreadAndVector ThreadAndVector>
struct HostBasedNestLevel {
  static constexpr bool is_direction_left =
      (Rank::outer_direction == Iterate::Left);
  static constexpr int par_rt  = is_direction_left ? Rank::rank - 1 : 0;
  static constexpr int par_rv  = is_direction_left ? 0 : Rank::rank - 1;
  static constexpr int invalid = -2;
};

template <typename Rank, TeamMDRangeThreadAndVector ThreadAndVector>
struct AcceleratorBasedNestLevel {
  static constexpr bool is_direction_left =
      (Rank::outer_direction == Iterate::Left);

  // If vector parallelism is in use, deploy thread parallelism on
  // the second to the last nested level; otherwise, thread parallelism on the
  // last nested level
  static constexpr int left_par_rt =
      (ThreadAndVector == TeamMDRangeThreadAndVector::Both) ? 1 : 0;

  static constexpr int right_par_rt =
      (ThreadAndVector == TeamMDRangeThreadAndVector::Both) ? Rank::rank - 2
                                                            : Rank::rank - 1;

  static constexpr int par_rt = is_direction_left ? left_par_rt : right_par_rt;

  // Vector parallelism will always be on the last index
  static constexpr int par_rv  = is_direction_left ? 0 : Rank::rank - 1;
  static constexpr int invalid = -2;
};

template <typename TeamHandle>
KOKKOS_INLINE_FUNCTION auto nested_policy(
    TeamMDRangeMode<TeamMDRangeLastNestLevel::NotLastNestLevel,
                    TeamMDRangeParThread::ParThread,
                    TeamMDRangeParVector::NotParVector>,
    TeamHandle const team, int count) {
  return TeamThreadRange(team, count);
}

template <typename TeamHandle>
KOKKOS_INLINE_FUNCTION auto nested_policy(
    TeamMDRangeMode<TeamMDRangeLastNestLevel::NotLastNestLevel,
                    TeamMDRangeParThread::NotParThread,
                    TeamMDRangeParVector::ParVector>,
    TeamHandle const team, int count) {
  return ThreadVectorRange(team, count);
}

template <typename TeamHandle>
KOKKOS_INLINE_FUNCTION auto nested_policy(
    TeamMDRangeMode<TeamMDRangeLastNestLevel::NotLastNestLevel,
                    TeamMDRangeParThread::ParThread,
                    TeamMDRangeParVector::ParVector>,
    TeamHandle const team, int count) {
  return TeamVectorRange(team, count);
}

// TeamMDRangeNestingTracker is only needed to deduce template parameters
template <typename Rank, int ParThreadNestLevel, int ParVectorNestLevel,
          int CurrentNestLevel, typename Policy, typename Lambda,
          typename... Args>
KOKKOS_INLINE_FUNCTION void nested_loop(
    TeamMDRangeMode<TeamMDRangeLastNestLevel::LastNestLevel,
                    TeamMDRangeParThread::NotParThread,
                    TeamMDRangeParVector::NotParVector> const,
    TeamMDRangeNestingTracker<Rank, ParThreadNestLevel, ParVectorNestLevel,
                              CurrentNestLevel>,
    Policy const&, Lambda const& lambda, Impl::NoReductionTag&&, Args... args) {
  lambda(args...);
}

template <typename Rank, int ParThreadNestLevel, int ParVectorNestLevel,
          int CurrentNestLevel, typename Policy, typename Lambda,
          typename ReducerValueType, typename... Args>
KOKKOS_INLINE_FUNCTION void nested_loop(
    TeamMDRangeMode<TeamMDRangeLastNestLevel::LastNestLevel,
                    TeamMDRangeParThread::NotParThread,
                    TeamMDRangeParVector::NotParVector> const,
    TeamMDRangeNestingTracker<Rank, ParThreadNestLevel, ParVectorNestLevel,
                              CurrentNestLevel>,
    Policy const&, Lambda const& lambda, ReducerValueType& val, Args... args) {
  lambda(args..., val);
}

// Nested loop for serial iteration
template <typename Rank, int ParThreadNestLevel, int ParVectorNestLevel,
          int CurrentNestLevel, typename Policy, typename Lambda,
          typename ReducerValueType, typename... Args>
KOKKOS_INLINE_FUNCTION void nested_loop(
    TeamMDRangeMode<TeamMDRangeLastNestLevel::NotLastNestLevel,
                    TeamMDRangeParThread::NotParThread,
                    TeamMDRangeParVector::NotParVector> const,
    TeamMDRangeNestingTracker<Rank, ParThreadNestLevel, ParVectorNestLevel,
                              CurrentNestLevel>,
    Policy const& policy, Lambda const& lambda, ReducerValueType&& val,
    Args... args) {
  constexpr int next_nest_level =
      CurrentNestLevel + (Rank::outer_direction == Iterate::Right ? 1 : -1);
  using NextNestingTracker =
      TeamMDRangeNestingTracker<Rank, ParThreadNestLevel, ParVectorNestLevel,
                                next_nest_level>;
  using TeamMDNextMode = typename NextNestingTracker::RangeMode;

  for (int i = 0; i != policy.boundaries[CurrentNestLevel]; ++i) {
    if constexpr (Rank::outer_direction == Iterate::Right) {
      nested_loop(TeamMDNextMode(), NextNestingTracker(), policy, lambda,
                  std::forward<ReducerValueType>(val), args..., i);
    } else {
      nested_loop(TeamMDNextMode(), NextNestingTracker(), policy, lambda,
                  std::forward<ReducerValueType>(val), i, args...);
    }
  }
}

template <TeamMDRangeParThread ParThread, TeamMDRangeParVector ParVector,
          typename Rank, int ParThreadNestLevel, int ParVectorNestLevel,
          int CurrentNestLevel, typename Policy, typename Lambda,
          typename ReducerValueType, typename... Args>
KOKKOS_INLINE_FUNCTION void nested_loop(
    TeamMDRangeMode<TeamMDRangeLastNestLevel::NotLastNestLevel, ParThread,
                    ParVector> const mode,
    TeamMDRangeNestingTracker<Rank, ParThreadNestLevel, ParVectorNestLevel,
                              CurrentNestLevel>,
    Policy const& policy, Lambda const& lambda, ReducerValueType&& val,
    Args... args) {
  constexpr int next_nest_level =
      CurrentNestLevel + (Rank::outer_direction == Iterate::Right ? 1 : -1);
  using NextNestingTracker =
      TeamMDRangeNestingTracker<Rank, ParThreadNestLevel, ParVectorNestLevel,
                                next_nest_level>;
  using TeamMDNextMode = typename NextNestingTracker::RangeMode;

  // This recursively processes ranks from [0..TotalNestLevel-1]
  // args... is passed by value because it should always be ints
  parallel_for(
      nested_policy(mode, policy.team, policy.boundaries[CurrentNestLevel]),
      [&](int const& i) {
        if constexpr (Rank::outer_direction == Iterate::Right) {
          nested_loop(TeamMDNextMode(), NextNestingTracker(), policy, lambda,
                      std::forward<ReducerValueType>(val), args..., i);
        } else {
          nested_loop(TeamMDNextMode(), NextNestingTracker(), policy, lambda,
                      std::forward<ReducerValueType>(val), i, args...);
        }
      });
}

template <typename Rank, typename TeamMDPolicy, typename Lambda,
          typename ReductionValueType>
KOKKOS_INLINE_FUNCTION void md_parallel_impl(TeamMDPolicy const& policy,
                                             Lambda const& lambda,
                                             ReductionValueType&& val) {
  static_assert(TeamMDPolicy::total_nest_level >= 2 &&
                TeamMDPolicy::total_nest_level <= 8);

  using TeamHandle = typename TeamMDPolicy::TeamHandleType;

  constexpr auto total_nest_level = TeamMDPolicy::total_nest_level;
  constexpr auto iter             = TeamMDPolicy::iter;
  constexpr auto thread_and_vector =
      ((TeamMDPolicy::par_thread == Impl::TeamMDRangeParThread::ParThread) &&
       (TeamMDPolicy::par_vector == Impl::TeamMDRangeParVector::ParVector))
          ? Impl::TeamMDRangeThreadAndVector::Both
          : Impl::TeamMDRangeThreadAndVector::NotBoth;
  constexpr auto begin_rank =
      (iter == Iterate::Right) ? 0 : (total_nest_level - 1);

  using ThreadAndVectorNestLevel =
      Impl::ThreadAndVectorNestLevel<Rank, typename TeamHandle::execution_space,
                                     thread_and_vector>;

  constexpr auto par_thread_nest_level =
      (TeamMDPolicy::par_thread == TeamMDRangeParThread::ParThread)
          ? ThreadAndVectorNestLevel::par_rt
          : ThreadAndVectorNestLevel::invalid;
  constexpr auto par_vector_nest_level =
      (TeamMDPolicy::par_vector == TeamMDRangeParVector::ParVector)
          ? ThreadAndVectorNestLevel::par_rv
          : ThreadAndVectorNestLevel::invalid;

  using InitNestingTracker =
      TeamMDRangeNestingTracker<Rank, par_thread_nest_level,
                                par_vector_nest_level, begin_rank>;

  using InitTeamMDMode = typename InitNestingTracker::RangeMode;

  nested_loop(InitTeamMDMode(), InitNestingTracker(), policy, lambda,
              std::forward<ReductionValueType>(val));
}

}  // namespace Impl

}  // namespace Kokkos

#endif
