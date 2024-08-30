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

#ifndef KOKKOS_TASKTEAMMEMBER_HPP
#define KOKKOS_TASKTEAMMEMBER_HPP

//----------------------------------------------------------------------------

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_TaskScheduler_fwd.hpp>
//----------------------------------------------------------------------------

#include <Kokkos_MemoryPool.hpp>

#include <Kokkos_Future.hpp>
#include <impl/Kokkos_TaskQueue.hpp>
#include <impl/Kokkos_SingleTaskQueue.hpp>
#include <impl/Kokkos_TaskQueueMultiple.hpp>
#include <impl/Kokkos_TaskPolicyData.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class TeamMember, class Scheduler>
class TaskTeamMemberAdapter : public TeamMember {
 private:
  Scheduler m_scheduler;

 public:
  //----------------------------------------

  // Forward everything but the Scheduler to the constructor of the TeamMember
  // type that we're adapting
  template <typename... Args>
  KOKKOS_INLINE_FUNCTION explicit TaskTeamMemberAdapter(
      std::enable_if_t<std::is_constructible<TeamMember, Args...>::value,
                       Scheduler>
          arg_scheduler,
      Args&&... args)  // TODO @tasking @minor DSH noexcept specification
      : TeamMember(std::forward<Args>(args)...),
        m_scheduler(
            std::move(arg_scheduler).get_team_scheduler(this->league_rank())) {}

  // (rule of 6 constructors)

  KOKKOS_DEFAULTED_FUNCTION
  TaskTeamMemberAdapter() = default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskTeamMemberAdapter(TaskTeamMemberAdapter const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskTeamMemberAdapter(TaskTeamMemberAdapter&&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskTeamMemberAdapter& operator=(TaskTeamMemberAdapter const&) = default;

  KOKKOS_DEFAULTED_FUNCTION
  TaskTeamMemberAdapter& operator=(TaskTeamMemberAdapter&&) = default;

  KOKKOS_DEFAULTED_FUNCTION ~TaskTeamMemberAdapter() = default;

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  Scheduler const& scheduler() const noexcept { return m_scheduler; }

  KOKKOS_INLINE_FUNCTION
  Scheduler& scheduler() noexcept { return m_scheduler; }

  //----------------------------------------
};

}  // end namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_TASKTEAMMEMBER_HPP */
