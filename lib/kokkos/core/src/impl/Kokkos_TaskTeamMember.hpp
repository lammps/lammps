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
