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

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_OPENMPTARGET) && defined(KOKKOS_ENABLE_TASKPOLICY)

#include <impl/Kokkos_TaskQueue_impl.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template class TaskQueue<Kokkos::Experimental::OpenMPTarget>;

//----------------------------------------------------------------------------

TaskExec<Kokkos::Experimental::OpenMPTarget>::TaskExec()
    : m_self_exec(0),
      m_team_exec(0),
      m_sync_mask(0),
      m_sync_value(0),
      m_sync_step(0),
      m_group_rank(0),
      m_team_rank(0),
      m_team_size(1) {}

TaskExec<Kokkos::Experimental::OpenMPTarget>::TaskExec(
    Kokkos::Impl::OpenMPTargetExec &arg_exec, int const arg_team_size)
    : m_self_exec(&arg_exec),
      m_team_exec(arg_exec.pool_rev(arg_exec.pool_rank_rev() / arg_team_size)),
      m_sync_mask(0),
      m_sync_value(0),
      m_sync_step(0),
      m_group_rank(arg_exec.pool_rank_rev() / arg_team_size),
      m_team_rank(arg_exec.pool_rank_rev() % arg_team_size),
      m_team_size(arg_team_size) {
  // This team spans
  //    m_self_exec->pool_rev( team_size * group_rank )
  //    m_self_exec->pool_rev( team_size * ( group_rank + 1 ) - 1 )

  int64_t volatile *const sync = (int64_t *)m_self_exec->scratch_reduce();

  sync[0] = int64_t(0);
  sync[1] = int64_t(0);

  for (int i = 0; i < m_team_size; ++i) {
    m_sync_value |= int64_t(1) << (8 * i);
    m_sync_mask |= int64_t(3) << (8 * i);
  }

  Kokkos::memory_fence();
}

void TaskExec<Kokkos::Experimental::OpenMPTarget>::team_barrier_impl() const {
  if (m_team_exec->scratch_reduce_size() < int(2 * sizeof(int64_t))) {
    Kokkos::abort("TaskQueue<OpenMPTarget> scratch_reduce memory too small");
  }

  // Use team shared memory to synchronize.
  // Alternate memory locations between barriers to avoid a sequence
  // of barriers overtaking one another.

  int64_t volatile *const sync =
      ((int64_t *)m_team_exec->scratch_reduce()) + (m_sync_step & 0x01);

  // This team member sets one byte within the sync variable
  int8_t volatile *const sync_self = ((int8_t *)sync) + m_team_rank;

  *sync_self = int8_t(m_sync_value & 0x03);  // signal arrival

  while (m_sync_value != *sync)
    ;  // wait for team to arrive

  ++m_sync_step;

  if (0 == (0x01 & m_sync_step)) {  // Every other step
    m_sync_value ^= m_sync_mask;
    if (1000 < m_sync_step) m_sync_step = 0;
  }
}

//----------------------------------------------------------------------------

void TaskQueueSpecialization<Kokkos::Experimental::OpenMPTarget>::execute(
    TaskQueue<Kokkos::Experimental::OpenMPTarget> *const queue) {
  using execution_space = Kokkos::Experimental::OpenMPTarget;
  using queue_type      = TaskQueue<execution_space>;
  using task_root_type  = TaskBase<execution_space, void, void>;
  using PoolExec        = Kokkos::Impl::OpenMPTargetExec;
  using Member          = TaskExec<execution_space>;

  task_root_type *const end = (task_root_type *)task_root_type::EndTag;

  // Required:  team_size <= 8

  const int team_size = PoolExec::pool_size(2);  // Threads per core
  // const int team_size = PoolExec::pool_size(1); // Threads per NUMA

  if (8 < team_size) {
    Kokkos::abort("TaskQueue<OpenMPTarget> unsupported team size");
  }

#pragma omp parallel
  {
    PoolExec &self = *PoolExec::get_thread_omp();

    Member single_exec;
    Member team_exec(self, team_size);

    // Team shared memory
    task_root_type *volatile *const task_shared =
        (task_root_type **)team_exec.m_team_exec->scratch_thread();

// Barrier across entire OpenMPTarget thread pool to insure initialization
#pragma omp barrier

    // Loop until all queues are empty and no tasks in flight

    do {
      task_root_type *task = 0;

      // Each team lead attempts to acquire either a thread team task
      // or a single thread task for the team.

      if (0 == team_exec.team_rank()) {
        task = 0 < *((volatile int *)&queue->m_ready_count) ? end : 0;

        // Loop by priority and then type
        for (int i = 0; i < queue_type::NumQueue && end == task; ++i) {
          for (int j = 0; j < 2 && end == task; ++j) {
            task = queue_type::pop_task(&queue->m_ready[i][j]);
          }
        }
      }

      // Team lead broadcast acquired task to team members:

      if (1 < team_exec.team_size()) {
        if (0 == team_exec.team_rank()) *task_shared = task;

        // Fence to be sure task_shared is stored before the barrier
        Kokkos::memory_fence();

        // Whole team waits for every team member to reach this statement
        team_exec.team_barrier();

        // Fence to be sure task_shared is stored
        Kokkos::memory_fence();

        task = *task_shared;
      }

      if (0 == task) break;  // 0 == m_ready_count

      if (end == task) {
        // All team members wait for whole team to reach this statement.
        // Is necessary to prevent task_shared from being updated
        // before it is read by all threads.
        team_exec.team_barrier();
      } else if (task_root_type::TaskTeam == task->m_task_type) {
        // Thread Team Task
        (*task->m_apply)(task, &team_exec);

        // The m_apply function performs a barrier

        if (0 == team_exec.team_rank()) {
          // team member #0 completes the task, which may delete the task
          queue->complete(task);
        }
      } else {
        // Single Thread Task

        if (0 == team_exec.team_rank()) {
          (*task->m_apply)(task, &single_exec);

          queue->complete(task);
        }

        // All team members wait for whole team to reach this statement.
        // Not necessary to complete the task.
        // Is necessary to prevent task_shared from being updated
        // before it is read by all threads.
        team_exec.team_barrier();
      }
    } while (1);
  }
  // END #pragma omp parallel
}

void TaskQueueSpecialization<Kokkos::Experimental::OpenMPTarget>::
    iff_single_thread_recursive_execute(
        TaskQueue<Kokkos::Experimental::OpenMPTarget> *const queue) {
  using execution_space = Kokkos::Experimental::OpenMPTarget;
  using queue_type      = TaskQueue<execution_space>;
  using task_root_type  = TaskBase<execution_space, void, void>;
  using Member          = TaskExec<execution_space>;

  if (1 == omp_get_num_threads()) {
    task_root_type *const end = (task_root_type *)task_root_type::EndTag;

    Member single_exec;

    task_root_type *task = end;

    do {
      task = end;

      // Loop by priority and then type
      for (int i = 0; i < queue_type::NumQueue && end == task; ++i) {
        for (int j = 0; j < 2 && end == task; ++j) {
          task = queue_type::pop_task(&queue->m_ready[i][j]);
        }
      }

      if (end == task) break;

      (*task->m_apply)(task, &single_exec);

      queue->complete(task);

    } while (1);
  }
}

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_OPENMPTARGET ) && defined( \
          KOKKOS_ENABLE_TASKPOLICY ) */
