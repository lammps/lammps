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

#ifndef KOKKOS_IMPL_OPENMP_TASK_HPP
#define KOKKOS_IMPL_OPENMP_TASK_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_OPENMP) && defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_Atomic.hpp>
#include <Kokkos_TaskScheduler_fwd.hpp>

#include <impl/Kokkos_HostThreadTeam.hpp>
#include <OpenMP/Kokkos_OpenMP.hpp>

#include <type_traits>
#include <cassert>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

class HostThreadTeamDataSingleton : private HostThreadTeamData {
 private:
  HostThreadTeamDataSingleton();
  ~HostThreadTeamDataSingleton();

 public:
  static HostThreadTeamData& singleton();
};

// Hack this as a partial specialization for now
// TODO @tasking @cleanup DSH Make this the general class template and make the
// old code the partial specialization
template <class QueueType>
class TaskQueueSpecialization<SimpleTaskScheduler<Kokkos::OpenMP, QueueType>> {
 public:
  using execution_space = Kokkos::OpenMP;
  using scheduler_type  = SimpleTaskScheduler<Kokkos::OpenMP, QueueType>;
  using member_type =
      TaskTeamMemberAdapter<Kokkos::Impl::HostThreadTeamMember<execution_space>,
                            scheduler_type>;
  using memory_space = Kokkos::HostSpace;

  enum : int { max_league_size = HostThreadTeamData::max_pool_members };

  // Must provide task queue execution function
  static void execute(scheduler_type const& scheduler) {
    using task_base_type = typename scheduler_type::task_base_type;

    // Unused; ChaseLev queue still needs worker ID even in single case (so we
    // need to use the thread data from inside of the parallel region.  Team
    // size is fixed at 1 for now anyway
    // HostThreadTeamData& team_data_single =
    // HostThreadTeamDataSingleton::singleton();

    Impl::OpenMPInternal* instance =
        execution_space().impl_internal_space_instance();
    const int pool_size = get_max_team_count(scheduler.get_execution_space());

    instance->acquire_lock();

    // TODO @tasking @new_feature DSH allow team sizes other than 1
    const int team_size = 1;                      // Threads per core
    instance->resize_thread_data(0,               /* global reduce buffer */
                                 512 * team_size, /* team reduce buffer */
                                 0,               /* team shared buffer */
                                 0                /* thread local buffer */
    );
    assert(pool_size % team_size == 0);

    auto& queue = scheduler.queue();

    // queue.initialize_team_queues(pool_size / team_size);

#pragma omp parallel num_threads(pool_size)
    {
      Impl::HostThreadTeamData& self = *(instance->get_thread_data());

      // Organizing threads into a team performs a barrier across the
      // entire pool to insure proper initialization of the team
      // rendezvous mechanism before a team rendezvous can be performed.

      // organize_team() returns true if this is an active team member
      if (self.organize_team(team_size)) {
        member_type single_exec(scheduler, self);
        member_type team_exec(scheduler, self);

        auto& team_scheduler = team_exec.scheduler();

        auto current_task = OptionalRef<task_base_type>(nullptr);

        while (!queue.is_done()) {
          // Each team lead attempts to acquire either a thread team task
          // or a single thread task for the team.
          if (team_exec.team_rank() == 0) {
            // loop while both:
            //   - the queue is not done
            //   - the most recently popped task is a single task or empty
            while (!queue.is_done()) {
              current_task =
                  queue.pop_ready_task(team_scheduler.team_scheduler_info());

              if (current_task) {
                if (current_task->is_team_runnable()) {
                  // break out of the team leader loop to run the team task
                  break;
                } else {
                  KOKKOS_ASSERT(current_task->is_single_runnable());
                  current_task->as_runnable_task().run(single_exec);
                  // Respawns are handled in the complete function
                  queue.complete((*std::move(current_task)).as_runnable_task(),
                                 team_scheduler.team_scheduler_info());
                }

              }  // end if current_task is not null

              current_task = nullptr;

            }  // end team leader loop
          }

          // Otherwise, make sure everyone in the team has the same task
          team_exec.team_broadcast(current_task, 0);

          if (current_task) {
            KOKKOS_ASSERT(current_task->is_team_runnable());
            current_task->as_runnable_task().run(team_exec);

            if (team_exec.team_rank() == 0) {
              // Respawns are handled in the complete function
              queue.complete((*std::move(current_task)).as_runnable_task(),
                             team_scheduler.team_scheduler_info());
            }
          }
        }
      }
      self.disband_team();
    }  // end pragma omp parallel

    instance->release_lock();
  }

  static uint32_t get_max_team_count(execution_space const& espace) {
    return static_cast<uint32_t>(espace.impl_thread_pool_size());
  }

  // TODO @tasking @optimization DSH specialize this for trivially destructible
  // types
  template <typename TaskType>
  static void get_function_pointer(typename TaskType::function_type& ptr,
                                   typename TaskType::destroy_type& dtor) {
    ptr  = TaskType::apply;
    dtor = TaskType::destroy;
  }
};

template <class Scheduler>
class TaskQueueSpecializationConstrained<
    Scheduler,
    std::enable_if_t<std::is_same<typename Scheduler::execution_space,
                                  Kokkos::OpenMP>::value>> {
 public:
  using execution_space = Kokkos::OpenMP;
  using scheduler_type  = Scheduler;
  using member_type =
      TaskTeamMemberAdapter<Kokkos::Impl::HostThreadTeamMember<execution_space>,
                            scheduler_type>;
  using memory_space = Kokkos::HostSpace;

  enum : int { max_league_size = HostThreadTeamData::max_pool_members };

  static void iff_single_thread_recursive_execute(
      scheduler_type const& scheduler) {
    using task_base_type = typename scheduler_type::task_base;
    using queue_type     = typename scheduler_type::queue_type;

    execution_space exec;
    if (1 == exec.impl_thread_pool_size()) {
      task_base_type* const end = (task_base_type*)task_base_type::EndTag;

      HostThreadTeamData& team_data_single =
          HostThreadTeamDataSingleton::singleton();

      member_type single_exec(scheduler, team_data_single);

      task_base_type* task = end;

      do {
        task = end;

        // Loop by priority and then type
        for (int i = 0; i < queue_type::NumQueue && end == task; ++i) {
          for (int j = 0; j < 2 && end == task; ++j) {
            task =
                queue_type::pop_ready_task(&scheduler.m_queue->m_ready[i][j]);
          }
        }

        if (end == task) break;

        (*task->m_apply)(task, &single_exec);

        scheduler.m_queue->complete(task);

      } while (1);
    }
  }

  // Must provide task queue execution function
  static void execute(scheduler_type const& scheduler) {
    using task_base_type = typename scheduler_type::task_base;
    using queue_type     = typename scheduler_type::queue_type;

    static task_base_type* const end = (task_base_type*)task_base_type::EndTag;

    constexpr task_base_type* no_more_tasks_sentinel = nullptr;

    HostThreadTeamData& team_data_single =
        HostThreadTeamDataSingleton::singleton();

    Impl::OpenMPInternal* instance =
        execution_space().impl_internal_space_instance();
    const int pool_size = instance->thread_pool_size();

    instance->acquire_lock();

    const int team_size = 1;       // Threads per core
    instance->resize_thread_data(0 /* global reduce buffer */
                                 ,
                                 512 * team_size /* team reduce buffer */
                                 ,
                                 0 /* team shared buffer */
                                 ,
                                 0 /* thread local buffer */
    );
    assert(pool_size % team_size == 0);
    auto& queue = scheduler.queue();
    queue.initialize_team_queues(pool_size / team_size);

#pragma omp parallel num_threads(pool_size)
    {
      Impl::HostThreadTeamData& self = *(instance->get_thread_data());

      // Organizing threads into a team performs a barrier across the
      // entire pool to insure proper initialization of the team
      // rendezvous mechanism before a team rendezvous can be performed.

      // organize_team() returns true if this is an active team member
      if (self.organize_team(team_size)) {
        member_type single_exec(scheduler, team_data_single);
        member_type team_exec(scheduler, self);

        auto& team_queue = team_exec.scheduler().queue();

        // Loop until all queues are empty and no tasks in flight

        task_base_type* task = no_more_tasks_sentinel;

        do {
          // Each team lead attempts to acquire either a thread team task
          // or a single thread task for the team.

          if (0 == team_exec.team_rank()) {
            bool leader_loop = false;

            do {
              if (task != no_more_tasks_sentinel && task != end) {
                // team member #0 completes the previously executed task,
                // completion may delete the task
                team_queue.complete(task);
              }

              // If 0 == m_ready_count then set task = 0

              if (desul::atomic_load(&team_queue.m_ready_count,
                                     desul::MemoryOrderAcquire(),
                                     desul::MemoryScopeDevice()) > 0) {
                task = end;
                // Attempt to acquire a task
                // Loop by priority and then type
                for (int i = 0; i < queue_type::NumQueue && end == task; ++i) {
                  for (int j = 0; j < 2 && end == task; ++j) {
                    task =
                        queue_type::pop_ready_task(&team_queue.m_ready[i][j]);
                  }
                }
              } else {
                // returns nullptr if and only if all other queues have a ready
                // count of 0 also. Otherwise, returns a task from another queue
                // or `end` if one couldn't be popped
                task = team_queue.attempt_to_steal_task();
              }

              // If still tasks are still executing
              // and no task could be acquired
              // then continue this leader loop
              if (task == end) {
                // this means that the ready task count was not zero, but we
                // couldn't pop a task (because, for instance, someone else
                // got there before us
                leader_loop = true;
              } else if ((task != no_more_tasks_sentinel) &&
                         (task_base_type::TaskSingle == task->m_task_type)) {
                // if a single thread task then execute now

                (*task->m_apply)(task, &single_exec);

                leader_loop = true;
              } else {
                leader_loop = false;
              }
            } while (leader_loop);
          }

          // Team lead either found 0 == m_ready_count or a team task
          // Team lead broadcast acquired task:

          team_exec.team_broadcast(task, 0);

          if (task != no_more_tasks_sentinel) {  // Thread Team Task

            (*task->m_apply)(task, &team_exec);

            // The m_apply function performs a barrier
          }
        } while (task != no_more_tasks_sentinel);
      }
      self.disband_team();
    }  // end pragma omp parallel

    instance->release_lock();
  }

  template <typename TaskType>
  static void get_function_pointer(typename TaskType::function_type& ptr,
                                   typename TaskType::destroy_type& dtor) {
    ptr  = TaskType::apply;
    dtor = TaskType::destroy;
  }
};

extern template class TaskQueue<Kokkos::OpenMP,
                                typename Kokkos::OpenMP::memory_space>;

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_OPENMP_TASK_HPP */
