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

#ifndef KOKKOS_HPX_TASK_HPP
#define KOKKOS_HPX_TASK_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_HPX) && defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_TaskScheduler_fwd.hpp>

#include <HPX/Kokkos_HPX.hpp>

#include <hpx/local/execution.hpp>
#include <hpx/local/future.hpp>

#include <type_traits>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class QueueType>
class TaskQueueSpecialization<
    SimpleTaskScheduler<Kokkos::Experimental::HPX, QueueType>> {
 public:
  void setup() const {
    const int num_worker_threads = Kokkos::Experimental::HPX().concurrency();

    hpx_thread_buffer &buffer = Kokkos::Experimental::HPX().impl_get_buffer();
    buffer.resize(num_worker_threads, 512);
  }

  void execute_range(int t) const {
    // NOTE: This implementation has been simplified based on the
    // assumption that team_size = 1. The HPX backend currently only
    // supports a team size of 1.
    const int num_worker_threads = Kokkos::Experimental::HPX().concurrency();

    hpx_thread_buffer &buffer = Kokkos::Experimental::HPX().impl_get_buffer();

    buffer.get(t);
    HPXTeamMember member(
        TeamPolicyInternal<Kokkos::Experimental::HPX>(
            Kokkos::Experimental::HPX(), num_worker_threads, 1),
        0, t, buffer.get(t), 512);

    member_type single_exec(*scheduler, member);
    member_type &team_exec = single_exec;

    auto &queue          = scheduler->queue();
    auto &team_scheduler = team_exec.scheduler();

    using task_base_type = typename scheduler_type::task_base_type;
    auto current_task    = OptionalRef<task_base_type>(nullptr);

    while (!queue.is_done()) {
      current_task = queue.pop_ready_task(team_scheduler.team_scheduler_info());

      if (current_task) {
        KOKKOS_EXPECTS(current_task->is_single_runnable() ||
                       current_task->is_team_runnable());
        current_task->as_runnable_task().run(single_exec);
        queue.complete((*std::move(current_task)).as_runnable_task(),
                       team_scheduler.team_scheduler_info());
      }
    }
  }

  void finalize() const {}

  using execution_space = Kokkos::Experimental::HPX;
  using scheduler_type =
      SimpleTaskScheduler<Kokkos::Experimental::HPX, QueueType>;
  using member_type =
      TaskTeamMemberAdapter<Kokkos::Impl::HPXTeamMember, scheduler_type>;
  using memory_space = Kokkos::HostSpace;

  static void execute(scheduler_type const &scheduler) {
    // NOTE: We create an instance so that we can use impl_bulk_setup_finalize.
    // This is not necessarily the most efficient, but can be improved later.
    TaskQueueSpecialization<scheduler_type> task_queue;
    task_queue.scheduler         = &scheduler;
    const int num_worker_threads = Kokkos::Experimental::HPX().concurrency();
    Kokkos::Experimental::HPX().impl_bulk_setup_finalize(
        true, false, task_queue, num_worker_threads,
        hpx::threads::thread_stacksize::nostack);
  }

  static uint32_t get_max_team_count(execution_space const &espace) {
    return static_cast<uint32_t>(espace.concurrency());
  }

  template <typename TaskType>
  static void get_function_pointer(typename TaskType::function_type &ptr,
                                   typename TaskType::destroy_type &dtor) {
    ptr  = TaskType::apply;
    dtor = TaskType::destroy;
  }

 private:
  const scheduler_type *scheduler;
};

template <class Scheduler>
class TaskQueueSpecializationConstrained<
    Scheduler,
    std::enable_if_t<std::is_same<typename Scheduler::execution_space,
                                  Kokkos::Experimental::HPX>::value>> {
 public:
  void setup() const {
    const int num_worker_threads = Kokkos::Experimental::HPX().concurrency();

    hpx_thread_buffer &buffer = Kokkos::Experimental::HPX().impl_get_buffer();
    buffer.resize(num_worker_threads, 512);

    auto &queue = scheduler->queue();
    queue.initialize_team_queues(num_worker_threads);
  }

  void execute_range(int t) const {
    // NOTE: This implementation has been simplified based on the
    // assumption that team_size = 1. The HPX backend currently only
    // supports a team size of 1.
    const int num_worker_threads = Kokkos::Experimental::HPX().concurrency();

    hpx_thread_buffer &buffer = Kokkos::Experimental::HPX().impl_get_buffer();

    buffer.get(Kokkos::Experimental::HPX::impl_hardware_thread_id());
    HPXTeamMember member(
        TeamPolicyInternal<Kokkos::Experimental::HPX>(
            Kokkos::Experimental::HPX(), num_worker_threads, 1),
        0, t, buffer.get(t), 512);

    using task_base_type = typename scheduler_type::task_base;
    using queue_type     = typename scheduler_type::queue_type;

    static task_base_type *const end = (task_base_type *)task_base_type::EndTag;
    constexpr task_base_type *no_more_tasks_sentinel = nullptr;

    member_type single_exec(*scheduler, member);
    member_type &team_exec = single_exec;

    auto &team_queue     = team_exec.scheduler().queue();
    task_base_type *task = no_more_tasks_sentinel;

    do {
      if (task != no_more_tasks_sentinel && task != end) {
        team_queue.complete(task);
      }

      if (*((volatile int *)&team_queue.m_ready_count) > 0) {
        task = end;
        for (int i = 0; i < queue_type::NumQueue && end == task; ++i) {
          for (int j = 0; j < 2 && end == task; ++j) {
            task = queue_type::pop_ready_task(&team_queue.m_ready[i][j]);
          }
        }
      } else {
        task = team_queue.attempt_to_steal_task();
      }

      if (task != no_more_tasks_sentinel && task != end) {
        (*task->m_apply)(task, &single_exec);
      }
    } while (task != no_more_tasks_sentinel);
  }

  void finalize() const {}

  using execution_space = Kokkos::Experimental::HPX;
  using scheduler_type  = Scheduler;
  using member_type =
      TaskTeamMemberAdapter<Kokkos::Impl::HPXTeamMember, scheduler_type>;
  using memory_space = Kokkos::HostSpace;

  static void iff_single_thread_recursive_execute(
      scheduler_type const &scheduler) {
    using task_base_type = typename scheduler_type::task_base;
    using queue_type     = typename scheduler_type::queue_type;

    if (1 == Kokkos::Experimental::HPX().concurrency()) {
      task_base_type *const end = (task_base_type *)task_base_type::EndTag;
      task_base_type *task      = end;

      HPXTeamMember member(TeamPolicyInternal<Kokkos::Experimental::HPX>(
                               Kokkos::Experimental::HPX(), 1, 1),
                           0, 0, nullptr, 0);
      member_type single_exec(scheduler, member);

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

      } while (true);
    }
  }

  static void execute(scheduler_type const &scheduler) {
    // NOTE: We create an instance so that we can use impl_bulk_setup_finalize.
    // This is not necessarily the most efficient, but can be improved later.
    TaskQueueSpecializationConstrained<scheduler_type> task_queue;
    task_queue.scheduler         = &scheduler;
    const int num_worker_threads = Kokkos::Experimental::HPX().concurrency();
    Kokkos::Experimental::HPX().impl_bulk_setup_finalize(
        true, false, task_queue, num_worker_threads,
        hpx::threads::thread_stacksize::nostack);
  }

  template <typename TaskType>
  static void get_function_pointer(typename TaskType::function_type &ptr,
                                   typename TaskType::destroy_type &dtor) {
    ptr  = TaskType::apply;
    dtor = TaskType::destroy;
  }

 private:
  const scheduler_type *scheduler;
};

extern template class TaskQueue<
    Kokkos::Experimental::HPX,
    typename Kokkos::Experimental::HPX::memory_space>;

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_HPX_TASK_HPP */
