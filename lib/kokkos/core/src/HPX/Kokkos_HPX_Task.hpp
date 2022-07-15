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

#ifndef KOKKOS_HPX_TASK_HPP
#define KOKKOS_HPX_TASK_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_HPX) && defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_TaskScheduler_fwd.hpp>

#include <Kokkos_HPX.hpp>

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
  using execution_space = Kokkos::Experimental::HPX;
  using scheduler_type =
      SimpleTaskScheduler<Kokkos::Experimental::HPX, QueueType>;
  using member_type =
      TaskTeamMemberAdapter<Kokkos::Impl::HPXTeamMember, scheduler_type>;
  using memory_space = Kokkos::HostSpace;

  static void execute(scheduler_type const &scheduler) {
    // NOTE: We create an instance so that we can use dispatch_execute_task.
    // This is not necessarily the most efficient, but can be improved later.
    TaskQueueSpecialization<scheduler_type> task_queue;
    task_queue.scheduler = &scheduler;
    Kokkos::Impl::dispatch_execute_task(&task_queue,
                                        Kokkos::Experimental::HPX());
    Kokkos::Experimental::HPX().fence(
        "Kokkos::Impl::TaskQueueSpecialization<SimpleTask>::execute: fence "
        "after task execution");
  }

  // Must provide task queue execution function
  void execute_task() const {
    // See [note 1] in Kokkos_HPX.hpp for an explanation. The work graph policy
    // does not store an execution space instance, so we only need to reset the
    // parallel region count here.
    Kokkos::Experimental::HPX::reset_count_on_exit_parallel reset_count_on_exit;

    using hpx::for_loop;
    using hpx::execution::par;
    using hpx::execution::static_chunk_size;
    using task_base_type = typename scheduler_type::task_base_type;

    const int num_worker_threads = Kokkos::Experimental::HPX::concurrency();

    thread_buffer &buffer = Kokkos::Experimental::HPX().impl_get_buffer();
    buffer.resize(num_worker_threads, 512);

    auto &queue = scheduler->queue();

    for_loop(par.with(static_chunk_size(1)), 0, num_worker_threads,
             [this, &queue, &buffer, num_worker_threads](int) {
               // NOTE: This implementation has been simplified based on the
               // assumption that team_size = 1. The HPX backend currently only
               // supports a team size of 1.
               std::size_t t =
                   Kokkos::Experimental::HPX::impl_hardware_thread_id();

               buffer.get(t);
               HPXTeamMember member(
                   TeamPolicyInternal<Kokkos::Experimental::HPX>(
                       Kokkos::Experimental::HPX(), num_worker_threads, 1),
                   0, t, buffer.get(t), 512);

               member_type single_exec(*scheduler, member);
               member_type &team_exec = single_exec;

               auto &team_scheduler = team_exec.scheduler();
               auto current_task    = OptionalRef<task_base_type>(nullptr);

               while (!queue.is_done()) {
                 current_task =
                     queue.pop_ready_task(team_scheduler.team_scheduler_info());

                 if (current_task) {
                   KOKKOS_ASSERT(current_task->is_single_runnable() ||
                                 current_task->is_team_runnable());
                   current_task->as_runnable_task().run(single_exec);
                   queue.complete((*std::move(current_task)).as_runnable_task(),
                                  team_scheduler.team_scheduler_info());
                 }
               }
             });
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
    Scheduler, typename std::enable_if<
                   std::is_same<typename Scheduler::execution_space,
                                Kokkos::Experimental::HPX>::value>::type> {
 public:
  using execution_space = Kokkos::Experimental::HPX;
  using scheduler_type  = Scheduler;
  using member_type =
      TaskTeamMemberAdapter<Kokkos::Impl::HPXTeamMember, scheduler_type>;
  using memory_space = Kokkos::HostSpace;

  static void iff_single_thread_recursive_execute(
      scheduler_type const &scheduler) {
    using task_base_type = typename scheduler_type::task_base;
    using queue_type     = typename scheduler_type::queue_type;

    if (1 == Kokkos::Experimental::HPX::concurrency()) {
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
    // NOTE: We create an instance so that we can use dispatch_execute_task.
    // This is not necessarily the most efficient, but can be improved later.
    TaskQueueSpecializationConstrained<scheduler_type> task_queue;
    task_queue.scheduler = &scheduler;
    Kokkos::Impl::dispatch_execute_task(&task_queue,
                                        Kokkos::Experimental::HPX());
    Kokkos::Experimental::HPX().fence(
        "Kokkos::Impl::TaskQueueSpecialization<SimpleTask>::execute: fence "
        "after task execution");
  }

  // Must provide task queue execution function
  void execute_task() const {
    // See [note 1] in Kokkos_HPX.hpp for an explanation. The work graph policy
    // does not store an execution space instance, so we only need to reset the
    // parallel region count here.
    Kokkos::Experimental::HPX::reset_count_on_exit_parallel reset_count_on_exit;

    using hpx::for_loop;
    using hpx::execution::par;
    using hpx::execution::static_chunk_size;

    using task_base_type = typename scheduler_type::task_base;
    using queue_type     = typename scheduler_type::queue_type;

    const int num_worker_threads     = Kokkos::Experimental::HPX::concurrency();
    static task_base_type *const end = (task_base_type *)task_base_type::EndTag;
    constexpr task_base_type *no_more_tasks_sentinel = nullptr;

    thread_buffer &buffer = Kokkos::Experimental::HPX().impl_get_buffer();
    buffer.resize(num_worker_threads, 512);

    auto &queue = scheduler->queue();
    queue.initialize_team_queues(num_worker_threads);

    auto exec = Kokkos::Experimental::HPX::impl_get_executor();

    for_loop(
        par.on(exec).with(static_chunk_size(1)), 0, num_worker_threads,
        [this, &buffer, num_worker_threads](int t) {
          // NOTE: This implementation has been simplified based on the
          // assumption that team_size = 1. The HPX backend currently only
          // supports a team size of 1.
          buffer.get(Kokkos::Experimental::HPX::impl_hardware_thread_id());
          HPXTeamMember member(
              TeamPolicyInternal<Kokkos::Experimental::HPX>(
                  Kokkos::Experimental::HPX(), num_worker_threads, 1),
              0, t, buffer.get(t), 512);

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
        });
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
