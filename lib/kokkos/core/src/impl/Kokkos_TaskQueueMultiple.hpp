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

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_IMPL_TASKQUEUEMULTIPLE_HPP
#define KOKKOS_IMPL_TASKQUEUEMULTIPLE_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_TaskScheduler_fwd.hpp>
#include <Kokkos_Core_fwd.hpp>

#include <Kokkos_MemoryPool.hpp>

#include <impl/Kokkos_TaskBase.hpp>
#include <impl/Kokkos_TaskResult.hpp>
#include <impl/Kokkos_TaskQueue.hpp>

#include <Kokkos_Atomic.hpp>

#include <string>
#include <typeinfo>
#include <cassert>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <typename ExecSpace,
          typename MemorySpace = typename ExecSpace::memory_space>
class LeagueQueueCollection;

template <class ExecSpace, class MemorySpace>
class TaskQueueMultiple : public TaskQueue<ExecSpace, MemorySpace> {
 private:
  using base_t             = TaskQueue<ExecSpace, MemorySpace>;
  using queue_collection_t = LeagueQueueCollection<ExecSpace, MemorySpace>;

  int m_league_rank = static_cast<int>(KOKKOS_INVALID_INDEX);

  // This pointer is owning only if m_league_rank == 0
  queue_collection_t* m_other_queues = nullptr;

 public:
  struct Destroy {
    TaskQueueMultiple* m_queue;
    void destroy_shared_allocation();
  };

  using team_queue_type = TaskQueueMultiple;

  TaskQueueMultiple(int arg_league_rank, queue_collection_t* arg_other_queues,
                    typename base_t::memory_pool const& arg_memory_pool)
      : base_t(arg_memory_pool),
        m_league_rank(arg_league_rank),
        m_other_queues(arg_other_queues) {}

  explicit TaskQueueMultiple(
      typename base_t::memory_pool const& arg_memory_pool)
      : base_t(arg_memory_pool), m_league_rank(0) {
    void* other_queues_buffer =
        typename base_t::memory_space{}.allocate(sizeof(queue_collection_t));
    m_other_queues = new (other_queues_buffer) queue_collection_t(this);
  }

  ~TaskQueueMultiple() {
    if (m_league_rank == 0 && m_other_queues != nullptr) {
      m_other_queues->~queue_collection_t();
      typename base_t::memory_space{}.deallocate(m_other_queues,
                                                 sizeof(queue_collection_t));
    }
    // rest of destruction is handled in the base class
  }

  //----------------------------------------

  void initialize_team_queues(int arg_league_size) const noexcept {
    m_other_queues->initialize_team_queues(arg_league_size, this->m_memory);
  }

  KOKKOS_INLINE_FUNCTION
  team_queue_type& get_team_queue(int arg_league_rank) noexcept {
    if (arg_league_rank == m_league_rank)
      return *this;
    else
      return m_other_queues->get_team_queue(arg_league_rank);
  }

  KOKKOS_INLINE_FUNCTION
  typename base_t::task_root_type* attempt_to_steal_task() noexcept {
    TaskBase* rv        = nullptr;
    auto* const end_tag = reinterpret_cast<TaskBase*>(TaskBase::EndTag);

    if (m_other_queues == nullptr) {
      Kokkos::abort("attempted to steal task before queues were initialized!");
    }

    // Loop by priority and then type, and then team
    for (int i = 0; i < base_t::NumQueue; ++i) {
      for (int j = 0; j < 2; ++j) {
        // for now, always start by trying to steal from team zero
        for (int iteam = 0; iteam < m_other_queues->size(); ++iteam) {
          if (iteam == m_league_rank) continue;
          auto& steal_from = get_team_queue(iteam);
          if (*((volatile int*)&steal_from.m_ready_count) > 0) {
            // we've found at least one queue that's not done, so even if we
            // can't pop something off of it we shouldn't return a nullptr
            // indicating completion.  rv will be end_tag when the pop fails
            rv = base_t::pop_ready_task(&steal_from.m_ready[i][j]);
            if (rv != end_tag) {
              // task stolen.
              // first increment our ready count, then decrement the ready count
              // on the other queue:
              desul::atomic_inc(
                  &this->m_ready_count, desul::MemoryOrderSeqCst(),
                  desul::MemoryScopeDevice());  // TODO? memory_order_relaxed
              desul::atomic_dec(
                  &steal_from.m_ready_count, desul::MemoryOrderSeqCst(),
                  desul::MemoryScopeDevice());  // TODO? memory_order_relaxed
              return rv;
            }
          }
        }
      }
    }

    // at this point, rv will only be nullptr if *all* of the queues had an
    // m_ready_count of 0.  This indicates quiescence.  If at least some of them
    // had non-zero, there would have been at least one pop_ready_task that
    // was called and returned end_tag if it couldn't pop a task
    return rv;
  }
};

template <typename ExecSpace, typename MemorySpace>
class LeagueQueueCollection {
 private:
  using execution_space     = ExecSpace;
  using memory_space        = MemorySpace;
  using device_type         = Kokkos::Device<execution_space, memory_space>;
  using memory_pool         = Kokkos::MemoryPool<device_type>;
  using team_queue_type     = TaskQueueMultiple<execution_space, memory_space>;
  using team_scheduler_type = BasicTaskScheduler<ExecSpace, team_queue_type>;
  using specialization      = TaskQueueSpecialization<team_scheduler_type>;

  enum : long { max_num_queues = 6 };  // specialization::max_league_size };

  // this is a non-owning pointer
  team_queue_type* m_rank_zero_queue = nullptr;
  // This really needs to be an optional<TaskQueue<ExecSpace>>
  union optional_queue {
    KOKKOS_INLINE_FUNCTION
    optional_queue() : uninitialized(0) {}
    KOKKOS_INLINE_FUNCTION
    ~optional_queue() { uninitialized = 0; }
    char uninitialized;
    team_queue_type initialized;
  } m_queues[max_num_queues];
  int m_size = static_cast<int>(KOKKOS_INVALID_INDEX);

 public:
  LeagueQueueCollection()                             = delete;
  LeagueQueueCollection(LeagueQueueCollection const&) = delete;
  LeagueQueueCollection(LeagueQueueCollection&&)      = delete;
  LeagueQueueCollection& operator=(LeagueQueueCollection const&) = delete;
  LeagueQueueCollection& operator=(LeagueQueueCollection&&) = delete;

  ~LeagueQueueCollection() {
    // destroy only the initialized queues that we own
    for (int iteam = 0; iteam < m_size - 1; ++iteam) {
      m_queues[iteam].initialized.~team_queue_type();
      m_queues[iteam].uninitialized = 0;
    }
  }

  KOKKOS_INLINE_FUNCTION
  explicit LeagueQueueCollection(team_queue_type* arg_rank_zero_queue)
      : m_rank_zero_queue(arg_rank_zero_queue), m_size(1) {}

  void initialize_team_queues(int arg_count,
                              memory_pool const& arg_memory_pool) noexcept {
    arg_count = std::min((int)max_num_queues, arg_count);
    // assert(arg_count <= max_num_queues);
    if (arg_count > m_size) {
      for (int i = m_size; i < arg_count; ++i) {
        new (&m_queues[i - 1].initialized)
            team_queue_type(i, this, arg_memory_pool);
      }
      m_size = arg_count;
    }
  }

  KOKKOS_INLINE_FUNCTION
  constexpr int size() const noexcept { return m_size; }

  KOKKOS_INLINE_FUNCTION
  constexpr bool initialized() const noexcept {
    return m_size != int(KOKKOS_INVALID_INDEX);
  }

  KOKKOS_INLINE_FUNCTION
  team_queue_type& get_team_queue(int iteam) {
    iteam %= max_num_queues;
#if !defined(__HIP_DEVICE_COMPILE__) && !defined(__CUDA_ARCH__)
    assert(initialized());
    assert(iteam < m_size);
    assert(iteam >= 0);
#endif
    if (iteam == 0)
      return *m_rank_zero_queue;
    else
      return m_queues[iteam - 1].initialized;
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <impl/Kokkos_TaskQueueMultiple_impl.hpp>

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_TASKQUEUEMULTIPLE_HPP */
