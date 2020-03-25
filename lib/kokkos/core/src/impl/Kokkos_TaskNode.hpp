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

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_IMPL_TASKNODE_HPP
#define KOKKOS_IMPL_TASKNODE_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

#include <Kokkos_TaskScheduler_fwd.hpp>
#include <Kokkos_Core_fwd.hpp>

#include <Kokkos_PointerOwnership.hpp>

#include <impl/Kokkos_VLAEmulation.hpp>
#include <impl/Kokkos_LIFO.hpp>
#include <impl/Kokkos_ChaseLev.hpp>
#include <impl/Kokkos_EBO.hpp>
#include <Kokkos_Concepts.hpp>

#include <string>
#include <typeinfo>
#include <stdexcept>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_COMPILER_PGI
// Bizzarely, an extra jump instruction forces the PGI compiler to not have a
// bug related to (probably?) empty base optimization and/or aggregate
// construction.  This must be defined out-of-line to generate a jump
// jump instruction
void _kokkos_pgi_compiler_bug_workaround();
#endif

enum TaskType : int16_t {
  TaskTeam    = 0,
  TaskSingle  = 1,
  Aggregate   = 2,
  TaskSpecial = -1
};

//==============================================================================

/** Intrusive base class for things allocated with a Kokkos::MemoryPool
 *
 *  @warning Memory pools assume that the address of this class is the same
 *           as the address of the most derived type that was allocated to
 *           have the given size.  As a consequence, when interacting with
 *           multiple inheritance, this must always be the first base class
 *           of any derived class that uses it!
 *  @todo Consider inverting inheritance structure to avoid this problem?
 *
 *  @tparam CountType type of integer used to store the allocation size
 */
template <class CountType = int32_t>
class alignas(void*) PoolAllocatedObjectBase {
 public:
  using pool_allocation_size_type = CountType;

 private:
  pool_allocation_size_type m_alloc_size;

 public:
  KOKKOS_INLINE_FUNCTION
  constexpr explicit PoolAllocatedObjectBase(
      pool_allocation_size_type allocation_size)
      : m_alloc_size(allocation_size) {}

  KOKKOS_INLINE_FUNCTION
  CountType get_allocation_size() const noexcept { return m_alloc_size; }
};

//==============================================================================

// TODO @tasking @cleanup DSH move this?
template <class CountType = int32_t>
class ReferenceCountedBase {
 public:
  using reference_count_size_type = CountType;

 private:
  reference_count_size_type m_ref_count = 0;

 public:
  KOKKOS_INLINE_FUNCTION
#ifndef KOKKOS_COMPILER_PGI
  constexpr
#endif
      explicit ReferenceCountedBase(
          reference_count_size_type initial_reference_count)
      : m_ref_count(initial_reference_count) {
    // This can't be here because it breaks constexpr
    // KOKKOS_EXPECTS(initial_reference_count > 0);
#ifdef KOKKOS_COMPILER_PGI
    Impl::_kokkos_pgi_compiler_bug_workaround();
#endif
  }

  /** Decrement the reference count,
   *  and return true iff this decrement caused
   *  the reference count to become zero
   */
  KOKKOS_INLINE_FUNCTION
  bool decrement_and_check_reference_count() {
    // TODO @tasking @memory_order DSH memory order
    auto old_count = Kokkos::atomic_fetch_add(&m_ref_count, -1);

    KOKKOS_ASSERT(old_count > 0 && "reference count greater less than zero!");

    return (old_count == 1);
  }

  KOKKOS_INLINE_FUNCTION
  void increment_reference_count() { Kokkos::atomic_increment(&m_ref_count); }
};

template <class TaskQueueTraits, class SchedulingInfo>
class AggregateTask;

template <class TaskQueueTraits>
class RunnableTaskBase;

//==============================================================================

template <class TaskQueueTraits>
class TaskNode
    : public PoolAllocatedObjectBase<int32_t>,  // size 4, must be first!
      public ReferenceCountedBase<int32_t>,     // size 4
      public TaskQueueTraits::template intrusive_task_base_type<
          TaskNode<TaskQueueTraits>>  // size 8+
{
 public:
  using priority_type = int16_t;

 private:
  using task_base_type              = TaskNode<TaskQueueTraits>;
  using pool_allocated_base_type    = PoolAllocatedObjectBase<int32_t>;
  using reference_counted_base_type = ReferenceCountedBase<int32_t>;
  using task_queue_traits           = TaskQueueTraits;
  using waiting_queue_type =
      typename task_queue_traits::template waiting_queue_type<TaskNode>;

  waiting_queue_type m_wait_queue;  // size 8+

  // TODO @tasking @cleanup DSH eliminate this, or make its purpose a bit more
  // clear.  It's only used in BasicFuture, and only for deallocation purposes
  TaskQueueBase* m_ready_queue_base;

  TaskType m_task_type;      // size 2
  priority_type m_priority;  // size 2
  bool m_is_respawning = false;

 public:
  KOKKOS_INLINE_FUNCTION
  constexpr TaskNode(TaskType task_type, TaskPriority priority,
                     TaskQueueBase* queue_base,
                     reference_count_size_type initial_reference_count,
                     pool_allocation_size_type allocation_size)
      : pool_allocated_base_type(
            /* allocation_size = */ allocation_size),
        reference_counted_base_type(
            /* initial_reference_count = */ initial_reference_count),
        m_wait_queue(),
        m_ready_queue_base(queue_base),
        m_task_type(task_type),
        m_priority(static_cast<priority_type>(priority)),
        m_is_respawning(false) {}

  TaskNode()                = delete;
  TaskNode(TaskNode const&) = delete;
  TaskNode(TaskNode&&)      = delete;
  TaskNode& operator=(TaskNode const&) = delete;
  TaskNode& operator=(TaskNode&&) = delete;

  KOKKOS_INLINE_FUNCTION
  bool is_aggregate() const noexcept {
    return m_task_type == TaskType::Aggregate;
  }

  KOKKOS_INLINE_FUNCTION
  bool is_runnable() const noexcept {
    return m_task_type != TaskType::Aggregate;
  }

  KOKKOS_INLINE_FUNCTION
  bool is_runnable() const volatile noexcept {
    return m_task_type != TaskType::Aggregate;
  }

  KOKKOS_INLINE_FUNCTION
  bool is_single_runnable() const noexcept {
    return m_task_type == TaskType::TaskSingle;
  }

  KOKKOS_INLINE_FUNCTION
  bool is_team_runnable() const noexcept {
    return m_task_type == TaskType::TaskTeam;
  }

  KOKKOS_INLINE_FUNCTION
  TaskType get_task_type() const noexcept { return m_task_type; }

  KOKKOS_INLINE_FUNCTION
  RunnableTaskBase<TaskQueueTraits>& as_runnable_task() & {
    KOKKOS_EXPECTS(this->is_runnable());
    return static_cast<RunnableTaskBase<TaskQueueTraits>&>(*this);
  }

  KOKKOS_INLINE_FUNCTION
  RunnableTaskBase<TaskQueueTraits> const& as_runnable_task() const& {
    KOKKOS_EXPECTS(this->is_runnable());
    return static_cast<RunnableTaskBase<TaskQueueTraits> const&>(*this);
  }

  KOKKOS_INLINE_FUNCTION
  RunnableTaskBase<TaskQueueTraits> volatile& as_runnable_task() volatile& {
    KOKKOS_EXPECTS(this->is_runnable());
    return static_cast<RunnableTaskBase<TaskQueueTraits> volatile&>(*this);
  }

  KOKKOS_INLINE_FUNCTION
  RunnableTaskBase<TaskQueueTraits> const volatile& as_runnable_task() const
      volatile& {
    KOKKOS_EXPECTS(this->is_runnable());
    return static_cast<RunnableTaskBase<TaskQueueTraits> const volatile&>(
        *this);
  }

  KOKKOS_INLINE_FUNCTION
  RunnableTaskBase<TaskQueueTraits>&& as_runnable_task() && {
    KOKKOS_EXPECTS(this->is_runnable());
    return static_cast<RunnableTaskBase<TaskQueueTraits>&&>(*this);
  }

  template <class SchedulingInfo>
  KOKKOS_INLINE_FUNCTION AggregateTask<TaskQueueTraits, SchedulingInfo>&
  as_aggregate() & {
    KOKKOS_EXPECTS(this->is_aggregate());
    return static_cast<AggregateTask<TaskQueueTraits, SchedulingInfo>&>(*this);
  }

  template <class SchedulingInfo>
  KOKKOS_INLINE_FUNCTION AggregateTask<TaskQueueTraits, SchedulingInfo> const&
  as_aggregate() const& {
    KOKKOS_EXPECTS(this->is_aggregate());
    return static_cast<AggregateTask<TaskQueueTraits, SchedulingInfo> const&>(
        *this);
  }

  template <class SchedulingInfo>
  KOKKOS_INLINE_FUNCTION AggregateTask<TaskQueueTraits, SchedulingInfo>&&
  as_aggregate() && {
    KOKKOS_EXPECTS(this->is_aggregate());
    return static_cast<AggregateTask<TaskQueueTraits, SchedulingInfo>&&>(*this);
  }

  KOKKOS_INLINE_FUNCTION
  bool try_add_waiting(task_base_type& depends_on_this) {
    return m_wait_queue.try_push(depends_on_this);
  }

  template <class Function>
  KOKKOS_INLINE_FUNCTION void consume_wait_queue(Function&& f) {
    KOKKOS_EXPECTS(not m_wait_queue.is_consumed());
    m_wait_queue.consume(std::forward<Function>(f));
  }

  KOKKOS_INLINE_FUNCTION
  bool wait_queue_is_consumed() const noexcept {
    // TODO @tasking @memory_order DSH memory order
    return m_wait_queue.is_consumed();
  }

  KOKKOS_INLINE_FUNCTION
  TaskQueueBase* ready_queue_base_ptr() const noexcept {
    return m_ready_queue_base;
  }

  KOKKOS_INLINE_FUNCTION
  void set_priority(TaskPriority priority) noexcept {
    KOKKOS_EXPECTS(!this->is_enqueued());
    m_priority = (priority_type)priority;
  }

  KOKKOS_INLINE_FUNCTION
  void set_priority(TaskPriority priority) volatile noexcept {
    KOKKOS_EXPECTS(!this->is_enqueued());
    m_priority = (priority_type)priority;
  }

  KOKKOS_INLINE_FUNCTION
  TaskPriority get_priority() const noexcept {
    return (TaskPriority)m_priority;
  }

  KOKKOS_INLINE_FUNCTION
  bool get_respawn_flag() const { return m_is_respawning; }

  KOKKOS_INLINE_FUNCTION
  void set_respawn_flag(bool value = true) { m_is_respawning = value; }

  KOKKOS_INLINE_FUNCTION
  void set_respawn_flag(bool value = true) volatile { m_is_respawning = value; }
};

//==============================================================================

template <class BaseClass, class SchedulingInfo>
class SchedulingInfoStorage;

//==============================================================================

template <class BaseType, class SchedulingInfo>
class SchedulingInfoStorage
    : public BaseType,  // must be first base class for allocation reasons!!!
      private NoUniqueAddressMemberEmulation<SchedulingInfo> {
 private:
  using base_t                    = BaseType;
  using task_scheduling_info_type = SchedulingInfo;

 public:
  // Can't just do using base_t::base_t because of stupid stuff with clang cuda
  template <class... Args>
  // requires std::is_constructible_v<base_t, Args&&...>
  KOKKOS_INLINE_FUNCTION constexpr explicit SchedulingInfoStorage(
      Args&&... args)
      : base_t(std::forward<Args>(args)...) {}

  KOKKOS_INLINE_FUNCTION
  task_scheduling_info_type& scheduling_info() & {
    return this->no_unique_address_data_member();
  }

  KOKKOS_INLINE_FUNCTION
  task_scheduling_info_type const& scheduling_info() const& {
    return this->no_unique_address_data_member();
  }

  KOKKOS_INLINE_FUNCTION
  task_scheduling_info_type&& scheduling_info() && {
    return std::move(*this).no_unique_address_data_member();
  }
};

//==============================================================================

template <class TaskQueueTraits, class SchedulingInfo>
class alignas(16) AggregateTask final
    : public SchedulingInfoStorage<TaskNode<TaskQueueTraits>,
                                   SchedulingInfo>,  // must be first base class
                                                     // for allocation
                                                     // reasons!!!
      public ObjectWithVLAEmulation<
          AggregateTask<TaskQueueTraits, SchedulingInfo>,
          OwningRawPtr<TaskNode<TaskQueueTraits>>> {
 private:
  using base_t =
      SchedulingInfoStorage<TaskNode<TaskQueueTraits>, SchedulingInfo>;
  using vla_base_t =
      ObjectWithVLAEmulation<AggregateTask<TaskQueueTraits, SchedulingInfo>,
                             OwningRawPtr<TaskNode<TaskQueueTraits>>>;

  using task_base_type = TaskNode<TaskQueueTraits>;

 public:
  using aggregate_task_type = AggregateTask;  // concept marker

  template <class... Args>
  // requires std::is_constructible_v<base_t, Args&&...>
  KOKKOS_INLINE_FUNCTION constexpr explicit AggregateTask(
      int32_t aggregate_predecessor_count, Args&&... args)
      : base_t(TaskType::Aggregate,
               TaskPriority::Regular,  // all aggregates are regular priority
               std::forward<Args>(args)...),
        vla_base_t(aggregate_predecessor_count) {}

  KOKKOS_INLINE_FUNCTION
  int32_t dependence_count() const { return this->n_vla_entries(); }
};

// KOKKOS_IMPL_IS_CONCEPT(aggregate_task);

//==============================================================================

template <class TaskQueueTraits>
class RunnableTaskBase
    : public TaskNode<TaskQueueTraits>  // must be first base class for
                                        // allocation reasons!!!
{
 private:
  using base_t = TaskNode<TaskQueueTraits>;

 public:
  using task_base_type     = TaskNode<TaskQueueTraits>;
  using function_type      = void (*)(task_base_type*, void*);
  using destroy_type       = void (*)(task_base_type*);
  using runnable_task_type = RunnableTaskBase;

 private:
  function_type m_apply;
  task_base_type* m_predecessor = nullptr;

 public:
  template <class... Args>
  // requires std::is_constructible_v<base_t, Args&&...>
  KOKKOS_INLINE_FUNCTION constexpr explicit RunnableTaskBase(
      function_type apply_function_ptr, Args&&... args)
      : base_t(std::forward<Args>(args)...), m_apply(apply_function_ptr) {}

  KOKKOS_INLINE_FUNCTION
  bool has_predecessor() const { return m_predecessor != nullptr; }

  KOKKOS_INLINE_FUNCTION
  void clear_predecessor() { m_predecessor = nullptr; }

  KOKKOS_INLINE_FUNCTION
  void clear_predecessor() volatile { m_predecessor = nullptr; }

  template <class SchedulingInfo>
  KOKKOS_INLINE_FUNCTION SchedulingInfo& scheduling_info_as() {
    using info_storage_type =
        SchedulingInfoStorage<RunnableTaskBase, SchedulingInfo>;

    return static_cast<info_storage_type*>(this)->scheduling_info();
  }

  template <class SchedulingInfo>
  KOKKOS_INLINE_FUNCTION SchedulingInfo const& scheduling_info_as() const {
    using info_storage_type =
        SchedulingInfoStorage<RunnableTaskBase, SchedulingInfo>;

    return static_cast<info_storage_type const*>(this)->scheduling_info();
  }

  KOKKOS_INLINE_FUNCTION
  task_base_type& get_predecessor() const {
    KOKKOS_EXPECTS(m_predecessor != nullptr);
    return *m_predecessor;
  }

  KOKKOS_INLINE_FUNCTION
  void set_predecessor(task_base_type& predecessor) {
    KOKKOS_EXPECTS(m_predecessor == nullptr);
    // Increment the reference count so that predecessor doesn't go away
    // before this task is enqueued.
    // (should be memory order acquire)
    predecessor.increment_reference_count();
    m_predecessor = &predecessor;
  }

  KOKKOS_INLINE_FUNCTION
  void acquire_predecessor_from(runnable_task_type& other) {
    KOKKOS_EXPECTS(m_predecessor == nullptr ||
                   other.m_predecessor == m_predecessor);
    // since we're transferring, no need to modify the reference count
    m_predecessor       = other.m_predecessor;
    other.m_predecessor = nullptr;
  }

  KOKKOS_INLINE_FUNCTION
  void acquire_predecessor_from(runnable_task_type& other) volatile {
    KOKKOS_EXPECTS(m_predecessor == nullptr ||
                   other.m_predecessor == m_predecessor);
    // since we're transferring, no need to modify the reference count
    m_predecessor       = other.m_predecessor;
    other.m_predecessor = nullptr;
  }

  template <class TeamMember>
  KOKKOS_INLINE_FUNCTION void run(TeamMember& member) {
    (*m_apply)(this, &member);
  }
};

// KOKKOS_IMPL_IS_CONCEPT(runnable_task);

//==============================================================================

template <class ResultType, class Base>
class TaskResultStorage : public Base {
 private:
  using base_t = Base;

  alignas(Base) ResultType m_value = ResultType{};

 public:
  // using base_t::base_t;
  // Can't just do using base_t::base_t because of stupid stuff with clang cuda
  template <class... Args>
  // requires std::is_constructible_v<base_t, Args&&...>
  KOKKOS_INLINE_FUNCTION constexpr explicit TaskResultStorage(Args&&... args)
      : base_t(std::forward<Args>(args)...) {}

  KOKKOS_INLINE_FUNCTION
  ResultType* value_pointer() {
    // Over-alignment makes this a non-standard-layout class,
    // so alignas() doesn't work
    // static_assert(
    //  offsetof(TaskResultStorage, m_value) == sizeof(Base),
    //  "TaskResultStorage must be POD for layout purposes"
    //);
    return &m_value;
  }

  KOKKOS_INLINE_FUNCTION
  ResultType& value_reference() { return m_value; }
};

// TODO @tasking @optimization DSH optimization for empty types (in addition to
// void)
template <class Base>
class TaskResultStorage<void, Base> : public Base {
 private:
  using base_t = Base;

 public:
  // using base_t::base_t;
  // Can't just do using base_t::base_t because of stupid stuff with clang cuda
  template <class... Args>
  // requires std::is_constructible_v<base_t, Args&&...>
  KOKKOS_INLINE_FUNCTION constexpr explicit TaskResultStorage(Args&&... args)
      : base_t(std::forward<Args>(args)...) {}

  KOKKOS_INLINE_FUNCTION
  void* value_pointer() noexcept { return nullptr; }

  KOKKOS_INLINE_FUNCTION
  void value_reference() noexcept {}
};

//==============================================================================

template <class TaskQueueTraits, class Scheduler, class ResultType,
          class FunctorType>
class alignas(16) RunnableTask
    :  // using nesting of base classes to control layout; multiple empty base
       // classes may not be ABI compatible with CUDA on Windows
       public TaskResultStorage<
           ResultType,
           SchedulingInfoStorage<RunnableTaskBase<TaskQueueTraits>,
                                 typename Scheduler::task_queue_type::
                                     task_scheduling_info_type>>,  // must be
                                                                   // first base
                                                                   // class
       public FunctorType {
 private:
  using base_t = TaskResultStorage<
      ResultType,
      SchedulingInfoStorage<
          RunnableTaskBase<TaskQueueTraits>,
          typename Scheduler::task_queue_type::task_scheduling_info_type>>;

  using runnable_task_base_type = RunnableTaskBase<TaskQueueTraits>;
  using scheduler_type          = Scheduler;
  using scheduling_info_type =
      typename scheduler_type::task_scheduling_info_type;
  using scheduling_info_storage_base = base_t;

  using task_base_type = TaskNode<TaskQueueTraits>;
  using specialization = TaskQueueSpecialization<scheduler_type>;
  using member_type    = typename specialization::member_type;
  using result_type    = ResultType;
  using functor_type   = FunctorType;

 public:
  template <class... Args>
  // requires std::is_constructible_v<base_t, Args&&...>
  KOKKOS_INLINE_FUNCTION constexpr explicit RunnableTask(FunctorType&& functor,
                                                         Args&&... args)
      : base_t(std::forward<Args>(args)...), functor_type(std::move(functor)) {}

  KOKKOS_INLINE_FUNCTION
  ~RunnableTask() = delete;

  KOKKOS_INLINE_FUNCTION
  void update_scheduling_info(member_type& member) {
    // TODO @tasking @generalization DSH call a queue-specific hook here; for
    // now, this info is already updated elsewhere this->scheduling_info() =
    // member.scheduler().scheduling_info();
  }

  KOKKOS_INLINE_FUNCTION
  void apply_functor(member_type* member, void*) {
    update_scheduling_info(*member);
    this->functor_type::operator()(*member);
  }

  template <typename T>
  KOKKOS_INLINE_FUNCTION void apply_functor(member_type* member, T* val) {
    update_scheduling_info(*member);
    // this->functor_type::operator()(*member, *val);
    this->functor_type::operator()(*member, *val);
  }

  KOKKOS_FUNCTION static void destroy(task_base_type* root) {
    // TaskResult<result_type>::destroy(root);
  }

  KOKKOS_FUNCTION static void apply(task_base_type* self,
                                    void* member_as_void) {
    using task_type = Impl::RunnableTask<TaskQueueTraits, Scheduler, ResultType,
                                         FunctorType>*;
    auto* const task   = static_cast<task_type>(self);
    auto* const member = reinterpret_cast<member_type*>(member_as_void);

    // Now that we're over-aligning the result storage, this isn't a problem any
    // more
    // static_assert(std::is_standard_layout<task_type>::value,
    //  "Tasks must be standard layout"
    //);
    // static_assert(std::is_pod<task_type>::value,
    //  "Tasks must be PODs"
    //);

    // Task may be serial or team.
    // If team then must synchronize before querying if respawn was requested.
    // If team then only one thread calls destructor.

    const bool only_one_thread =
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
        0 == threadIdx.x && 0 == threadIdx.y;
#else
        0 == member->team_rank();
#endif

    // Ensure that the respawn flag is set to zero
    self->set_respawn_flag(false);

    // task->apply_functor(member, TaskResult<result_type>::ptr(task));
    task->apply_functor(member, task->value_pointer());

    member->team_barrier();

    if (only_one_thread && !(task->get_respawn_flag())) {
      // Did not respawn, destroy the functor to free memory.
      task->functor_type::~functor_type();
      // Cannot destroy and deallocate the task until its dependences
      // have been processed.
    }
  }
};

} /* namespace Impl */

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_TASKNODE_HPP */
