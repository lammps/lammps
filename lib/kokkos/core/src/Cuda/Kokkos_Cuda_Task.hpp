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

#ifndef KOKKOS_IMPL_CUDA_TASK_HPP
#define KOKKOS_IMPL_CUDA_TASK_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_TASKDAG)

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <Kokkos_Core_fwd.hpp>

#include <impl/Kokkos_TaskBase.hpp>
#include <Cuda/Kokkos_Cuda_Error.hpp>  // KOKKOS_IMPL_CUDA_SAFE_CALL
#include <impl/Kokkos_TaskTeamMember.hpp>

//----------------------------------------------------------------------------

#if defined(__CUDA_ARCH__)
#define KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN(MSG)                           \
  {                                                                        \
    __syncwarp();                                                          \
    const unsigned b = __activemask();                                     \
    if (b != 0xffffffff) {                                                 \
      printf(" SYNCWARP AT %s (%d,%d,%d) (%d,%d,%d) failed %x\n", MSG,     \
             blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, \
             threadIdx.z, b);                                              \
      return;                                                              \
    }                                                                      \
  }
#else
#define KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN(MSG)
#endif

namespace Kokkos {
namespace Impl {
namespace {

template <typename TaskType>
__global__ void set_cuda_task_base_apply_function_pointer(
    typename TaskType::function_type* ptr,
    typename TaskType::destroy_type* dtor) {
  *ptr  = TaskType::apply;
  *dtor = TaskType::destroy;
}

template <typename Scheduler>
__global__ void cuda_task_queue_execute(Scheduler scheduler,
                                        int32_t shmem_size) {
  TaskQueueSpecialization<Scheduler>::driver(std::move(scheduler), shmem_size);
}

}  // namespace

template <class, class>
class TaskExec;

template <class QueueType>
class TaskQueueSpecialization<SimpleTaskScheduler<Kokkos::Cuda, QueueType>> {
 public:
  using scheduler_type  = SimpleTaskScheduler<Kokkos::Cuda, QueueType>;
  using execution_space = Kokkos::Cuda;
  using memory_space    = Kokkos::CudaUVMSpace;
  using member_type     = TaskExec<Kokkos::Cuda, scheduler_type>;

  enum : long { max_league_size = 16 };
  enum : int { warps_per_block = 4 };

  KOKKOS_INLINE_FUNCTION
  static void iff_single_thread_recursive_execute(scheduler_type const&) {}

  static int get_max_team_count(execution_space const&) {
    return Kokkos::Impl::cuda_internal_multiprocessor_count() * warps_per_block;
  }

  __device__ static void driver(scheduler_type scheduler,
                                int32_t shmem_per_warp) {
    using queue_type     = typename scheduler_type::task_queue_type;
    using task_base_type = typename scheduler_type::task_base_type;
    using runnable_task_base_type =
        typename scheduler_type::runnable_task_base_type;
    using scheduling_info_storage_type = SchedulingInfoStorage<
        runnable_task_base_type,
        typename scheduler_type::task_scheduling_info_type>;

    extern __shared__ int32_t shmem_all[];

    int32_t* const warp_shmem =
        shmem_all + (threadIdx.z * shmem_per_warp) / sizeof(int32_t);

    task_base_type* const shared_memory_task_copy = (task_base_type*)warp_shmem;

    const int warp_lane = threadIdx.x + threadIdx.y * blockDim.x;

    member_type single_exec(scheduler, warp_shmem, 1);
    member_type team_exec(scheduler, warp_shmem, blockDim.y);

    auto& queue          = scheduler.queue();
    auto& team_scheduler = team_exec.scheduler();

    auto current_task = OptionalRef<task_base_type>();

    // Loop until all queues are empty and no tasks in flight
    while (!queue.is_done()) {
      if (warp_lane == 0) {  // should be (?) same as team_exec.team_rank() == 0
        // pop off a task
        current_task =
            queue.pop_ready_task(team_scheduler.team_scheduler_info());
      }

      // Broadcast task pointer:

      // Sync before the broadcast
      __syncwarp(0xffffffff);

      // pretend it's an int* for shuffle purposes
      ((int*)&current_task)[0] =
          __shfl_sync(0xffffffff, ((int*)&current_task)[0], 0, 32);
      ((int*)&current_task)[1] =
          __shfl_sync(0xffffffff, ((int*)&current_task)[1], 0, 32);

      if (current_task) {
        KOKKOS_ASSERT(!current_task->as_runnable_task().get_respawn_flag());

        int32_t b = sizeof(scheduling_info_storage_type) / sizeof(int32_t);
        static_assert(
            sizeof(scheduling_info_storage_type) % sizeof(int32_t) == 0,
            "bad task size");
        int32_t const e = current_task->get_allocation_size() / sizeof(int32_t);
        KOKKOS_ASSERT(current_task->get_allocation_size() % sizeof(int32_t) ==
                      0);

        int32_t volatile* const task_mem =
            (int32_t volatile*)current_task.get();

        // do a coordinated copy of the task closure from global to shared
        // memory:
        for (int32_t i = warp_lane; i < e; i += CudaTraits::WarpSize) {
          warp_shmem[i] = task_mem[i];
        }

        // Synchronize threads of the warp and insure memory
        // writes are visible to all threads in the warp.
        __syncwarp(0xffffffff);

        if (shared_memory_task_copy->is_team_runnable()) {
          // Thread Team Task
          shared_memory_task_copy->as_runnable_task().run(team_exec);
        } else if (threadIdx.y == 0) {
          // TODO @tasking @optimization DSH Change this to warp_lane == 0 when
          // we allow blockDim.x to be more than 1 Single Thread Task
          shared_memory_task_copy->as_runnable_task().run(single_exec);
        }

        // Synchronize threads of the warp and insure memory
        // writes are visible to all threads in the warp.

        __syncwarp(0xffffffff);

        // if(warp_lane < b % CudaTraits::WarpSize) b += CudaTraits::WarpSize;
        // b -= b % CudaTraits::WarpSize;

        // copy task closure from shared to global memory:
        for (int32_t i = b + warp_lane; i < e; i += CudaTraits::WarpSize) {
          task_mem[i] = warp_shmem[i];
        }

        // Synchronize threads of the warp and insure memory
        // writes are visible to root thread of the warp for
        // respawn or completion.

        __syncwarp(0xffffffff);

        if (warp_lane == 0) {
          // If respawn requested copy respawn data back to main memory
          if (shared_memory_task_copy->as_runnable_task().get_respawn_flag()) {
            if (shared_memory_task_copy->as_runnable_task().has_predecessor()) {
              // It's not necessary to make this a volatile write because
              // the next read of the predecessor is on this thread in complete,
              // and the predecessor is cleared there (using a volatile write)
              current_task->as_runnable_task().acquire_predecessor_from(
                  shared_memory_task_copy->as_runnable_task());
            }

            // It may not necessary to make this a volatile write, since the
            // next read will be done by this thread in complete where the
            // rescheduling occurs, but since the task could be stolen later
            // before this is written again, we should do the volatile write
            // here.  (It might not be necessary though because I don't know
            // where else the priority would be read after it is scheduled
            // by this thread; for now, we leave it volatile, but we should
            // benchmark the cost of this.)
            current_task.as_volatile()->set_priority(
                shared_memory_task_copy->get_priority());

            // It's not necessary to make this a volatile write, since the
            // next read of it (if true) will be by this thread in `complete()`,
            // which will unset the flag (using volatile) once it has handled
            // the respawn
            current_task->as_runnable_task().set_respawn_flag();
          }

          queue.complete((*std::move(current_task)).as_runnable_task(),
                         team_scheduler.team_scheduler_info());
        }
      }
    }
  }

  static void execute(scheduler_type const& scheduler) {
    const int shared_per_warp = 2048;
    const dim3 grid(Kokkos::Impl::cuda_internal_multiprocessor_count(), 1, 1);
    const dim3 block(1, Kokkos::Impl::CudaTraits::WarpSize, warps_per_block);
    const int shared_total    = shared_per_warp * warps_per_block;
    const cudaStream_t stream = nullptr;

    KOKKOS_ASSERT(
        static_cast<long>(grid.x * grid.y * grid.z * block.x * block.y *
                          block.z) ==
        static_cast<long>(get_max_team_count(scheduler.get_execution_space()) *
                          Kokkos::Impl::CudaTraits::WarpSize));

    auto& queue = scheduler.queue();

    Impl::cuda_device_synchronize(
        "Kokkos::Impl::TaskQueueSpecialization<SimpleTaskScheduler<Kokkos::"
        "Cuda>::execute: Pre Task Execution");

    // Query the stack size, in bytes:

    size_t previous_stack_size = 0;
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        cudaDeviceGetLimit(&previous_stack_size, cudaLimitStackSize));

    // If not large enough then set the stack size, in bytes:

    const size_t larger_stack_size = 1 << 11;

    if (previous_stack_size < larger_stack_size) {
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          cudaDeviceSetLimit(cudaLimitStackSize, larger_stack_size));
    }

    cuda_task_queue_execute<<<grid, block, shared_total, stream>>>(
        scheduler, shared_per_warp);

    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetLastError());

    Impl::cuda_device_synchronize(
        "Kokkos::Impl::TaskQueueSpecialization<SimpleTaskScheduler<Kokkos::"
        "Cuda>::execute: Post Task Execution");

    if (previous_stack_size < larger_stack_size) {
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          cudaDeviceSetLimit(cudaLimitStackSize, previous_stack_size));
    }
  }

  template <typename TaskType>
  static
      // TODO @tasking @optimiazation DSH specialize this for trivially
      // destructible types
      void
      get_function_pointer(typename TaskType::function_type& ptr,
                           typename TaskType::destroy_type& dtor) {
    using function_type = typename TaskType::function_type;
    using destroy_type  = typename TaskType::destroy_type;

    // TODO @tasking @minor DSH make sure there aren't any alignment concerns?
    void* storage = cuda_internal_scratch_unified(
        Kokkos::Cuda(), sizeof(function_type) + sizeof(destroy_type));
    function_type* ptr_ptr = (function_type*)storage;
    destroy_type* dtor_ptr =
        (destroy_type*)((char*)storage + sizeof(function_type));

    Impl::cuda_device_synchronize(
        "Kokkos::Impl::TaskQueueSpecialization<SimpleTaskScheduler<Kokkos::"
        "Cuda>::execute: Pre Get Function Pointer for Tasks");

    set_cuda_task_base_apply_function_pointer<TaskType>
        <<<1, 1>>>(ptr_ptr, dtor_ptr);

    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetLastError());
    Impl::cuda_device_synchronize(
        "Kokkos::Impl::TaskQueueSpecialization<SimpleTaskScheduler<Kokkos::"
        "Cuda>::execute: Post Get Function Pointer for Tasks");

    ptr  = *ptr_ptr;
    dtor = *dtor_ptr;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class Scheduler>
class TaskQueueSpecializationConstrained<
    Scheduler, std::enable_if_t<std::is_same<
                   typename Scheduler::execution_space, Kokkos::Cuda>::value>> {
 public:
  using scheduler_type  = Scheduler;
  using execution_space = Kokkos::Cuda;
  using memory_space    = Kokkos::CudaUVMSpace;
  using member_type     = TaskExec<Kokkos::Cuda, Scheduler>;

  enum : long { max_league_size = 16 };

  KOKKOS_INLINE_FUNCTION
  static void iff_single_thread_recursive_execute(scheduler_type const&) {}

  __device__ static void driver(scheduler_type scheduler,
                                int32_t shmem_per_warp) {
    using queue_type     = typename scheduler_type::queue_type;
    using task_root_type = TaskBase;

    extern __shared__ int32_t shmem_all[];

    task_root_type* const end = (task_root_type*)task_root_type::EndTag;
    task_root_type* const no_more_tasks_sentinel = nullptr;

    int32_t* const warp_shmem =
        shmem_all + (threadIdx.z * shmem_per_warp) / sizeof(int32_t);

    task_root_type* const task_shmem = (task_root_type*)warp_shmem;

    const int warp_lane = threadIdx.x + threadIdx.y * blockDim.x;

    member_type single_exec(scheduler, warp_shmem, 1);
    member_type team_exec(scheduler, warp_shmem, blockDim.y);

    auto& team_queue = team_exec.scheduler().queue();

    task_root_type* task_ptr = no_more_tasks_sentinel;

    // Loop until all queues are empty and no tasks in flight

    do {
      // Each team lead attempts to acquire either a thread team task
      // or collection of single thread tasks for the team.

      if (0 == warp_lane) {
        if (*((volatile int*)&team_queue.m_ready_count) > 0) {
          task_ptr = end;
          // Attempt to acquire a task
          // Loop by priority and then type
          for (int i = 0; i < queue_type::NumQueue && end == task_ptr; ++i) {
            for (int j = 0; j < 2 && end == task_ptr; ++j) {
              task_ptr = queue_type::pop_ready_task(&team_queue.m_ready[i][j]);
            }
          }
        } else {
          // returns nullptr if and only if all other queues have a ready
          // count of 0 also. Otherwise, returns a task from another queue
          // or `end` if one couldn't be popped
          task_ptr = team_queue.attempt_to_steal_task();
        }
      }

      // Synchronize warp with memory fence before broadcasting task pointer:

      // KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN( "A" );
      __syncwarp(0xffffffff);

      // Broadcast task pointer:

      ((int*)&task_ptr)[0] =
          __shfl_sync(0xffffffff, ((int*)&task_ptr)[0], 0, 32);
      ((int*)&task_ptr)[1] =
          __shfl_sync(0xffffffff, ((int*)&task_ptr)[1], 0, 32);

#if defined(KOKKOS_ENABLE_DEBUG)
      KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN("TaskQueue CUDA task_ptr");
#endif

      if (0 == task_ptr) break;  // 0 == queue->m_ready_count

      if (end != task_ptr) {
        // Whole warp copy task's closure to/from shared memory.
        // Use all threads of warp for coalesced read/write.

        int32_t const b = sizeof(task_root_type) / sizeof(int32_t);
        int32_t const e =
            *((int32_t volatile*)(&task_ptr->m_alloc_size)) / sizeof(int32_t);

        int32_t volatile* const task_mem = (int32_t volatile*)task_ptr;

        KOKKOS_ASSERT(e * sizeof(int32_t) < shmem_per_warp);

        // copy task closure from global to shared memory:

        for (int32_t i = warp_lane; i < e; i += CudaTraits::WarpSize) {
          warp_shmem[i] = task_mem[i];
        }

        // Synchronize threads of the warp and insure memory
        // writes are visible to all threads in the warp.

        // KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN( "B" );
        __syncwarp(0xffffffff);

        if (task_root_type::TaskTeam == task_shmem->m_task_type) {
          // Thread Team Task
          (*task_shmem->m_apply)(task_shmem, &team_exec);
        } else if (0 == threadIdx.y) {
          // Single Thread Task
          (*task_shmem->m_apply)(task_shmem, &single_exec);
        }

        // Synchronize threads of the warp and insure memory
        // writes are visible to all threads in the warp.

        // KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN( "C" );
        __syncwarp(0xffffffff);

        // copy task closure from shared to global memory:

        for (int32_t i = b + warp_lane; i < e; i += CudaTraits::WarpSize) {
          task_mem[i] = warp_shmem[i];
        }

        // Synchronize threads of the warp and insure memory
        // writes are visible to root thread of the warp for
        // respawn or completion.

        // KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN( "D" );
        __syncwarp(0xffffffff);

        // If respawn requested copy respawn data back to main memory

        if (0 == warp_lane) {
          if (((task_root_type*)task_root_type::LockTag) !=
              task_shmem->m_next) {
            ((volatile task_root_type*)task_ptr)->m_next = task_shmem->m_next;
            ((volatile task_root_type*)task_ptr)->m_priority =
                task_shmem->m_priority;
          }

          team_queue.complete(task_ptr);
        }
      }
    } while (1);
  }

  static void execute(scheduler_type const& scheduler) {
    const int shared_per_warp = 2048;
    const int warps_per_block = 4;
    const dim3 grid(Kokkos::Impl::cuda_internal_multiprocessor_count(), 1, 1);
    // const dim3 grid( 1 , 1 , 1 );
    const dim3 block(1, Kokkos::Impl::CudaTraits::WarpSize, warps_per_block);
    const int shared_total    = shared_per_warp * warps_per_block;
    const cudaStream_t stream = 0;

    auto& queue = scheduler.queue();
    queue.initialize_team_queues(warps_per_block * grid.x);

    Impl::cuda_device_synchronize(
        "Kokkos::Impl::TaskQueueSpecializationConstrained<SimpleTaskScheduler<"
        "Kokkos::Cuda>::execute: Pre Execute Task");

    // Query the stack size, in bytes:

    size_t previous_stack_size = 0;
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        cudaDeviceGetLimit(&previous_stack_size, cudaLimitStackSize));

    // If not large enough then set the stack size, in bytes:

    const size_t larger_stack_size = 2048;

    if (previous_stack_size < larger_stack_size) {
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          cudaDeviceSetLimit(cudaLimitStackSize, larger_stack_size));
    }

    cuda_task_queue_execute<<<grid, block, shared_total, stream>>>(
        scheduler, shared_per_warp);

    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetLastError());

    Impl::cuda_device_synchronize(
        "Kokkos::Impl::TaskQueueSpecializationConstrained<SimpleTaskScheduler<"
        "Kokkos::Cuda>::execute: Post Execute Task");

    if (previous_stack_size < larger_stack_size) {
      KOKKOS_IMPL_CUDA_SAFE_CALL(
          cudaDeviceSetLimit(cudaLimitStackSize, previous_stack_size));
    }
  }

  template <typename TaskType>
  static void get_function_pointer(typename TaskType::function_type& ptr,
                                   typename TaskType::destroy_type& dtor) {
    using function_type = typename TaskType::function_type;
    using destroy_type  = typename TaskType::destroy_type;

    void* storage = cuda_internal_scratch_unified(
        Kokkos::Cuda(), sizeof(function_type) + sizeof(destroy_type));
    function_type* ptr_ptr = (function_type*)storage;
    destroy_type* dtor_ptr =
        (destroy_type*)((char*)storage + sizeof(function_type));

    Impl::cuda_device_synchronize(
        "Kokkos::Impl::TaskQueueSpecializationConstrained<SimpleTaskScheduler<"
        "Kokkos::Cuda>::get_function_pointer: Pre Get Function Pointer");

    set_cuda_task_base_apply_function_pointer<TaskType>
        <<<1, 1>>>(ptr_ptr, dtor_ptr);

    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetLastError());
    Impl::cuda_device_synchronize(
        "Kokkos::Impl::TaskQueueSpecializationConstrained<SimpleTaskScheduler<"
        "Kokkos::Cuda>::get_function_pointer: Post Get Function Pointer");

    ptr  = *ptr_ptr;
    dtor = *dtor_ptr;
  }
};

extern template class TaskQueue<
    Kokkos::Cuda,
    default_tasking_memory_space_for_execution_space_t<Kokkos::Cuda>>;

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/**\brief  Impl::TaskExec<Cuda> is the TaskScheduler<Cuda>::member_type
 *         passed to tasks running in a Cuda space.
 *
 *  Cuda thread blocks for tasking are dimensioned:
 *    blockDim.x == vector length
 *    blockDim.y == team size
 *    blockDim.z == number of teams
 *  where
 *    blockDim.x * blockDim.y == WarpSize
 *
 *  Current implementation requires blockDim.x == 1.
 *  Vector level parallelism with blockDim.y > 1 on Volta will
 *  require a vector-level synchronization mask for vector-level
 *  collective operaitons.
 *
 *  Both single thread and thread team tasks are run by a full Cuda warp.
 *  A single thread task is called by warp lane #0 and the remaining
 *  lanes of the warp are idle.
 *
 *  When executing a single thread task the syncwarp or other
 *  warp synchronizing functions must not be called.
 */
template <class Scheduler>
class TaskExec<Kokkos::Cuda, Scheduler> {
 private:
  enum : int { WarpSize = Kokkos::Impl::CudaTraits::WarpSize };

  TaskExec(TaskExec&&)      = delete;
  TaskExec(TaskExec const&) = delete;
  TaskExec& operator=(TaskExec&&) = delete;
  TaskExec& operator=(TaskExec const&) = delete;

  friend class Kokkos::Impl::TaskQueue<
      Kokkos::Cuda,
      default_tasking_memory_space_for_execution_space_t<Kokkos::Cuda>>;
  template <class, class>
  friend class Kokkos::Impl::TaskQueueSpecializationConstrained;
  template <class>
  friend class Kokkos::Impl::TaskQueueSpecialization;

  int32_t* m_team_shmem;
  const int m_team_size;
  Scheduler m_scheduler;

  // If constructed with arg_team_size == 1 the object
  // can only be used by 0 == threadIdx.y.
  KOKKOS_INLINE_FUNCTION
  TaskExec(Scheduler const& parent_scheduler, int32_t* arg_team_shmem,
           int arg_team_size = blockDim.y)
      : m_team_shmem(arg_team_shmem),
        m_team_size(arg_team_size),
        m_scheduler(parent_scheduler.get_team_scheduler(league_rank())) {}

 public:
  using thread_team_member = TaskExec;

#if defined(__CUDA_ARCH__)
  __device__ int team_rank() const { return threadIdx.y; }
  __device__ int team_size() const { return m_team_size; }
  //__device__ int league_rank() const { return threadIdx.z; }
  __device__ int league_rank() const {
    return blockIdx.x * blockDim.z + threadIdx.z;
  }
  __device__ int league_size() const { return blockDim.z * gridDim.x; }

  __device__ void team_barrier() const {
    if (1 < m_team_size) {
      __syncwarp(0xffffffff);
    }
  }

  template <class ValueType>
  __device__ void team_broadcast(ValueType& val, const int thread_id) const {
    if (1 < m_team_size) {
      // WarpSize = blockDim.X * blockDim.y
      // thread_id < blockDim.y
      ValueType tmp(val);  // input might not be register variable
      Impl::in_place_shfl(val, tmp, blockDim.x * thread_id, WarpSize);
    }
  }

#else
  __host__ int team_rank() const { return 0; }
  __host__ int team_size() const { return 0; }
  __host__ int league_rank() const { return 0; }
  __host__ int league_size() const { return 0; }
  __host__ void team_barrier() const {}
  template <class ValueType>
  __host__ void team_broadcast(ValueType&, const int) const {}
#endif

  KOKKOS_INLINE_FUNCTION Scheduler const& scheduler() const noexcept {
    return m_scheduler;
  }
  KOKKOS_INLINE_FUNCTION Scheduler& scheduler() noexcept { return m_scheduler; }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <typename iType, typename Scheduler>
struct TeamThreadRangeBoundariesStruct<iType,
                                       TaskExec<Kokkos::Cuda, Scheduler>> {
  using index_type  = iType;
  using member_type = TaskExec<Kokkos::Cuda, Scheduler>;

  const iType start;
  const iType end;
  const iType increment;
  member_type const& thread;

#if defined(__CUDA_ARCH__)

  __device__ inline TeamThreadRangeBoundariesStruct(
      member_type const& arg_thread, const iType& arg_count)
      : start(threadIdx.y),
        end(arg_count),
        increment(blockDim.y),
        thread(arg_thread) {}

  __device__ inline TeamThreadRangeBoundariesStruct(
      member_type const& arg_thread, const iType& arg_start,
      const iType& arg_end)
      : start(arg_start + threadIdx.y),
        end(arg_end),
        increment(blockDim.y),
        thread(arg_thread) {}

#else

  TeamThreadRangeBoundariesStruct(member_type const& arg_thread,
                                  const iType& arg_count);

  TeamThreadRangeBoundariesStruct(member_type const& arg_thread,
                                  const iType& arg_start, const iType& arg_end);

#endif
};

//----------------------------------------------------------------------------

template <typename iType, typename Scheduler>
struct ThreadVectorRangeBoundariesStruct<iType,
                                         TaskExec<Kokkos::Cuda, Scheduler>> {
  using index_type  = iType;
  using member_type = TaskExec<Kokkos::Cuda, Scheduler>;

  const index_type start;
  const index_type end;
  const index_type increment;
  const member_type& thread;

#if defined(__CUDA_ARCH__)

  __device__ inline ThreadVectorRangeBoundariesStruct(
      member_type const& arg_thread, const index_type& arg_count)
      : start(threadIdx.x),
        end(arg_count),
        increment(blockDim.x),
        thread(arg_thread) {}

  __device__ inline ThreadVectorRangeBoundariesStruct(
      member_type const& arg_thread, const index_type& arg_begin,
      const index_type& arg_end)
      : start(arg_begin + threadIdx.x),
        end(arg_end),
        increment(blockDim.x),
        thread(arg_thread) {}

#else

  ThreadVectorRangeBoundariesStruct(member_type const& arg_thread,
                                    const index_type& arg_count);

  ThreadVectorRangeBoundariesStruct(member_type const& arg_thread,
                                    const index_type& arg_begin,
                                    const index_type& arg_end);

#endif
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {

// template<typename iType>
// KOKKOS_INLINE_FUNCTION
// Impl::TeamThreadRangeBoundariesStruct< iType, Impl::TaskExec< Kokkos::Cuda >
// > TeamThreadRange( const Impl::TaskExec< Kokkos::Cuda > & thread, const iType
// & count )
//{
//  return Impl::TeamThreadRangeBoundariesStruct< iType, Impl::TaskExec<
//  Kokkos::Cuda > >( thread, count );
//}
//
// template<typename iType1, typename iType2>
// KOKKOS_INLINE_FUNCTION
// Impl::TeamThreadRangeBoundariesStruct
//  < std::common_type_t<iType1,iType2>
//  , Impl::TaskExec< Kokkos::Cuda > >
// TeamThreadRange( const Impl::TaskExec< Kokkos::Cuda > & thread
//               , const iType1 & begin, const iType2 & end )
//{
//  using iType = std::common_type_t< iType1, iType2 >;
//  return Impl::TeamThreadRangeBoundariesStruct< iType, Impl::TaskExec<
//  Kokkos::Cuda > >(
//           thread, iType(begin), iType(end) );
//}
//
// template<typename iType>
// KOKKOS_INLINE_FUNCTION
// Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Cuda >
// > ThreadVectorRange( const Impl::TaskExec< Kokkos::Cuda > & thread
//                 , const iType & count )
//{
//  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec<
//  Kokkos::Cuda > >(thread,count);
//}
//
// template<typename iType>
// KOKKOS_INLINE_FUNCTION
// Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Cuda >
// > ThreadVectorRange( const Impl::TaskExec< Kokkos::Cuda > & thread
//                 , const iType & arg_begin
//                 , const iType & arg_end )
//{
//  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec<
//  Kokkos::Cuda > >(thread,arg_begin,arg_end);
//}

// KOKKOS_INLINE_FUNCTION
// Impl::ThreadSingleStruct<Impl::TaskExec< Kokkos::Cuda > >
// PerTeam(const Impl::TaskExec< Kokkos::Cuda >& thread)
// {
//   return Impl::ThreadSingleStruct<Impl::TaskExec< Kokkos::Cuda > >(thread);
// }

// KOKKOS_INLINE_FUNCTION
// Impl::VectorSingleStruct<Impl::TaskExec< Kokkos::Cuda > >
// PerThread(const Impl::TaskExec< Kokkos::Cuda >& thread)
// {
//   return Impl::VectorSingleStruct<Impl::TaskExec< Kokkos::Cuda > >(thread);
// }

/** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each
 * i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team.
 */
template <typename iType, class Lambda, class Scheduler>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Cuda, Scheduler>>& loop_boundaries,
    const Lambda& lambda) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i);
  }
}

template <typename iType, class Lambda, class Scheduler>
KOKKOS_INLINE_FUNCTION void parallel_for(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Cuda, Scheduler>>& loop_boundaries,
    const Lambda& lambda) {
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i);
  }
}

// reduce across corresponding lanes between team members within warp
// assume stride*team_size == warp_size
template <typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void strided_shfl_warp_reduction(const JoinType& join,
                                                        ValueType& val,
                                                        int team_size,
                                                        int stride) {
  for (int lane_delta = (team_size * stride) >> 1; lane_delta >= stride;
       lane_delta >>= 1) {
    join(val, Kokkos::shfl_down(val, lane_delta, team_size * stride));
  }
}

// multiple within-warp non-strided reductions
template <typename ValueType, class JoinType>
KOKKOS_INLINE_FUNCTION void multi_shfl_warp_reduction(const JoinType& join,
                                                      ValueType& val,
                                                      int vec_length) {
  for (int lane_delta = vec_length >> 1; lane_delta; lane_delta >>= 1) {
    join(val, Kokkos::shfl_down(val, lane_delta, vec_length));
  }
}

// broadcast within warp
template <class ValueType>
KOKKOS_INLINE_FUNCTION ValueType shfl_warp_broadcast(ValueType& val,
                                                     int src_lane, int width) {
  if (1 < width) {
    return Kokkos::shfl(val, src_lane, width);
  } else {
    return val;
  }
}

/*// all-reduce across corresponding vector lanes between team members within
warp
// assume vec_length*team_size == warp_size
// blockDim.x == vec_length == stride
// blockDim.y == team_size
// threadIdx.x == position in vec
// threadIdx.y == member number
template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce
  (const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::TaskExec<
Kokkos::Cuda > >& loop_boundaries, const Lambda & lambda, const JoinType& join,
   ValueType& initialized_result) {

  ValueType result = initialized_result;
  for( iType i = loop_boundaries.start; i < loop_boundaries.end;
i+=loop_boundaries.increment) { lambda(i,result);
  }
  initialized_result = result;

  strided_shfl_warp_reduction<ValueType, JoinType>(
                          join,
                          initialized_result,
                          loop_boundaries.thread.team_size(),
                          blockDim.x);
  initialized_result = shfl_warp_broadcast<ValueType>( initialized_result,
threadIdx.x, Impl::CudaTraits::WarpSize );
}*/

// all-reduce across corresponding vector lanes between team members within warp
// if no join() provided, use sum
// assume vec_length*team_size == warp_size
// blockDim.x == vec_length == stride
// blockDim.y == team_size
// threadIdx.x == position in vec
// threadIdx.y == member number
template <typename iType, class Lambda, typename ValueType, class Scheduler>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Cuda, Scheduler>>& loop_boundaries,
    const Lambda& lambda, ValueType& initialized_result) {
  // TODO @internal_documentation what is the point of creating this temporary?
  ValueType result = initialized_result;
  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }
  initialized_result = result;

  if (1 < loop_boundaries.thread.team_size()) {
    strided_shfl_warp_reduction(
        [&](ValueType& val1, const ValueType& val2) { val1 += val2; },
        initialized_result, loop_boundaries.thread.team_size(), blockDim.x);

    initialized_result = shfl_warp_broadcast<ValueType>(
        initialized_result, threadIdx.x, Impl::CudaTraits::WarpSize);
  }
}

template <typename iType, class Lambda, typename ReducerType, class Scheduler>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Cuda, Scheduler>>& loop_boundaries,
    const Lambda& lambda, const ReducerType& reducer) {
  using ValueType = typename ReducerType::value_type;
  // TODO @internal_documentation what is the point of creating this temporary?
  ValueType result = ValueType();
  reducer.init(result);

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }

  if (1 < loop_boundaries.thread.team_size()) {
    strided_shfl_warp_reduction(
        [&](ValueType& val1, const ValueType& val2) {
          reducer.join(val1, val2);
        },
        result, loop_boundaries.thread.team_size(), blockDim.x);

    reducer.reference() = shfl_warp_broadcast<ValueType>(
        result, threadIdx.x, Impl::CudaTraits::WarpSize);
  } else {
    reducer.reference() = result;
  }
}
// all-reduce within team members within warp
// assume vec_length*team_size == warp_size
// blockDim.x == vec_length == stride
// blockDim.y == team_size
// threadIdx.x == position in vec
// threadIdx.y == member number
/*template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce
  (const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec<
Kokkos::Cuda > >& loop_boundaries, const Lambda & lambda, const JoinType& join,
   ValueType& initialized_result) {

  ValueType result = initialized_result;
  for( iType i = loop_boundaries.start; i < loop_boundaries.end;
i+=loop_boundaries.increment) { lambda(i,result);
  }
  initialized_result = result;

  multi_shfl_warp_reduction<ValueType, JoinType>(join, initialized_result,
blockDim.x); initialized_result = shfl_warp_broadcast<ValueType>(
initialized_result, 0, blockDim.x );
}*/

// all-reduce within team members within warp
// if no join() provided, use sum
// assume vec_length*team_size == warp_size
// blockDim.x == vec_length == stride
// blockDim.y == team_size
// threadIdx.x == position in vec
// threadIdx.y == member number
template <typename iType, class Lambda, typename ValueType, class Scheduler>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Cuda, Scheduler>>& loop_boundaries,
    const Lambda& lambda, ValueType& initialized_result) {
  ValueType result = initialized_result;

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }

  initialized_result = result;

  if (1 < loop_boundaries.thread.team_size()) {
    // initialized_result = multi_shfl_warp_reduction(
    multi_shfl_warp_reduction(
        [&](ValueType& val1, const ValueType& val2) { val1 += val2; },
        initialized_result, blockDim.x);

    initialized_result =
        shfl_warp_broadcast<ValueType>(initialized_result, 0, blockDim.x);
  }
}

template <typename iType, class Lambda, typename ReducerType, class Scheduler>
KOKKOS_INLINE_FUNCTION void parallel_reduce(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Cuda, Scheduler>>& loop_boundaries,
    const Lambda& lambda, const ReducerType& reducer) {
  using ValueType = typename ReducerType::value_type;

  ValueType result = ValueType();
  reducer.init(result);

  for (iType i = loop_boundaries.start; i < loop_boundaries.end;
       i += loop_boundaries.increment) {
    lambda(i, result);
  }

  if (1 < loop_boundaries.thread.team_size()) {
    multi_shfl_warp_reduction(
        [&](ValueType& val1, const ValueType& val2) {
          reducer.join(val1, val2);
        },
        result, blockDim.x);

    reducer.reference() = shfl_warp_broadcast<ValueType>(result, 0, blockDim.x);
  } else {
    reducer.reference() = result;
  }
}
// scan across corresponding vector lanes between team members within warp
// assume vec_length*team_size == warp_size
// blockDim.x == vec_length == stride
// blockDim.y == team_size
// threadIdx.x == position in vec
// threadIdx.y == member number
template <typename iType, class Closure, class Scheduler>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::TeamThreadRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Cuda, Scheduler>>& loop_boundaries,
    const Closure& closure) {
  // Extract value_type from closure

  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, Closure>::value_type;

  if (1 < loop_boundaries.thread.team_size()) {
    // make sure all threads perform all loop iterations
    const iType bound = loop_boundaries.end + loop_boundaries.start;
    const int lane    = threadIdx.y * blockDim.x;

    value_type accum = 0;
    value_type val, y, local_total;

    for (iType i = loop_boundaries.start; i < bound;
         i += loop_boundaries.increment) {
      val = 0;
      if (i < loop_boundaries.end) closure(i, val, false);

      // intra-blockDim.y exclusive scan on 'val'
      // accum = accumulated, sum in total for this iteration

      // INCLUSIVE scan
      for (int offset = blockDim.x; offset < Impl::CudaTraits::WarpSize;
           offset <<= 1) {
        y = Kokkos::shfl_up(val, offset, Impl::CudaTraits::WarpSize);
        if (lane >= offset) {
          val += y;
        }
      }

      // pass accum to all threads
      local_total = shfl_warp_broadcast<value_type>(
          val, threadIdx.x + Impl::CudaTraits::WarpSize - blockDim.x,
          Impl::CudaTraits::WarpSize);

      // make EXCLUSIVE scan by shifting values over one
      val = Kokkos::shfl_up(val, blockDim.x, Impl::CudaTraits::WarpSize);
      if (threadIdx.y == 0) {
        val = 0;
      }

      val += accum;
      if (i < loop_boundaries.end) closure(i, val, true);
      accum += local_total;
    }
  } else {
    value_type accum = 0;
    for (iType i = loop_boundaries.start; i < loop_boundaries.end;
         i += loop_boundaries.increment) {
      closure(i, accum, true);
    }
  }
}

// scan within team member (vector) within warp
// assume vec_length*team_size == warp_size
// blockDim.x == vec_length == stride
// blockDim.y == team_size
// threadIdx.x == position in vec
// threadIdx.y == member number
template <typename iType, class Closure, class Scheduler>
KOKKOS_INLINE_FUNCTION void parallel_scan(
    const Impl::ThreadVectorRangeBoundariesStruct<
        iType, Impl::TaskExec<Kokkos::Cuda, Scheduler>>& loop_boundaries,
    const Closure& closure) {
  // Extract value_type from closure

  using value_type = typename Kokkos::Impl::FunctorAnalysis<
      Kokkos::Impl::FunctorPatternInterface::SCAN, void, Closure>::value_type;

  if (1 < loop_boundaries.thread.team_size()) {
    // make sure all threads perform all loop iterations
    const iType bound = loop_boundaries.end + loop_boundaries.start;

    value_type accum = 0;
    value_type val, y, local_total;

    for (iType i = loop_boundaries.start; i < bound;
         i += loop_boundaries.increment) {
      val = 0;
      if (i < loop_boundaries.end) closure(i, val, false);

      // intra-blockDim.x exclusive scan on 'val'
      // accum = accumulated, sum in total for this iteration

      // INCLUSIVE scan
      for (int offset = 1; offset < blockDim.x; offset <<= 1) {
        y = Kokkos::shfl_up(val, offset, blockDim.x);
        if (threadIdx.x >= offset) {
          val += y;
        }
      }

      // pass accum to all threads
      local_total =
          shfl_warp_broadcast<value_type>(val, blockDim.x - 1, blockDim.x);

      // make EXCLUSIVE scan by shifting values over one
      val = Kokkos::shfl_up(val, 1, blockDim.x);
      if (threadIdx.x == 0) {
        val = 0;
      }

      val += accum;
      if (i < loop_boundaries.end) closure(i, val, true);
      accum += local_total;
    }
  } else {
    value_type accum = 0;
    for (iType i = loop_boundaries.start; i < loop_boundaries.end;
         i += loop_boundaries.increment) {
      closure(i, accum, true);
    }
  }
}

} /* namespace Kokkos */

namespace Kokkos {

template <class FunctorType, class Scheduler>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::TaskExec<Kokkos::Cuda, Scheduler>>&,
    const FunctorType& lambda) {
#ifdef __CUDA_ARCH__
  if (threadIdx.x == 0) lambda();
#endif
}

template <class FunctorType, class Scheduler>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::TaskExec<Kokkos::Cuda, Scheduler>>&,
    const FunctorType& lambda) {
#ifdef __CUDA_ARCH__
  if (threadIdx.x == 0 && threadIdx.y == 0) lambda();
#endif
}

template <class FunctorType, class ValueType, class Scheduler>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::VectorSingleStruct<Impl::TaskExec<Kokkos::Cuda, Scheduler>>& s,
    const FunctorType& lambda, ValueType& val) {
#ifdef __CUDA_ARCH__
  if (threadIdx.x == 0) lambda(val);
  if (1 < s.team_member.team_size()) {
    val = shfl(val, 0, blockDim.x);
  }
#endif
}

template <class FunctorType, class ValueType, class Scheduler>
KOKKOS_INLINE_FUNCTION void single(
    const Impl::ThreadSingleStruct<Impl::TaskExec<Kokkos::Cuda, Scheduler>>&
        single_struct,
    const FunctorType& lambda, ValueType& val) {
#ifdef __CUDA_ARCH__
  if (threadIdx.x == 0 && threadIdx.y == 0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val, 0);
#endif
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#undef KOKKOS_IMPL_CUDA_SYNCWARP_OR_RETURN

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_CUDA_TASK_HPP */
