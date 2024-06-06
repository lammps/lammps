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
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include "Threads/Kokkos_Threads_Instance.hpp"
#endif

#include <Kokkos_Macros.hpp>

#include <utility>
#include <iostream>
#include <sstream>
#include <thread>

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_CPUDiscovery.hpp>
#include <impl/Kokkos_Tools.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace {

// std::thread compatible driver.
// Recovery from an exception would require constant intra-thread health
// verification; which would negatively impact runtime.  As such simply
// abort the process.
void internal_cppthread_driver() {
  try {
    ThreadsInternal::driver();
  } catch (const std::exception &x) {
    std::cerr << "Exception thrown from worker thread: " << x.what()
              << std::endl;
    std::cerr.flush();
    std::abort();
  } catch (...) {
    std::cerr << "Exception thrown from worker thread" << std::endl;
    std::cerr.flush();
    std::abort();
  }
}

ThreadsInternal s_threads_process;
ThreadsInternal *s_threads_exec[ThreadsInternal::MAX_THREAD_COUNT] = {nullptr};
std::thread::id s_threads_pid[ThreadsInternal::MAX_THREAD_COUNT];
std::pair<unsigned, unsigned>
    s_threads_coord[ThreadsInternal::MAX_THREAD_COUNT];

int s_thread_pool_size[3] = {0, 0, 0};

void (*volatile s_current_function)(ThreadsInternal &, const void *);
const void *volatile s_current_function_arg = nullptr;

inline unsigned fan_size(const unsigned rank, const unsigned size) {
  const unsigned rank_rev = size - (rank + 1);
  unsigned count          = 0;
  for (unsigned n = 1; (rank_rev + n < size) && !(rank_rev & n); n <<= 1) {
    ++count;
  }
  return count;
}

void wait_yield(volatile ThreadState &flag, const ThreadState value) {
  while (value == flag) {
    std::this_thread::yield();
  }
}

}  // namespace
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

bool ThreadsInternal::is_process() {
  static const std::thread::id master_pid = std::this_thread::get_id();

  return master_pid == std::this_thread::get_id();
}

//----------------------------------------------------------------------------

void execute_function_noop(ThreadsInternal &, const void *) {}

void ThreadsInternal::driver() {
  SharedAllocationRecord<void, void>::tracking_enable();

  ThreadsInternal this_thread;

  while (this_thread.m_pool_state == ThreadState::Active) {
    (*s_current_function)(this_thread, s_current_function_arg);

    // Deactivate thread and wait for reactivation
    this_thread.m_pool_state = ThreadState::Inactive;

    wait_yield(this_thread.m_pool_state, ThreadState::Inactive);
  }
}

ThreadsInternal::ThreadsInternal()
    : m_pool_base(nullptr),
      m_scratch(nullptr),
      m_scratch_reduce_end(0),
      m_scratch_thread_end(0),
      m_pool_rank(0),
      m_pool_size(0),
      m_pool_fan_size(0),
      m_pool_state(ThreadState::Terminating) {
  if (&s_threads_process != this) {
    // The code in the if is executed by a spawned thread not by the root
    // thread
    ThreadsInternal *const nil = nullptr;

    // Which entry in 's_threads_exec', possibly determined from hwloc binding
    const int entry = reinterpret_cast<size_t>(s_current_function_arg) <
                              size_t(s_thread_pool_size[0])
                          ? reinterpret_cast<size_t>(s_current_function_arg)
                          : size_t(Kokkos::hwloc::bind_this_thread(
                                s_thread_pool_size[0], s_threads_coord));

    // Given a good entry set this thread in the 's_threads_exec' array
    if (entry < s_thread_pool_size[0] &&
        nil == atomic_compare_exchange(s_threads_exec + entry, nil, this)) {
      m_pool_base     = s_threads_exec;
      m_pool_rank     = s_thread_pool_size[0] - (entry + 1);
      m_pool_rank_rev = s_thread_pool_size[0] - (pool_rank() + 1);
      m_pool_size     = s_thread_pool_size[0];
      m_pool_fan_size = fan_size(m_pool_rank, m_pool_size);
      m_pool_state    = ThreadState::Active;

      s_threads_pid[m_pool_rank] = std::this_thread::get_id();

      // Inform spawning process that the threads_exec entry has been set.
      s_threads_process.m_pool_state = ThreadState::Active;
    } else {
      // Inform spawning process that the threads_exec entry could not be set.
      s_threads_process.m_pool_state = ThreadState::Terminating;
    }
  } else {
    // Enables 'parallel_for' to execute on unitialized Threads device
    m_pool_rank  = 0;
    m_pool_size  = 1;
    m_pool_state = ThreadState::Inactive;

    s_threads_pid[m_pool_rank] = std::this_thread::get_id();
  }
}

ThreadsInternal::~ThreadsInternal() {
  const unsigned entry = m_pool_size - (m_pool_rank + 1);

  if (m_scratch) {
    Kokkos::kokkos_free<Kokkos::HostSpace>(m_scratch);
    m_scratch = nullptr;
  }

  m_pool_base          = nullptr;
  m_scratch_reduce_end = 0;
  m_scratch_thread_end = 0;
  m_pool_rank          = 0;
  m_pool_size          = 0;
  m_pool_fan_size      = 0;

  m_pool_state = ThreadState::Terminating;

  if (&s_threads_process != this && entry < MAX_THREAD_COUNT) {
    ThreadsInternal *const nil = nullptr;

    atomic_compare_exchange(s_threads_exec + entry, this, nil);

    s_threads_process.m_pool_state = ThreadState::Terminating;
  }
}

ThreadsInternal *ThreadsInternal::get_thread(const int init_thread_rank) {
  ThreadsInternal *const th =
      init_thread_rank < s_thread_pool_size[0]
          ? s_threads_exec[s_thread_pool_size[0] - (init_thread_rank + 1)]
          : nullptr;

  if (nullptr == th || th->m_pool_rank != init_thread_rank) {
    std::ostringstream msg;
    msg << "Kokkos::Impl::ThreadsInternal::get_thread ERROR : "
        << "thread " << init_thread_rank << " of " << s_thread_pool_size[0];
    if (nullptr == th) {
      msg << " does not exist";
    } else {
      msg << " has wrong thread_rank " << th->m_pool_rank;
    }
    Kokkos::Impl::throw_runtime_exception(msg.str());
  }

  return th;
}

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void ThreadsInternal::verify_is_process(const std::string &name,
                                        const bool initialized) {
  if (!is_process()) {
    std::string msg(name);
    msg.append(
        " FAILED : Called by a worker thread, can only be called by the master "
        "process.");
    Kokkos::Impl::throw_runtime_exception(msg);
  }

  if (initialized && 0 == s_thread_pool_size[0]) {
    std::string msg(name);
    msg.append(" FAILED : Threads not initialized.");
    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
KOKKOS_DEPRECATED int ThreadsInternal::in_parallel() {
  // A thread function is in execution and
  // the function argument is not the special threads process argument and
  // the master process is a worker or is not the master process.
  return s_current_function && (&s_threads_process != s_current_function_arg) &&
         (s_threads_process.m_pool_base || !is_process());
}
#endif
void ThreadsInternal::fence() {
  fence("Kokkos::ThreadsInternal::fence: Unnamed Instance Fence");
}
void ThreadsInternal::fence(const std::string &name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::Threads>(
      name, Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{1},
      internal_fence);
}

// Wait for root thread to become inactive
void ThreadsInternal::internal_fence() {
  if (s_thread_pool_size[0]) {
    // Wait for the root thread to complete:
    Impl::spinwait_while_equal(s_threads_exec[0]->m_pool_state,
                               ThreadState::Active);
  }

  s_current_function     = nullptr;
  s_current_function_arg = nullptr;

  // Make sure function and arguments are cleared before
  // potentially re-activating threads with a subsequent launch.
  memory_fence();
}

/** \brief  Begin execution of the asynchronous functor */
void ThreadsInternal::start(void (*func)(ThreadsInternal &, const void *),
                            const void *arg) {
  verify_is_process("ThreadsInternal::start", true);

  if (s_current_function || s_current_function_arg) {
    Kokkos::Impl::throw_runtime_exception(
        std::string("ThreadsInternal::start() FAILED : already executing"));
  }

  s_current_function     = func;
  s_current_function_arg = arg;

  // Make sure function and arguments are written before activating threads.
  memory_fence();

  // Activate threads. The spawned threads will start working on
  // s_current_function. The root thread is only set to active, we still need to
  // call s_current_function.
  for (int i = s_thread_pool_size[0]; 0 < i--;) {
    s_threads_exec[i]->m_pool_state = ThreadState::Active;
  }

  if (s_threads_process.m_pool_size) {
    // Master process is the root thread, run it:
    (*func)(s_threads_process, arg);
    s_threads_process.m_pool_state = ThreadState::Inactive;
  }
}

//----------------------------------------------------------------------------

void ThreadsInternal::execute_resize_scratch_in_serial() {
  const unsigned begin = s_threads_process.m_pool_base ? 1 : 0;

  auto deallocate_scratch_memory = [](ThreadsInternal &exec) {
    if (exec.m_scratch) {
      Kokkos::kokkos_free<Kokkos::HostSpace>(exec.m_scratch);
      exec.m_scratch = nullptr;
    }
  };
  if (s_threads_process.m_pool_base) {
    for (unsigned i = s_thread_pool_size[0]; begin < i;) {
      deallocate_scratch_memory(*s_threads_exec[--i]);
    }
  }

  s_current_function     = &first_touch_allocate_thread_private_scratch;
  s_current_function_arg = &s_threads_process;

  // Make sure function and arguments are written before activating threads.
  memory_fence();

  for (unsigned i = s_thread_pool_size[0]; begin < i;) {
    ThreadsInternal &th = *s_threads_exec[--i];

    th.m_pool_state = ThreadState::Active;

    wait_yield(th.m_pool_state, ThreadState::Active);
  }

  if (s_threads_process.m_pool_base) {
    deallocate_scratch_memory(s_threads_process);
    s_threads_process.m_pool_state = ThreadState::Active;
    first_touch_allocate_thread_private_scratch(s_threads_process, nullptr);
    s_threads_process.m_pool_state = ThreadState::Inactive;
  }

  s_current_function_arg = nullptr;
  s_current_function     = nullptr;

  // Make sure function and arguments are cleared before proceeding.
  memory_fence();
}

//----------------------------------------------------------------------------

void *ThreadsInternal::root_reduce_scratch() {
  return s_threads_process.reduce_memory();
}

void ThreadsInternal::first_touch_allocate_thread_private_scratch(
    ThreadsInternal &exec, const void *) {
  exec.m_scratch_reduce_end = s_threads_process.m_scratch_reduce_end;
  exec.m_scratch_thread_end = s_threads_process.m_scratch_thread_end;

  if (s_threads_process.m_scratch_thread_end) {
    // Allocate tracked memory:
    {
      exec.m_scratch = Kokkos::kokkos_malloc<Kokkos::HostSpace>(
          "Kokkos::thread_scratch", s_threads_process.m_scratch_thread_end);
    }

    unsigned *ptr = reinterpret_cast<unsigned *>(exec.m_scratch);

    unsigned *const end =
        ptr + s_threads_process.m_scratch_thread_end / sizeof(unsigned);

    // touch on this thread
    while (ptr < end) *ptr++ = 0;
  }
}

void *ThreadsInternal::resize_scratch(size_t reduce_size, size_t thread_size) {
  enum { ALIGN_MASK = Kokkos::Impl::MEMORY_ALIGNMENT - 1 };

  fence();

  const size_t old_reduce_size = s_threads_process.m_scratch_reduce_end;
  const size_t old_thread_size = s_threads_process.m_scratch_thread_end -
                                 s_threads_process.m_scratch_reduce_end;

  reduce_size = (reduce_size + ALIGN_MASK) & ~ALIGN_MASK;
  thread_size = (thread_size + ALIGN_MASK) & ~ALIGN_MASK;

  // Increase size or deallocate completely.

  if ((old_reduce_size < reduce_size) || (old_thread_size < thread_size) ||
      ((reduce_size == 0 && thread_size == 0) &&
       (old_reduce_size != 0 || old_thread_size != 0))) {
    verify_is_process("ThreadsInternal::resize_scratch", true);

    s_threads_process.m_scratch_reduce_end = reduce_size;
    s_threads_process.m_scratch_thread_end = reduce_size + thread_size;

    execute_resize_scratch_in_serial();

    s_threads_process.m_scratch = s_threads_exec[0]->m_scratch;
  }

  return s_threads_process.m_scratch;
}

//----------------------------------------------------------------------------

void ThreadsInternal::print_configuration(std::ostream &s, const bool detail) {
  verify_is_process("ThreadsInternal::print_configuration", false);

  fence();

  s << "Kokkos::Threads";

#if defined(KOKKOS_ENABLE_THREADS)
  s << " KOKKOS_ENABLE_THREADS";
#endif
#if defined(KOKKOS_ENABLE_HWLOC)
  const unsigned numa_count     = Kokkos::hwloc::get_available_numa_count();
  const unsigned cores_per_numa = Kokkos::hwloc::get_available_cores_per_numa();
  const unsigned threads_per_core =
      Kokkos::hwloc::get_available_threads_per_core();

  s << " hwloc[" << numa_count << "x" << cores_per_numa << "x"
    << threads_per_core << "]";
#endif

  if (s_thread_pool_size[0]) {
    s << " threads[" << s_thread_pool_size[0] << "]"
      << " threads_per_numa[" << s_thread_pool_size[1] << "]"
      << " threads_per_core[" << s_thread_pool_size[2] << "]";
    if (nullptr == s_threads_process.m_pool_base) {
      s << " Asynchronous";
    }
    s << std::endl;

    if (detail) {
      for (int i = 0; i < s_thread_pool_size[0]; ++i) {
        ThreadsInternal *const th = s_threads_exec[i];

        if (th) {
          const int rank_rev = th->m_pool_size - (th->m_pool_rank + 1);

          s << " Thread[ " << th->m_pool_rank << " ]";

          s << " Fan{";
          for (int j = 0; j < th->m_pool_fan_size; ++j) {
            ThreadsInternal *const thfan = th->m_pool_base[rank_rev + (1 << j)];
            s << " [ " << thfan->m_pool_rank << " ]";
          }
          s << " }";

          if (th == &s_threads_process) {
            s << " is_process";
          }
        }
        s << std::endl;
      }
    }
  } else {
    s << " not initialized" << std::endl;
  }
}

//----------------------------------------------------------------------------

int ThreadsInternal::is_initialized() { return nullptr != s_threads_exec[0]; }

void ThreadsInternal::initialize(int thread_count_arg) {
  unsigned thread_count = thread_count_arg == -1 ? 0 : thread_count_arg;

  const bool is_initialized = 0 != s_thread_pool_size[0];

  unsigned thread_spawn_failed = 0;

  for (int i = 0; i < ThreadsInternal::MAX_THREAD_COUNT; i++)
    s_threads_exec[i] = nullptr;

  if (!is_initialized) {
    // If thread_count is zero then it will be given default values based upon
    // hwloc detection.
    const bool hwloc_avail = Kokkos::hwloc::available();
    const bool hwloc_can_bind =
        hwloc_avail && Kokkos::hwloc::can_bind_threads();

    if (thread_count == 0) {
      thread_count = hwloc_avail
                         ? Kokkos::hwloc::get_available_numa_count() *
                               Kokkos::hwloc::get_available_cores_per_numa() *
                               Kokkos::hwloc::get_available_threads_per_core()
                         : 1;
    }

    const bool allow_asynchronous_threadpool = false;
    unsigned use_numa_count                  = 0;
    unsigned use_cores_per_numa              = 0;
    hwloc::thread_mapping("Kokkos::Threads::initialize",
                          allow_asynchronous_threadpool, thread_count,
                          use_numa_count, use_cores_per_numa, s_threads_coord);

    const std::pair<unsigned, unsigned> proc_coord = s_threads_coord[0];

    // Synchronous with s_threads_coord[0] as the process core
    // Claim entry #0 for binding the process core.
    s_threads_coord[0] = std::pair<unsigned, unsigned>(~0u, ~0u);

    s_thread_pool_size[0] = thread_count;
    s_thread_pool_size[1] = s_thread_pool_size[0] / use_numa_count;
    s_thread_pool_size[2] = s_thread_pool_size[1] / use_cores_per_numa;
    s_current_function =
        &execute_function_noop;  // Initialization work function

    for (unsigned ith = 1; ith < thread_count; ++ith) {
      s_threads_process.m_pool_state = ThreadState::Inactive;

      // If hwloc available then spawned thread will
      // choose its own entry in 's_threads_coord'
      // otherwise specify the entry.
      s_current_function_arg =
          reinterpret_cast<void *>(hwloc_can_bind ? ~0u : ith);

      // Make sure all outstanding memory writes are complete
      // before spawning the new thread.
      memory_fence();

      // Spawn thread executing the 'driver()' function.
      // Wait until spawned thread has attempted to initialize.
      // If spawning and initialization is successful then
      // an entry in 's_threads_exec' will be assigned.
      std::thread t(internal_cppthread_driver);
      t.detach();
      wait_yield(s_threads_process.m_pool_state, ThreadState::Inactive);
      if (s_threads_process.m_pool_state == ThreadState::Terminating) break;
    }

    // Wait for all spawned threads to deactivate before zeroing the function.

    for (unsigned ith = 1; ith < thread_count; ++ith) {
      // Try to protect against cache coherency failure by casting to volatile.
      ThreadsInternal *const th =
          ((ThreadsInternal * volatile *)s_threads_exec)[ith];
      if (th) {
        wait_yield(th->m_pool_state, ThreadState::Active);
      } else {
        ++thread_spawn_failed;
      }
    }

    s_current_function             = nullptr;
    s_current_function_arg         = nullptr;
    s_threads_process.m_pool_state = ThreadState::Inactive;

    memory_fence();

    if (!thread_spawn_failed) {
      // Bind process to the core on which it was located before spawning
      // occurred
      if (hwloc_can_bind) {
        Kokkos::hwloc::bind_this_thread(proc_coord);
      }

      s_threads_exec[0]             = &s_threads_process;
      s_threads_process.m_pool_base = s_threads_exec;
      s_threads_process.m_pool_rank =
          thread_count - 1;  // Reversed for scan-compatible reductions
      s_threads_process.m_pool_size     = thread_count;
      s_threads_process.m_pool_fan_size = fan_size(
          s_threads_process.m_pool_rank, s_threads_process.m_pool_size);
      s_threads_pid[s_threads_process.m_pool_rank] = std::this_thread::get_id();

      // Initial allocations:
      ThreadsInternal::resize_scratch(1024, 1024);
    } else {
      s_thread_pool_size[0] = 0;
      s_thread_pool_size[1] = 0;
      s_thread_pool_size[2] = 0;
    }
  }

  if (is_initialized || thread_spawn_failed) {
    std::ostringstream msg;

    msg << "Kokkos::Threads::initialize ERROR";

    if (is_initialized) {
      msg << " : already initialized";
    }
    if (thread_spawn_failed) {
      msg << " : failed to spawn " << thread_spawn_failed << " threads";
    }

    Kokkos::Impl::throw_runtime_exception(msg.str());
  }

  // Check for over-subscription
  auto const reported_ranks = mpi_ranks_per_node();
  auto const mpi_local_size = reported_ranks < 0 ? 1 : reported_ranks;
  int const procs_per_node  = std::thread::hardware_concurrency();
  if (Kokkos::show_warnings() &&
      (mpi_local_size * long(thread_count) > procs_per_node)) {
    std::cerr << "Kokkos::Threads::initialize WARNING: You are likely "
                 "oversubscribing your CPU cores."
              << std::endl;
    std::cerr << "                                    Detected: "
              << procs_per_node << " cores per node." << std::endl;
    std::cerr << "                                    Detected: "
              << mpi_local_size << " MPI_ranks per node." << std::endl;
    std::cerr << "                                    Requested: "
              << thread_count << " threads per process." << std::endl;
  }

  Impl::SharedAllocationRecord<void, void>::tracking_enable();
}

//----------------------------------------------------------------------------

void ThreadsInternal::finalize() {
  verify_is_process("ThreadsInternal::finalize", false);

  fence();

  resize_scratch(0, 0);

  const unsigned begin = s_threads_process.m_pool_base ? 1 : 0;

  for (unsigned i = s_thread_pool_size[0]; begin < i--;) {
    if (s_threads_exec[i]) {
      s_threads_exec[i]->m_pool_state = ThreadState::Terminating;

      wait_yield(s_threads_process.m_pool_state, ThreadState::Inactive);

      s_threads_process.m_pool_state = ThreadState::Inactive;
    }

    s_threads_pid[i] = std::thread::id();
  }

  if (s_threads_process.m_pool_base) {
    (&s_threads_process)->~ThreadsInternal();
    s_threads_exec[0] = nullptr;
  }

  if (Kokkos::hwloc::can_bind_threads()) {
    Kokkos::hwloc::unbind_this_thread();
  }

  s_thread_pool_size[0] = 0;
  s_thread_pool_size[1] = 0;
  s_thread_pool_size[2] = 0;

  // Reset master thread to run solo.
  s_threads_process.m_pool_base     = nullptr;
  s_threads_process.m_pool_rank     = 0;
  s_threads_process.m_pool_size     = 1;
  s_threads_process.m_pool_fan_size = 0;
  s_threads_process.m_pool_state    = ThreadState::Inactive;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
int Threads::concurrency() { return impl_thread_pool_size(0); }
#else
int Threads::concurrency() const { return impl_thread_pool_size(0); }
#endif

void Threads::fence(const std::string &name) const {
  Impl::ThreadsInternal::fence(name);
}

Threads &Threads::impl_instance(int) {
  static Threads t;
  return t;
}

int Threads::impl_thread_pool_rank_host() {
  const std::thread::id pid = std::this_thread::get_id();
  int i                     = 0;
  while ((i < Impl::s_thread_pool_size[0]) && (pid != Impl::s_threads_pid[i])) {
    ++i;
  }
  return i;
}

int Threads::impl_thread_pool_size(int depth) {
  return Impl::s_thread_pool_size[depth];
}

const char *Threads::name() { return "Threads"; }

namespace Impl {

int g_threads_space_factory_initialized =
    initialize_space_factory<Threads>("050_Threads");

}  // namespace Impl

} /* namespace Kokkos */
