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

#ifndef KOKKOS_THREADSEXEC_HPP
#define KOKKOS_THREADSEXEC_HPP

#include <Kokkos_Macros.hpp>

#include <cstdio>
#include <ostream>
#include <utility>

#include <impl/Kokkos_Spinwait.hpp>

#include <Kokkos_Atomic.hpp>
#include <Kokkos_Pair.hpp>

#include <impl/Kokkos_ConcurrentBitset.hpp>
#include <Threads/Kokkos_Threads.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
class ThreadsExec {
 public:
  // Fan array has log_2(NT) reduction threads plus 2 scan threads
  // Currently limited to 16k threads.
  enum { MAX_FAN_COUNT = 16 };
  enum { MAX_THREAD_COUNT = 1 << (MAX_FAN_COUNT - 2) };
  enum { VECTOR_LENGTH = 8 };

  /** \brief States of a worker thread */
  enum {
    Terminating  ///<  Termination in progress
    ,
    Inactive  ///<  Exists, waiting for work
    ,
    Active  ///<  Exists, performing work
    ,
    Rendezvous  ///<  Exists, waiting in a barrier or reduce

    ,
    ScanCompleted,
    ScanAvailable,
    ReductionAvailable
  };

 private:
  friend class Kokkos::Threads;

  // Fan-in operations' root is the highest ranking thread
  // to place the 'scan' reduction intermediate values on
  // the threads that need them.
  // For a simple reduction the thread location is arbitrary.

  ThreadsExec *const *m_pool_base;  ///< Base for pool fan-in

  void *m_scratch;
  int m_scratch_reduce_end;
  size_t m_scratch_thread_end;
  int m_numa_rank;
  int m_numa_core_rank;
  int m_pool_rank;
  int m_pool_rank_rev;
  int m_pool_size;
  int m_pool_fan_size;
  int volatile m_pool_state;  ///< State for global synchronizations

  // Members for dynamic scheduling
  // Which thread am I stealing from currently
  int m_current_steal_target;
  // This thread's owned work_range
  alignas(16) Kokkos::pair<long, long> m_work_range;
  // Team Offset if one thread determines work_range for others
  long m_team_work_index;

  // Is this thread stealing (i.e. its owned work_range is exhausted
  bool m_stealing;

  static void global_lock();
  static void global_unlock();
  static void spawn();

  static void first_touch_allocate_thread_private_scratch(ThreadsExec &,
                                                          const void *);
  static void execute_sleep(ThreadsExec &, const void *);

  ThreadsExec(const ThreadsExec &);
  ThreadsExec &operator=(const ThreadsExec &);

  static void execute_resize_scratch_in_serial();

 public:
  KOKKOS_INLINE_FUNCTION int pool_size() const { return m_pool_size; }
  KOKKOS_INLINE_FUNCTION int pool_rank() const { return m_pool_rank; }
  KOKKOS_INLINE_FUNCTION int numa_rank() const { return m_numa_rank; }
  KOKKOS_INLINE_FUNCTION int numa_core_rank() const { return m_numa_core_rank; }
  inline long team_work_index() const { return m_team_work_index; }

  static int get_thread_count();
  static ThreadsExec *get_thread(const int init_thread_rank);

  inline void *reduce_memory() const { return m_scratch; }
  KOKKOS_INLINE_FUNCTION void *scratch_memory() const {
    return reinterpret_cast<unsigned char *>(m_scratch) + m_scratch_reduce_end;
  }

  KOKKOS_INLINE_FUNCTION int volatile &state() { return m_pool_state; }
  KOKKOS_INLINE_FUNCTION ThreadsExec *const *pool_base() const {
    return m_pool_base;
  }

  static void driver(void);

  ~ThreadsExec();
  ThreadsExec();

  static void *resize_scratch(size_t reduce_size, size_t thread_size);

  static void *root_reduce_scratch();

  static bool is_process();

  static void verify_is_process(const std::string &, const bool initialized);

  static int is_initialized();

  static void initialize(int thread_count);

  static void finalize();

  /* Given a requested team size, return valid team size */
  static unsigned team_size_valid(unsigned);

  static void print_configuration(std::ostream &, const bool detail = false);

  //------------------------------------

  static void wait_yield(volatile int &, const int);

  //------------------------------------
  // All-thread functions:

  inline int all_reduce(const int value) {
    // Make sure there is enough scratch space:
    const int rev_rank = m_pool_size - (m_pool_rank + 1);

    *static_cast<volatile int *>(reduce_memory()) = value;

    memory_fence();

    // Fan-in reduction with highest ranking thread as the root
    for (int i = 0; i < m_pool_fan_size; ++i) {
      // Wait: Active -> Rendezvous
      Impl::spinwait_while_equal<int>(
          m_pool_base[rev_rank + (1 << i)]->m_pool_state, ThreadsExec::Active);
    }

    if (rev_rank) {
      m_pool_state = ThreadsExec::Rendezvous;
      // Wait: Rendezvous -> Active
      Impl::spinwait_while_equal<int>(m_pool_state, ThreadsExec::Rendezvous);
    } else {
      // Root thread does the reduction and broadcast

      int accum = 0;

      for (int rank = 0; rank < m_pool_size; ++rank) {
        accum +=
            *static_cast<volatile int *>(get_thread(rank)->reduce_memory());
      }

      for (int rank = 0; rank < m_pool_size; ++rank) {
        *static_cast<volatile int *>(get_thread(rank)->reduce_memory()) = accum;
      }

      memory_fence();

      for (int rank = 0; rank < m_pool_size; ++rank) {
        get_thread(rank)->m_pool_state = ThreadsExec::Active;
      }
    }

    return *static_cast<volatile int *>(reduce_memory());
  }

  inline void barrier() {
    // Make sure there is enough scratch space:
    const int rev_rank = m_pool_size - (m_pool_rank + 1);

    memory_fence();

    // Fan-in reduction with highest ranking thread as the root
    for (int i = 0; i < m_pool_fan_size; ++i) {
      // Wait: Active -> Rendezvous
      Impl::spinwait_while_equal<int>(
          m_pool_base[rev_rank + (1 << i)]->m_pool_state, ThreadsExec::Active);
    }

    if (rev_rank) {
      m_pool_state = ThreadsExec::Rendezvous;
      // Wait: Rendezvous -> Active
      Impl::spinwait_while_equal<int>(m_pool_state, ThreadsExec::Rendezvous);
    } else {
      // Root thread does the reduction and broadcast

      memory_fence();

      for (int rank = 0; rank < m_pool_size; ++rank) {
        get_thread(rank)->m_pool_state = ThreadsExec::Active;
      }
    }
  }

  //------------------------------------
  // All-thread functions:

  template <class FunctorType>
  inline void fan_in_reduce(const FunctorType &f) const {
    const int rev_rank = m_pool_size - (m_pool_rank + 1);

    for (int i = 0; i < m_pool_fan_size; ++i) {
      ThreadsExec &fan = *m_pool_base[rev_rank + (1 << i)];

      Impl::spinwait_while_equal<int>(fan.m_pool_state, ThreadsExec::Active);

      f.join(
          reinterpret_cast<typename FunctorType::value_type *>(reduce_memory()),
          reinterpret_cast<const typename FunctorType::value_type *>(
              fan.reduce_memory()));
    }

    if (!rev_rank) {
      f.final(reinterpret_cast<typename FunctorType::value_type *>(
          reduce_memory()));
    }

    //  This thread has updated 'reduce_memory()' and upon returning
    //  from this function will set 'm_pool_state' to inactive.
    //  If this is a non-root thread then setting 'm_pool_state'
    //  to inactive triggers another thread to exit a spinwait
    //  and read the 'reduce_memory'.
    //  Must 'memory_fence()' to guarantee that storing the update to
    //  'reduce_memory()' will complete before storing the the update to
    //  'm_pool_state'.

    memory_fence();
  }

  inline void fan_in() const {
    const int rev_rank = m_pool_size - (m_pool_rank + 1);

    for (int i = 0; i < m_pool_fan_size; ++i) {
      Impl::spinwait_while_equal<int>(
          m_pool_base[rev_rank + (1 << i)]->m_pool_state, ThreadsExec::Active);
    }
  }

  template <class FunctorType>
  inline void scan_large(const FunctorType &f) {
    // Sequence of states:
    //  0) Active             : entry and exit state
    //  1) ReductionAvailable : reduction value available
    //  2) ScanAvailable      : inclusive scan value available
    //  3) Rendezvous         : All threads inclusive scan value are available
    //  4) ScanCompleted      : exclusive scan value copied

    using scalar_type = typename FunctorType::value_type;

    const int rev_rank   = m_pool_size - (m_pool_rank + 1);
    const unsigned count = FunctorType::value_count(f);

    scalar_type *const work_value = (scalar_type *)reduce_memory();

    //--------------------------------
    // Fan-in reduction with highest ranking thread as the root
    for (int i = 0; i < m_pool_fan_size; ++i) {
      ThreadsExec &fan = *m_pool_base[rev_rank + (1 << i)];

      // Wait: Active -> ReductionAvailable (or ScanAvailable)
      Impl::spinwait_while_equal<int>(fan.m_pool_state, ThreadsExec::Active);
      f.join(work_value, fan.reduce_memory());
    }

    // Copy reduction value to scan value before releasing from this phase.
    for (unsigned i = 0; i < count; ++i) {
      work_value[i + count] = work_value[i];
    }

    if (rev_rank) {
      // Set: Active -> ReductionAvailable
      m_pool_state = ThreadsExec::ReductionAvailable;

      // Wait for contributing threads' scan value to be available.
      if ((1 << m_pool_fan_size) < (m_pool_rank + 1)) {
        ThreadsExec &th = *m_pool_base[rev_rank + (1 << m_pool_fan_size)];

        // Wait: Active             -> ReductionAvailable
        // Wait: ReductionAvailable -> ScanAvailable
        Impl::spinwait_while_equal<int>(th.m_pool_state, ThreadsExec::Active);
        Impl::spinwait_while_equal<int>(th.m_pool_state,
                                        ThreadsExec::ReductionAvailable);

        f.join(work_value + count, ((scalar_type *)th.reduce_memory()) + count);
      }

      // This thread has completed inclusive scan
      // Set: ReductionAvailable -> ScanAvailable
      m_pool_state = ThreadsExec::ScanAvailable;

      // Wait for all threads to complete inclusive scan
      // Wait: ScanAvailable -> Rendezvous
      Impl::spinwait_while_equal<int>(m_pool_state, ThreadsExec::ScanAvailable);
    }

    //--------------------------------

    for (int i = 0; i < m_pool_fan_size; ++i) {
      ThreadsExec &fan = *m_pool_base[rev_rank + (1 << i)];
      // Wait: ReductionAvailable -> ScanAvailable
      Impl::spinwait_while_equal<int>(fan.m_pool_state,
                                      ThreadsExec::ReductionAvailable);
      // Set: ScanAvailable -> Rendezvous
      fan.m_pool_state = ThreadsExec::Rendezvous;
    }

    // All threads have completed the inclusive scan.
    // All non-root threads are in the Rendezvous state.
    // Threads are free to overwrite their reduction value.
    //--------------------------------

    if ((rev_rank + 1) < m_pool_size) {
      // Exclusive scan: copy the previous thread's inclusive scan value

      ThreadsExec &th = *m_pool_base[rev_rank + 1];  // Not the root thread

      const scalar_type *const src_value =
          ((scalar_type *)th.reduce_memory()) + count;

      for (unsigned j = 0; j < count; ++j) {
        work_value[j] = src_value[j];
      }
    } else {
      f.init(work_value);
    }

    //--------------------------------
    // Wait for all threads to copy previous thread's inclusive scan value
    // Wait for all threads: Rendezvous -> ScanCompleted
    for (int i = 0; i < m_pool_fan_size; ++i) {
      Impl::spinwait_while_equal<int>(
          m_pool_base[rev_rank + (1 << i)]->m_pool_state,
          ThreadsExec::Rendezvous);
    }
    if (rev_rank) {
      // Set: ScanAvailable -> ScanCompleted
      m_pool_state = ThreadsExec::ScanCompleted;
      // Wait: ScanCompleted -> Active
      Impl::spinwait_while_equal<int>(m_pool_state, ThreadsExec::ScanCompleted);
    }
    // Set: ScanCompleted -> Active
    for (int i = 0; i < m_pool_fan_size; ++i) {
      m_pool_base[rev_rank + (1 << i)]->m_pool_state = ThreadsExec::Active;
    }
  }

  template <class FunctorType>
  inline void scan_small(const FunctorType &f) {
    using scalar_type = typename FunctorType::value_type;

    const int rev_rank   = m_pool_size - (m_pool_rank + 1);
    const unsigned count = f.length();

    scalar_type *const work_value = (scalar_type *)reduce_memory();

    //--------------------------------
    // Fan-in reduction with highest ranking thread as the root
    for (int i = 0; i < m_pool_fan_size; ++i) {
      // Wait: Active -> Rendezvous
      Impl::spinwait_while_equal<int>(
          m_pool_base[rev_rank + (1 << i)]->m_pool_state, ThreadsExec::Active);
    }

    for (unsigned i = 0; i < count; ++i) {
      work_value[i + count] = work_value[i];
    }

    if (rev_rank) {
      m_pool_state = ThreadsExec::Rendezvous;
      // Wait: Rendezvous -> Active
      Impl::spinwait_while_equal<int>(m_pool_state, ThreadsExec::Rendezvous);
    } else {
      // Root thread does the thread-scan before releasing threads

      scalar_type *ptr_prev = nullptr;

      for (int rank = 0; rank < m_pool_size; ++rank) {
        scalar_type *const ptr =
            (scalar_type *)get_thread(rank)->reduce_memory();
        if (rank) {
          for (unsigned i = 0; i < count; ++i) {
            ptr[i] = ptr_prev[i + count];
          }
          f.join(ptr + count, ptr);
        } else {
          f.init(ptr);
        }
        ptr_prev = ptr;
      }
    }

    for (int i = 0; i < m_pool_fan_size; ++i) {
      m_pool_base[rev_rank + (1 << i)]->m_pool_state = ThreadsExec::Active;
    }
  }

  //------------------------------------
  /** \brief  Wait for previous asynchronous functor to
   *          complete and release the Threads device.
   *          Acquire the Threads device and start this functor.
   */
  static void start(void (*)(ThreadsExec &, const void *), const void *);

  static int in_parallel();
  static void fence();
  static void fence(const std::string &);
  static void internal_fence(
      Impl::fence_is_static is_static = Impl::fence_is_static::yes);
  static void internal_fence(
      const std::string &,
      Impl::fence_is_static is_static = Impl::fence_is_static::yes);
  static bool sleep();
  static bool wake();

  /* Dynamic Scheduling related functionality */
  // Initialize the work range for this thread
  inline void set_work_range(const long &begin, const long &end,
                             const long &chunk_size) {
    m_work_range.first = (begin + chunk_size - 1) / chunk_size;
    m_work_range.second =
        end > 0 ? (end + chunk_size - 1) / chunk_size : m_work_range.first;
  }

  // Claim and index from this thread's range from the beginning
  inline long get_work_index_begin() {
    Kokkos::pair<long, long> work_range_new = m_work_range;
    Kokkos::pair<long, long> work_range_old = work_range_new;
    if (work_range_old.first >= work_range_old.second) return -1;

    work_range_new.first += 1;

    bool success = false;
    while (!success) {
      work_range_new = Kokkos::atomic_compare_exchange(
          &m_work_range, work_range_old, work_range_new);
      success        = ((work_range_new == work_range_old) ||
                 (work_range_new.first >= work_range_new.second));
      work_range_old = work_range_new;
      work_range_new.first += 1;
    }
    if (work_range_old.first < work_range_old.second)
      return work_range_old.first;
    else
      return -1;
  }

  // Claim and index from this thread's range from the end
  inline long get_work_index_end() {
    Kokkos::pair<long, long> work_range_new = m_work_range;
    Kokkos::pair<long, long> work_range_old = work_range_new;
    if (work_range_old.first >= work_range_old.second) return -1;
    work_range_new.second -= 1;
    bool success = false;
    while (!success) {
      work_range_new = Kokkos::atomic_compare_exchange(
          &m_work_range, work_range_old, work_range_new);
      success        = ((work_range_new == work_range_old) ||
                 (work_range_new.first >= work_range_new.second));
      work_range_old = work_range_new;
      work_range_new.second -= 1;
    }
    if (work_range_old.first < work_range_old.second)
      return work_range_old.second - 1;
    else
      return -1;
  }

  // Reset the steal target
  inline void reset_steal_target() {
    m_current_steal_target = (m_pool_rank + 1) % pool_size();
    m_stealing             = false;
  }

  // Reset the steal target
  inline void reset_steal_target(int team_size) {
    m_current_steal_target = (m_pool_rank_rev + team_size);
    if (m_current_steal_target >= pool_size())
      m_current_steal_target = 0;  // pool_size()-1;
    m_stealing = false;
  }

  // Get a steal target; start with my-rank + 1 and go round robin, until
  // arriving at this threads rank Returns -1 fi no active steal target
  // available
  inline int get_steal_target() {
    while ((m_pool_base[m_current_steal_target]->m_work_range.second <=
            m_pool_base[m_current_steal_target]->m_work_range.first) &&
           (m_current_steal_target != m_pool_rank)) {
      m_current_steal_target = (m_current_steal_target + 1) % pool_size();
    }
    if (m_current_steal_target == m_pool_rank)
      return -1;
    else
      return m_current_steal_target;
  }

  inline int get_steal_target(int team_size) {
    while ((m_pool_base[m_current_steal_target]->m_work_range.second <=
            m_pool_base[m_current_steal_target]->m_work_range.first) &&
           (m_current_steal_target != m_pool_rank_rev)) {
      if (m_current_steal_target + team_size < pool_size())
        m_current_steal_target = (m_current_steal_target + team_size);
      else
        m_current_steal_target = 0;
    }

    if (m_current_steal_target == m_pool_rank_rev)
      return -1;
    else
      return m_current_steal_target;
  }

  inline long steal_work_index(int team_size = 0) {
    long index = -1;
    int steal_target =
        team_size > 0 ? get_steal_target(team_size) : get_steal_target();
    while ((steal_target != -1) && (index == -1)) {
      index = m_pool_base[steal_target]->get_work_index_end();
      if (index == -1)
        steal_target =
            team_size > 0 ? get_steal_target(team_size) : get_steal_target();
    }
    return index;
  }

  // Get a work index. Claim from owned range until its exhausted, then steal
  // from other thread
  inline long get_work_index(int team_size = 0) {
    long work_index = -1;
    if (!m_stealing) work_index = get_work_index_begin();

    if (work_index == -1) {
      memory_fence();
      m_stealing = true;
      work_index = steal_work_index(team_size);
    }

    m_team_work_index = work_index;
    memory_fence();
    return work_index;
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

inline int Threads::in_parallel() { return Impl::ThreadsExec::in_parallel(); }

inline int Threads::impl_is_initialized() {
  return Impl::ThreadsExec::is_initialized();
}

inline void Threads::impl_initialize(InitializationSettings const &settings) {
  Impl::ThreadsExec::initialize(
      settings.has_num_threads() ? settings.get_num_threads() : -1);
}

inline void Threads::impl_finalize() { Impl::ThreadsExec::finalize(); }

inline void Threads::print_configuration(std::ostream &os, bool verbose) const {
  os << "Host Parallel Execution Space:\n";
  os << "  KOKKOS_ENABLE_THREADS: yes\n";

  os << "\nThreads Runtime Configuration:\n";
  Impl::ThreadsExec::print_configuration(os, verbose);
}

inline void Threads::impl_static_fence(const std::string &name) {
  Impl::ThreadsExec::internal_fence(name, Impl::fence_is_static::yes);
}
} /* namespace Kokkos */

#endif /* #define KOKKOS_THREADSEXEC_HPP */
