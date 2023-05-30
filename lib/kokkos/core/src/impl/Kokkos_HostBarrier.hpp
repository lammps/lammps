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

#ifndef KOKKOS_HOST_BARRIER_HPP
#define KOKKOS_HOST_BARRIER_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Atomic.hpp>

namespace Kokkos {
namespace Impl {

// class HostBarrier
//
// provides a static and member interface for a barrier shared between threads
// of execution.
//
// *buffer* is a shared resource between the threads of execution
// *step* should be a stack variable associated with the current thread of
// execution *size* is the number of threads which share the barrier
//
// before calling any arrive type function the buffer and step must have been
// initialized to 0 and one of the following conditions must be true
//
// 1) step == 0 (i.e. first arrive call to HostBarrier),
// 2) try_wait has returned true for the current thread of execution,
// 3) a wait type function has returned for the current thread of execution, or
// 4) split_arrive returned true on the current thread of execution and it has
//    called split_release
//
// The purporse of the split functions is to allow the last thread to arrive
// an opportunity to perform some actions before releasing the waiting threads
//
// If all threads have arrived (and split_release has been call if using
// split_arrive) before a wait type call, the wait may return quickly
class HostBarrier {
 public:
  using buffer_type                         = int;
  static constexpr int required_buffer_size = 128;
  static constexpr int required_buffer_length =
      required_buffer_size / sizeof(int);

 private:
  // fit the following 3 atomics within a 128 bytes while
  // keeping the arrive atomic at least 64 bytes away from
  // the wait atomic to reduce contention on the caches
  static constexpr int arrive_idx = 32 / sizeof(int);
  static constexpr int master_idx = 64 / sizeof(int);
  static constexpr int wait_idx   = 96 / sizeof(int);

  static constexpr int num_nops                   = 32;
  static constexpr int iterations_till_backoff    = 64;
  static constexpr int log2_iterations_till_yield = 4;
  static constexpr int log2_iterations_till_sleep = 6;

 public:
  // will return true if call is the last thread to arrive
  KOKKOS_INLINE_FUNCTION
  static bool split_arrive(int* buffer, const int size, int& step,
                           const bool master_wait = true) noexcept {
    if (size <= 1) return true;

    ++step;
    Kokkos::memory_fence();
    const bool result =
        Kokkos::atomic_fetch_add(buffer + arrive_idx, 1) == size - 1;

    if (master_wait && result) {
      Kokkos::atomic_fetch_add(buffer + master_idx, 1);
    }

    return result;
  }

  // release waiting threads
  // only the thread which received a return value of true from split_arrive
  // or the thread which calls split_master_wait may call split_release
  KOKKOS_INLINE_FUNCTION
  static void split_release(int* buffer, const int size, const int /*step*/
                            ) noexcept {
    if (size <= 1) return;
    Kokkos::memory_fence();
    Kokkos::atomic_fetch_sub(buffer + arrive_idx, size);
    Kokkos::atomic_fetch_add(buffer + wait_idx, 1);
  }

  // should only be called by the master thread, will allow the master thread to
  // resume after all threads have arrived
  KOKKOS_INLINE_FUNCTION
  static void split_master_wait(int* buffer, const int size, const int step,
                                const bool active_wait = true) noexcept {
    if (size <= 1) return;
    wait_until_equal(buffer + master_idx, step, active_wait);
  }

  // arrive, last thread automatically release waiting threads
  KOKKOS_INLINE_FUNCTION
  static void arrive(int* buffer, const int size, int& step) noexcept {
    if (size <= 1) return;
    if (split_arrive(buffer, size, step)) {
      split_release(buffer, size, step);
    }
  }

  // test if all threads have arrived
  KOKKOS_INLINE_FUNCTION
  static bool try_wait(int* buffer, const int size, const int step) noexcept {
    if (size <= 1) return true;
    return test_equal(buffer + wait_idx, step);
  }

  // wait for all threads to arrive
  KOKKOS_INLINE_FUNCTION
  static void wait(int* buffer, const int size, const int step,
                   bool active_wait = true) noexcept {
    if (size <= 1) return;
    wait_until_equal(buffer + wait_idx, step, active_wait);
  }

 public:
  KOKKOS_INLINE_FUNCTION
  bool split_arrive(const bool master_wait = true) const noexcept {
    return split_arrive(m_buffer, m_size, m_step, master_wait);
  }

  KOKKOS_INLINE_FUNCTION
  void split_release() const noexcept {
    split_release(m_buffer, m_size, m_step);
  }

  KOKKOS_INLINE_FUNCTION
  void split_master_wait(const bool active_wait = true) noexcept {
    split_master_wait(m_buffer, m_size, m_step, active_wait);
  }

  KOKKOS_INLINE_FUNCTION
  void arrive() const noexcept { return arrive(m_buffer, m_size, m_step); }

  KOKKOS_INLINE_FUNCTION
  bool try_wait() const noexcept { return try_wait(m_buffer, m_size, m_step); }

  KOKKOS_INLINE_FUNCTION
  void wait() const noexcept { wait(m_buffer, m_size, m_step); }

  HostBarrier()              = default;
  HostBarrier(HostBarrier&&) = default;
  HostBarrier& operator=(HostBarrier&&) = default;

  KOKKOS_INLINE_FUNCTION
  HostBarrier(int size, int* buffer)
      : m_size{size}, m_step{0u}, m_buffer{buffer} {}

  HostBarrier(const HostBarrier&) = delete;
  HostBarrier& operator=(const HostBarrier&) = delete;

 private:
  KOKKOS_INLINE_FUNCTION
  static bool test_equal(int* ptr, int v) noexcept {
    const bool result = Kokkos::atomic_fetch_add(ptr, 0) == v;
    if (result) {
      Kokkos::memory_fence();
    }
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  static void wait_until_equal(int* ptr, const int v,
                               bool active_wait = true) noexcept {
    KOKKOS_IF_ON_HOST((impl_wait_until_equal_host(ptr, v, active_wait);))

    KOKKOS_IF_ON_DEVICE(((void)active_wait; while (!test_equal(ptr, v)){}))
  }

  static void impl_wait_until_equal_host(int* ptr, const int v,
                                         bool active_wait = true) noexcept {
    bool result = test_equal(ptr, v);
    for (int i = 0; !result && i < iterations_till_backoff; ++i) {
#if defined(KOKKOS_ENABLE_ASM)
#if defined(_WIN32)
      for (int j = 0; j < num_nops; ++j) {
        __asm__ __volatile__("nop\n");
      }
      __asm__ __volatile__("pause\n" ::: "memory");
#elif defined(__PPC64__)
      for (int j = 0; j < num_nops; ++j) {
        asm volatile("nop\n");
      }
      asm volatile("or 27, 27, 27" ::: "memory");
#elif defined(__amd64) || defined(__amd64__) || defined(__x86_64) || \
    defined(__x86_64__)
      for (int j = 0; j < num_nops; ++j) {
        asm volatile("nop\n");
      }
      asm volatile("pause\n" ::: "memory");
#endif
#endif
      result = test_equal(ptr, v);
    }
    if (!result) {
      impl_backoff_wait_until_equal(ptr, v, active_wait);
    }
  }

  static void impl_backoff_wait_until_equal(int* ptr, const int v,
                                            const bool active_wait) noexcept;

 private:
  int m_size{0};
  mutable int m_step{0};
  int* m_buffer{nullptr};
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_HOST_BARRIER_HPP
