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
#endif

#include <Kokkos_Macros.hpp>

#include <Kokkos_Atomic.hpp>
#include <Threads/Kokkos_Threads_Spinwait.hpp>
#include <impl/Kokkos_BitOps.hpp>

#include <thread>
#if defined(_WIN32)
#include <process.h>
#include <winsock2.h>
#include <windows.h>
#endif

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

void host_thread_yield(const uint32_t i, const WaitMode mode) {
  static constexpr uint32_t sleep_limit = 1 << 13;
  static constexpr uint32_t yield_limit = 1 << 12;

  const int c = int_log2(i);

  if (WaitMode::ROOT != mode) {
    if (sleep_limit < i) {
      // Attempt to put the thread to sleep for 'c' microseconds
      std::this_thread::yield();
      std::this_thread::sleep_for(std::chrono::microseconds(c));
    }

    else if (mode == WaitMode::PASSIVE || yield_limit < i) {
      // Attempt to yield thread resources to runtime
      std::this_thread::yield();
    }
#if defined(KOKKOS_ENABLE_ASM)

    else if ((1u << 4) < i) {

      // Insert a few no-ops to quiet the thread:

      for (int k = 0; k < c; ++k) {
#if defined(__amd64) || defined(__amd64__) || defined(__x86_64) || \
    defined(__x86_64__)
#if !defined(_WIN32) /* IS NOT Microsoft Windows */
        asm volatile("nop\n");
#else
        __asm__ __volatile__("nop\n");
#endif
#elif defined(__PPC64__)
        asm volatile("nop\n");
#endif
      }
    }
#endif /* defined( KOKKOS_ENABLE_ASM ) */
  }
#if defined(KOKKOS_ENABLE_ASM)
  else if ((1u << 3) < i) {
    // no-ops for root thread
    for (int k = 0; k < c; ++k) {
#if defined(__amd64) || defined(__amd64__) || defined(__x86_64) || \
    defined(__x86_64__)
#if !defined(_WIN32) /* IS NOT Microsoft Windows */
      asm volatile("nop\n");
#else
      __asm__ __volatile__("nop\n");
#endif
#elif defined(__PPC64__)
      asm volatile("nop\n");
#endif
    }
  }

  {
    // Insert memory pause
#if defined(__amd64) || defined(__amd64__) || defined(__x86_64) || \
    defined(__x86_64__)
#if !defined(_WIN32) /* IS NOT Microsoft Windows */
    asm volatile("pause\n" ::: "memory");
#else
    __asm__ __volatile__("pause\n" ::: "memory");
#endif
#elif defined(__PPC64__)
    asm volatile("or 27, 27, 27" ::: "memory");
#endif
  }

#endif /* defined( KOKKOS_ENABLE_ASM ) */
}

void spinwait_while_equal(ThreadState const volatile& flag,
                          ThreadState const value) {
  Kokkos::store_fence();
  uint32_t i = 0;
  while (value == flag) {
    host_thread_yield(++i, WaitMode::ACTIVE);
  }
  Kokkos::load_fence();
}

}  // namespace Impl
}  // namespace Kokkos
