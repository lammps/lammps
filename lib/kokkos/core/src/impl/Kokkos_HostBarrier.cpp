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

#include <impl/Kokkos_HostBarrier.hpp>
#include <impl/Kokkos_BitOps.hpp>

#include <impl/Kokkos_HostBarrier.hpp>

#include <thread>
#if defined(_WIN32)
#include <process.h>
#include <winsock2.h>
#include <windows.h>
#endif

namespace Kokkos {
namespace Impl {

void HostBarrier::impl_backoff_wait_until_equal(
    int* ptr, const int v, const bool active_wait) noexcept {
  unsigned count = 0u;

  while (!test_equal(ptr, v)) {
    const int c = int_log2(++count);
    if (!active_wait || c > log2_iterations_till_sleep) {
      std::this_thread::sleep_for(
          std::chrono::nanoseconds(c < 16 ? 256 * c : 4096));
    } else if (c > log2_iterations_till_yield) {
      std::this_thread::yield();
    }
#if defined(KOKKOS_ENABLE_ASM)
#if defined(__PPC64__)
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
  }
}
}  // namespace Impl
}  // namespace Kokkos
