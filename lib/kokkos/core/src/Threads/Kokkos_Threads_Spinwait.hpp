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

#ifndef KOKKOS_THREADS_SPINWAIT_HPP
#define KOKKOS_THREADS_SPINWAIT_HPP

#include <Threads/Kokkos_Threads_State.hpp>

#include <cstdint>

namespace Kokkos {
namespace Impl {

enum class WaitMode : int {
  ACTIVE  // Used for tight loops to keep threads active longest
  ,
  PASSIVE  // Used to quickly yield the thread to quite down the system
  ,
  ROOT  // Never sleep or yield the root thread
};

void host_thread_yield(const uint32_t i, const WaitMode mode);

void spinwait_while_equal(ThreadState const volatile& flag,
                          ThreadState const value);

}  // namespace Impl
}  // namespace Kokkos

#endif
