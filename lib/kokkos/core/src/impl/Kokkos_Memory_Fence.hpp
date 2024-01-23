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

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_MEMORY_FENCE_HPP)
#define KOKKOS_MEMORY_FENCE_HPP
namespace Kokkos {

//////////////////////////////////////////////////////
// store_fence()
//
// If possible use a store fence on the architecture, if not run a full memory
// fence

KOKKOS_FORCEINLINE_FUNCTION
void store_fence() {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
  asm volatile("sfence" ::: "memory");
#else
  memory_fence();
#endif
}

//////////////////////////////////////////////////////
// load_fence()
//
// If possible use a load fence on the architecture, if not run a full memory
// fence

KOKKOS_FORCEINLINE_FUNCTION
void load_fence() {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
  asm volatile("lfence" ::: "memory");
#else
  memory_fence();
#endif
}

}  // namespace Kokkos

#endif
