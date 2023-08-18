/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_THREAD_FENCE_SCOPECALLER_HPP_
#define DESUL_ATOMICS_THREAD_FENCE_SCOPECALLER_HPP_

#include <desul/atomics/Common.hpp>

namespace desul {

// clang-format off
DESUL_INLINE_FUNCTION void atomic_thread_fence(MemoryOrderSeqCst , MemoryScopeCaller) {}
DESUL_INLINE_FUNCTION void atomic_thread_fence(MemoryOrderAcqRel , MemoryScopeCaller) {}
DESUL_INLINE_FUNCTION void atomic_thread_fence(MemoryOrderRelease, MemoryScopeCaller) {}
DESUL_INLINE_FUNCTION void atomic_thread_fence(MemoryOrderAcquire, MemoryScopeCaller) {}
DESUL_INLINE_FUNCTION void atomic_thread_fence(MemoryOrderRelaxed, MemoryScopeCaller) {}
// clang-format on

}  // namespace desul

#endif
