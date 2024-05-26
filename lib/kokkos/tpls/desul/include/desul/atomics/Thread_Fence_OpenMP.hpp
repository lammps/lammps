/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_THREAD_FENCE_OPENMP_HPP_
#define DESUL_ATOMICS_THREAD_FENCE_OPENMP_HPP_

#include <omp.h>

#include <desul/atomics/Common.hpp>

namespace desul {
namespace Impl {

// NVHPC compiler only supports the basic flush construct without the
// memory-order-clause.
#if _OPENMP > 201800 && !defined(__NVCOMPILER)

// There is no seq_cst flush in OpenMP, isn't it the same anyway for fence?
inline void host_atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeCore) {
#pragma omp flush acq_rel
}
inline void host_atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeCore) {
#pragma omp flush acq_rel
}
inline void host_atomic_thread_fence(MemoryOrderRelease, MemoryScopeCore) {
#pragma omp flush release
}
inline void host_atomic_thread_fence(MemoryOrderAcquire, MemoryScopeCore) {
#pragma omp flush acquire
}
inline void host_atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeDevice) {
#pragma omp flush acq_rel
}
inline void host_atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeDevice) {
#pragma omp flush acq_rel
}
inline void host_atomic_thread_fence(MemoryOrderRelease, MemoryScopeDevice) {
#pragma omp flush release
}
inline void host_atomic_thread_fence(MemoryOrderAcquire, MemoryScopeDevice) {
#pragma omp flush acquire
}

#else

inline void host_atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeCore) {
#pragma omp flush
}
inline void host_atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeCore) {
#pragma omp flush
}
inline void host_atomic_thread_fence(MemoryOrderRelease, MemoryScopeCore) {
#pragma omp flush
}
inline void host_atomic_thread_fence(MemoryOrderAcquire, MemoryScopeCore) {
#pragma omp flush
}
inline void host_atomic_thread_fence(MemoryOrderSeqCst, MemoryScopeDevice) {
#pragma omp flush
}
inline void host_atomic_thread_fence(MemoryOrderAcqRel, MemoryScopeDevice) {
#pragma omp flush
}
inline void host_atomic_thread_fence(MemoryOrderRelease, MemoryScopeDevice) {
#pragma omp flush
}
inline void host_atomic_thread_fence(MemoryOrderAcquire, MemoryScopeDevice) {
#pragma omp flush
}

#endif

}  // namespace Impl
}  // namespace desul

#endif
