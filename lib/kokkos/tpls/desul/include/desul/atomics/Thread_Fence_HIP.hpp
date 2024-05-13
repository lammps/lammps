/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_THREAD_FENCE_HIP_HPP_
#define DESUL_ATOMICS_THREAD_FENCE_HIP_HPP_

#include <desul/atomics/Common.hpp>

namespace desul {
namespace Impl {

// clang-format off
inline __device__ void device_atomic_thread_fence(MemoryOrderRelease, MemoryScopeDevice) { __threadfence();        }
inline __device__ void device_atomic_thread_fence(MemoryOrderAcquire, MemoryScopeDevice) { __threadfence();        }
inline __device__ void device_atomic_thread_fence(MemoryOrderAcqRel , MemoryScopeDevice) { __threadfence();        }
inline __device__ void device_atomic_thread_fence(MemoryOrderSeqCst , MemoryScopeDevice) { __threadfence();        }
inline __device__ void device_atomic_thread_fence(MemoryOrderRelease, MemoryScopeCore  ) { __threadfence_block();  }
inline __device__ void device_atomic_thread_fence(MemoryOrderAcquire, MemoryScopeCore  ) { __threadfence_block();  }
inline __device__ void device_atomic_thread_fence(MemoryOrderAcqRel , MemoryScopeCore  ) { __threadfence_block();  }
inline __device__ void device_atomic_thread_fence(MemoryOrderSeqCst , MemoryScopeCore  ) { __threadfence_block();  }
inline __device__ void device_atomic_thread_fence(MemoryOrderRelease, MemoryScopeNode  ) { __threadfence_system(); }
inline __device__ void device_atomic_thread_fence(MemoryOrderAcquire, MemoryScopeNode  ) { __threadfence_system(); }
inline __device__ void device_atomic_thread_fence(MemoryOrderAcqRel , MemoryScopeNode  ) { __threadfence_system(); }
inline __device__ void device_atomic_thread_fence(MemoryOrderSeqCst , MemoryScopeNode  ) { __threadfence_system(); }
// clang-format on

}  // namespace Impl
}  // namespace desul

#endif
